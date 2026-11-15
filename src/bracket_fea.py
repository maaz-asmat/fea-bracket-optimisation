"""
2D plane-stress finite element model of a steel L-bracket.

Units: mm, N, MPa (N/mm^2). Thickness t is out of plane.
Open-source stack: Triangle (meshing) + scikit-fem (FE solver), quadratic (P2) triangles.
"""
from dataclasses import dataclass, field
import time

import numpy as np
import triangle
from skfem import (Basis, ElementTriP2, ElementVector, Functional, LinearForm, MeshTri,
                   FacetBasis, asm, condense, solve)
from skfem.helpers import dot, sym_grad
from skfem.models.elasticity import lame_parameters, linear_elasticity

# ---------------------------------------------------------------- parameters

@dataclass
class Bracket:
    """L-bracket: vertical leg bolted to a wall (x = 0), horizontal arm loaded at its free end."""
    leg_width: float = 20.0      # mm, vertical leg thickness in x
    height: float = 100.0        # mm, total height
    arm_length: float = 120.0    # mm, total width (wall to load)
    arm_depth: float = 20.0      # mm, depth of horizontal arm
    fillet: float = 5.0          # mm, radius at the inner (re-entrant) corner; 0 = sharp
    thickness: float = 10.0      # mm, out-of-plane plate thickness

    @property
    def corner(self):
        return np.array([self.leg_width, self.height - self.arm_depth])


@dataclass
class Material:
    name: str = "Structural steel S235"
    E: float = 210_000.0         # MPa
    nu: float = 0.30
    yield_strength: float = 235.0  # MPa
    density: float = 7.85e-6     # kg/mm^3


@dataclass
class MeshSettings:
    h_far: float = 4.0           # mm, element size away from the corner
    h_fine: float = 0.5          # mm, element size at the fillet
    refine_radius: float = 15.0  # mm, zone around the corner that gets refined
    uniform: bool = False        # True = h_fine everywhere (for comparison)


# ---------------------------------------------------------------- geometry

def _points_along(p0, p1, h):
    n = max(1, int(np.ceil(np.linalg.norm(p1 - p0) / h)))
    t = np.linspace(0, 1, n + 1)[:-1]
    return p0 + np.outer(t, p1 - p0)


def outline(b: Bracket, ms: MeshSettings):
    """Closed boundary polygon (counter-clockwise) with local point spacing."""
    def h_at(p):
        if ms.uniform:
            return ms.h_fine
        d = np.linalg.norm(p - b.corner)
        return ms.h_fine if d < ms.refine_radius else ms.h_far

    cx, cy = b.corner
    r = b.fillet
    key = [np.array(p, float) for p in [(0, 0), (b.leg_width, 0), (cx, cy - r)]]
    pts = []
    # straight edges up to the fillet start, with finer spacing near the corner
    for p0, p1 in zip(key[:-1], key[1:]):
        seg = _points_along(p0, p1, h_at(p1) if np.linalg.norm(p1 - b.corner) < ms.refine_radius + 1 else ms.h_far)
        pts.extend(seg)
    # fillet arc: centre (cx + r, cy - r), from 180 deg to 90 deg
    if r > 0:
        n_arc = max(4, int(np.ceil(0.5 * np.pi * r / ms.h_fine)))
        ang = np.linspace(np.pi, np.pi / 2, n_arc + 1)[:-1]
        pts.extend(np.c_[cx + r + r * np.cos(ang), cy - r + r * np.sin(ang)])
    rest = [np.array(p, float) for p in [(cx + r, cy), (b.arm_length, cy), (b.arm_length, b.height), (0, b.height), (0, 0)]]
    for p0, p1 in zip(rest[:-1], rest[1:]):
        h = min(h_at(p0), h_at(p1))
        if not ms.uniform and np.linalg.norm(p0 - b.corner) < ms.refine_radius + r + 1:
            # grade from fine to coarse along the arm's lower edge
            d = np.linalg.norm(p1 - p0)
            s = [0.0]
            while s[-1] < d:
                x = p0 + (p1 - p0) * s[-1] / d
                s.append(s[-1] + h_at(x))
            s = np.array(s[:-1]) / d
            pts.extend(p0 + np.outer(s, p1 - p0))
        else:
            pts.extend(_points_along(p0, p1, h))
    pts = np.array(pts)
    seg = np.c_[np.arange(len(pts)), np.roll(np.arange(len(pts)), -1)]
    return pts, seg


def build_mesh(b: Bracket, ms: MeshSettings) -> MeshTri:
    pts, seg = outline(b, ms)
    area_far = np.sqrt(3) / 4 * ms.h_far ** 2
    area_fine = np.sqrt(3) / 4 * ms.h_fine ** 2
    t = triangle.triangulate({"vertices": pts, "segments": seg}, f"pq30a{(area_fine if ms.uniform else area_far):.6f}")
    if not ms.uniform:
        for _ in range(6):  # graded local refinement around the corner
            c = t["vertices"][t["triangles"]].mean(axis=1)
            d = np.linalg.norm(c - b.corner, axis=1)
            target = np.where(d < ms.refine_radius, area_fine,
                              np.minimum(area_far, area_fine * (1 + (d - ms.refine_radius) / ms.h_fine) ** 2))
            t["triangle_max_area"] = target
            t = triangle.triangulate(t, "rpq30a")
    return MeshTri(t["vertices"].T, t["triangles"].T)


def straight_cantilever_mesh(length, depth, h):
    pts = np.array([(0, 0), (length, 0), (length, depth), (0, depth)], float)
    seg = np.array([(0, 1), (1, 2), (2, 3), (3, 0)])
    t = triangle.triangulate({"vertices": pts, "segments": seg}, f"pq30a{np.sqrt(3) / 4 * h ** 2:.6f}")
    return MeshTri(t["vertices"].T, t["triangles"].T)


# ---------------------------------------------------------------- solver

@dataclass
class Result:
    n_elements: int
    n_dofs: int
    runtime_s: float
    tip_deflection: float            # mm, downward positive
    peak_von_mises: float            # MPa, in the region of interest
    peak_location: tuple
    mass_kg: float
    vm_qp: np.ndarray = field(repr=False, default=None)
    basis: object = field(repr=False, default=None)
    u: np.ndarray = field(repr=False, default=None)


def solve_plane_stress(mesh: MeshTri, mat: Material, thickness: float, load_N: float,
                       x_load: float, roi_centre=None, roi_radius=None) -> Result:
    """Fixed edge at x = 0, total load `load_N` (downward) spread over the edge x = x_load."""
    t0 = time.perf_counter()
    e = ElementVector(ElementTriP2())
    ib = Basis(mesh, e, intorder=4)

    lam, mu = lame_parameters(mat.E, mat.nu)
    lam_ps = 2 * lam * mu / (lam + 2 * mu)          # plane-stress modification
    K = asm(linear_elasticity(lam_ps, mu), ib)

    load_facets = mesh.facets_satisfying(lambda x: np.isclose(x[0], x_load))
    fb = FacetBasis(mesh, e, facets=load_facets, intorder=4)
    y = mesh.p[1, mesh.facets[:, load_facets]]
    loaded_height = y.max() - y.min()
    q = load_N / (thickness * loaded_height)     # MPa, uniform shear traction (downward)

    @LinearForm
    def traction(v, w):
        return -q * v[1]

    f = asm(traction, fb)
    fixed = ib.get_dofs(lambda x: np.isclose(x[0], 0.0)).all()
    u = solve(*condense(K, f, D=fixed))

    # stresses at quadrature points
    uh = ib.interpolate(u)
    eps = sym_grad(uh)
    tr = eps[0, 0] + eps[1, 1]
    sxx = 2 * mu * eps[0, 0] + lam_ps * tr
    syy = 2 * mu * eps[1, 1] + lam_ps * tr
    sxy = 2 * mu * eps[0, 1]
    vm = np.sqrt(sxx ** 2 - sxx * syy + syy ** 2 + 3 * sxy ** 2)
    xq = ib.global_coordinates().value  # (2, nelem, nqp)

    mask = np.ones_like(vm, bool)
    if roi_centre is not None:
        d = np.sqrt((xq[0] - roi_centre[0]) ** 2 + (xq[1] - roi_centre[1]) ** 2)
        mask = d < roi_radius
    i = np.unravel_index(np.argmax(np.where(mask, vm, -1)), vm.shape)
    peak = float(vm[i])
    loc = (float(xq[0][i]), float(xq[1][i]))

    # tip deflection: mean vertical displacement of the loaded edge
    tip_dofs = ib.get_dofs(lambda x: np.isclose(x[0], x_load)).nodal["u^2"]
    tip = float(-u[tip_dofs].mean())

    area = float(Functional(lambda w: 1.0 + 0 * w.x[0]).assemble(Basis(mesh, ElementTriP2())))  # mm^2
    mass = area * thickness * mat.density
    runtime = time.perf_counter() - t0
    return Result(mesh.t.shape[1], ib.N, runtime, tip, peak, loc, mass, vm, ib, u)
