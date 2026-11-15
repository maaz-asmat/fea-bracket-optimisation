"""Step 1 - Validate the solver on a straight cantilever against beam theory.

A 120 x 20 mm steel bar (10 mm thick) is fixed at x = 0 and loaded with 600 N at x = 120.
Compare the FE tip deflection and the bending stress at mid-span with
Euler-Bernoulli + Timoshenko (shear) beam theory.
"""
import sys, pathlib, csv
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
from bracket_fea import Material, straight_cantilever_mesh, solve_plane_stress

OUT = pathlib.Path(__file__).resolve().parents[1] / "results"
OUT.mkdir(exist_ok=True)

L, h, t, P = 120.0, 20.0, 10.0, 600.0
mat = Material()
I = t * h ** 3 / 12
G = mat.E / (2 * (1 + mat.nu))
kappa = 5 / 6
delta_bending = P * L ** 3 / (3 * mat.E * I)
delta_shear = P * L / (kappa * G * t * h)
delta_theory = delta_bending + delta_shear
x_mid = L / 2
sigma_theory = P * (L - x_mid) * (h / 2) / I

rows = []
for h_el in [8, 4, 2, 1]:
    mesh = straight_cantilever_mesh(L, h, h_el)
    r = solve_plane_stress(mesh, mat, t, P, x_load=L)
    # bending stress sigma_xx on the top fibre at mid-span
    xq = r.basis.global_coordinates().value
    uh = r.basis.interpolate(r.u)
    from skfem.helpers import sym_grad
    from skfem.models.elasticity import lame_parameters
    lam, mu = lame_parameters(mat.E, mat.nu); lam_ps = 2 * lam * mu / (lam + 2 * mu)
    eps = sym_grad(uh); tr = eps[0, 0] + eps[1, 1]
    sxx = 2 * mu * eps[0, 0] + lam_ps * tr
    near = (np.abs(xq[0] - x_mid) < h_el) & (xq[1] > h - h_el)
    # extrapolate linearly in y over the elements near the top fibre at mid-span
    yy, ss = xq[1][near], sxx[near]
    sigma_fe = float(np.polyval(np.polyfit(yy, ss, 1), h))
    rows.append([h_el, r.n_elements, r.n_dofs, round(r.tip_deflection, 4), round(sigma_fe, 2),
                 round(100 * (r.tip_deflection / delta_theory - 1), 2), round(100 * (sigma_fe / sigma_theory - 1), 2),
                 round(r.runtime_s, 3)])

hdr = ["element_size_mm", "elements", "dofs", "tip_deflection_mm", "sigma_xx_midspan_MPa",
       "deflection_error_%", "stress_error_%", "runtime_s"]
with open(OUT / "01_validation_cantilever.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(hdr); w.writerows(rows)

print(f"Beam theory: tip deflection = {delta_theory:.4f} mm "
      f"(bending {delta_bending:.4f} + shear {delta_shear:.4f}), sigma_xx(mid-span, top) = {sigma_theory:.2f} MPa")
for r in rows:
    print(dict(zip(hdr, r)))
