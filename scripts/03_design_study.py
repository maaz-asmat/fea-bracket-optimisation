"""Step 3 - Design-variant study: fillet radius and plate thickness.

Each variant uses the mesh settings found in step 2 (0.5 mm at the fillet, or
r/10 for small fillets). Outputs a results table, a comparison chart and
von Mises stress plots of the baseline and the recommended design.
"""
import sys, pathlib, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
from bracket_fea import Bracket, Material, MeshSettings, build_mesh, solve_plane_stress

OUT = pathlib.Path(__file__).resolve().parents[1] / "results"
LOAD = 600.0          # N at the free end of the arm
mat = Material()

variants = [
    ("A  r = 2 mm",            Bracket(fillet=2.0)),
    ("B  r = 5 mm (baseline)", Bracket(fillet=5.0)),
    ("C  r = 10 mm",           Bracket(fillet=10.0)),
    ("D  r = 15 mm",           Bracket(fillet=15.0)),
    ("E  r = 10 mm, t = 8 mm", Bracket(fillet=10.0, thickness=8.0)),
]

results = []
for name, b in variants:
    ms = MeshSettings(h_far=4.0, h_fine=min(0.5, b.fillet / 10), refine_radius=max(15.0, 1.5 * b.fillet))
    mesh = build_mesh(b, ms)
    r = solve_plane_stress(mesh, mat, b.thickness, LOAD, x_load=b.arm_length,
                           roi_centre=b.corner + [b.fillet * 0.3, -b.fillet * 0.3],
                           roi_radius=max(12.0, 1.5 * b.fillet))
    results.append((name, b, mesh, r))

base = results[1][3]
rows = []
for name, b, mesh, r in results:
    fos = mat.yield_strength / r.peak_von_mises
    rows.append([name, b.fillet, b.thickness, r.n_elements, round(r.peak_von_mises, 1),
                 round(100 * (r.peak_von_mises / base.peak_von_mises - 1), 1),
                 round(r.tip_deflection, 4), round(100 * (r.tip_deflection / base.tip_deflection - 1), 1),
                 round(r.mass_kg, 4), round(100 * (r.mass_kg / base.mass_kg - 1), 1), round(fos, 2)])
hdr = ["variant", "fillet_mm", "thickness_mm", "elements", "peak_von_mises_MPa", "stress_vs_baseline_%",
       "tip_deflection_mm", "deflection_vs_baseline_%", "mass_kg", "mass_vs_baseline_%", "safety_factor_vs_yield"]
with open(OUT / "03_design_study.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(hdr); w.writerows(rows)
for r in rows:
    print(dict(zip(hdr, r)))

# ---- comparison chart
names = [r[0] for r in rows]
fig, ax = plt.subplots(1, 3, figsize=(13, 4))
ax[0].barh(names, [r[4] for r in rows], color="#1f5f8b"); ax[0].set_xlabel("Peak von Mises [MPa]")
ax[0].axvline(mat.yield_strength, color="crimson", ls="--", lw=1); ax[0].text(mat.yield_strength, -0.6, " yield", color="crimson", fontsize=8)
ax[1].barh(names, [r[10] for r in rows], color="#0e7c7b"); ax[1].set_xlabel("Safety factor against yield")
ax[1].axvline(1.5, color="0.4", ls=":", lw=1); ax[1].text(1.5, -0.6, " target 1.5", color="0.3", fontsize=8)
ax[2].barh(names, [r[8] * 1000 for r in rows], color="#d97706"); ax[2].set_xlabel("Mass [g]")
for a in ax:
    a.invert_yaxis(); a.grid(axis="x", alpha=.3)
for a in ax[1:]:
    a.set_yticklabels([])
fig.suptitle(f"L-bracket design variants, S235 steel, {LOAD:.0f} N tip load")
fig.tight_layout()
fig.savefig(OUT / "03_design_comparison.png", dpi=150)

# ---- von Mises fields
def field_plot(idx, fname):
    name, b, mesh, r = results[idx]
    vm_el = r.vm_qp.mean(axis=1)
    tri = mtri.Triangulation(mesh.p[0], mesh.p[1], mesh.t.T)
    fig, ax = plt.subplots(figsize=(7, 5.4))
    tpc = ax.tripcolor(tri, facecolors=vm_el, cmap="turbo", vmin=0, vmax=160)
    fig.colorbar(tpc, ax=ax, label="von Mises stress [MPa]")
    ax.plot(*r.peak_location, "wx", ms=10, mew=2)
    ax.annotate(f"peak {r.peak_von_mises:.0f} MPa", r.peak_location, (r.peak_location[0] + 20, r.peak_location[1] - 25),
                color="k", arrowprops=dict(arrowstyle="->"))
    ax.set_aspect("equal"); ax.set(xlabel="x [mm]", ylabel="y [mm]", title=f"Variant {name}")
    ax.annotate("", (b.arm_length - 3, b.height + 1), (b.arm_length - 3, b.height + 16),
                arrowprops=dict(arrowstyle="->", color="crimson", lw=2))
    ax.text(b.arm_length - 6, b.height + 9, f"{LOAD:.0f} N", color="crimson", ha="right")
    ax.set_ylim(-5, b.height + 20)
    ax.plot([0, 0], [0, b.height], color="k", lw=4); ax.text(-3, 5, "fixed to wall", rotation=90, fontsize=8, ha="right")
    fig.tight_layout(); fig.savefig(OUT / fname, dpi=150)

field_plot(1, "03_von_mises_baseline_r5.png")
field_plot(2, "03_von_mises_r10.png")
