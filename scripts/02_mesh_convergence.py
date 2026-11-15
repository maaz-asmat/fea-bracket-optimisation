"""Step 2 - Mesh convergence at the fillet, and local vs uniform refinement.

Baseline bracket (5 mm fillet). The peak von Mises stress is read in a 12 mm zone
around the inner corner. We refine (a) only near the fillet and (b) everywhere,
and compare accuracy, model size and runtime. A sharp corner (r = 0) is included
to show why a zero-radius corner never converges (stress singularity).
"""
import sys, pathlib, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))
from bracket_fea import Bracket, Material, MeshSettings, build_mesh, solve_plane_stress

OUT = pathlib.Path(__file__).resolve().parents[1] / "results"
LOAD = 600.0
mat = Material()
ROI = 12.0


def run(b, ms):
    mesh = build_mesh(b, ms)
    best = None
    for _ in range(3):   # repeat and keep fastest run for a fair timing
        r = solve_plane_stress(mesh, mat, b.thickness, LOAD, x_load=b.arm_length,
                               roi_centre=b.corner + [b.fillet * 0.3, -b.fillet * 0.3], roi_radius=ROI)
        best = r if best is None or r.runtime_s < best.runtime_s else best
    return best


rows = []
base = Bracket(fillet=5.0)
for h in [4, 2, 1, 0.5, 0.25]:
    r = run(base, MeshSettings(h_far=4.0, h_fine=h))
    rows.append(["local", 5.0, h, r.n_elements, int(r.n_dofs), round(r.peak_von_mises, 2), round(r.runtime_s, 3)])
for h in [4, 2, 1, 0.5]:
    r = run(base, MeshSettings(h_far=h, h_fine=h, uniform=True))
    rows.append(["uniform", 5.0, h, r.n_elements, int(r.n_dofs), round(r.peak_von_mises, 2), round(r.runtime_s, 3)])
sharp = Bracket(fillet=0.0)
for h in [2, 1, 0.5, 0.25, 0.125]:
    r = run(sharp, MeshSettings(h_far=4.0, h_fine=h))
    rows.append(["local", 0.0, h, r.n_elements, int(r.n_dofs), round(r.peak_von_mises, 2), round(r.runtime_s, 3)])

hdr = ["refinement", "fillet_mm", "element_size_at_corner_mm", "elements", "dofs", "peak_von_mises_MPa", "runtime_s"]
with open(OUT / "02_mesh_convergence.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(hdr); w.writerows(rows)
for r in rows:
    print(dict(zip(hdr, r)))

fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
for kind, style, label in [("local", "o-", "local refinement at the fillet"), ("uniform", "s--", "uniform mesh")]:
    d = [r for r in rows if r[0] == kind and r[1] == 5.0]
    ax[0].plot([r[4] for r in d], [r[5] for r in d], style, label=label)
    ax[1].plot([r[6] for r in d], [r[5] for r in d], style, label=label)
d = [r for r in rows if r[1] == 0.0]
ax[2].plot([r[2] for r in d], [r[5] for r in d], "^:", color="crimson")
ax[0].set(xscale="log", xlabel="Degrees of freedom", ylabel="Peak von Mises [MPa]", title="(a) 5 mm fillet: convergence")
ax[1].set(xscale="log", xlabel="Solve time [s]", ylabel="Peak von Mises [MPa]", title="(b) 5 mm fillet: accuracy vs runtime")
ax[2].set(xscale="log", xlabel="Element size at corner [mm]", ylabel="Peak von Mises [MPa]", title="(c) Sharp corner never converges")
ax[2].invert_xaxis()
ax[2].set_xticks([2, 1, 0.5, 0.25, 0.125]); ax[2].set_xticklabels(["2", "1", "0.5", "0.25", "0.125"]); ax[2].minorticks_off()
for a in ax:
    a.grid(alpha=.3)
ax[0].legend(fontsize=8); ax[1].legend(fontsize=8)
fig.tight_layout()
fig.savefig(OUT / "02_mesh_convergence.png", dpi=150)
