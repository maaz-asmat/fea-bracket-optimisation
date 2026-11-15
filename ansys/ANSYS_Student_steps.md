# Repeating the study in ANSYS Mechanical (Student)

ANSYS Student is free and runs on Windows. The 2D model in this guide is small, so it runs on an ordinary laptop.

## 1. Geometry (SpaceClaim or DesignModeler)
1. New Static Structural system. In **Geometry → Properties**, set *Analysis Type* to **2D**.
2. In the XY plane, sketch the L-shape: outer corners (0,0), (20,0), (20,80), (120,80), (120,100), (0,100), all in mm.
3. Add a **fillet** at the inner corner (20,80). Start with r = 5 mm.
4. Create a surface from the sketch.

## 2. Material and section
1. Engineering Data: **Structural Steel**, or define S235 (E = 210 GPa, ν = 0.3, yield 235 MPa).
2. In Mechanical, select the geometry. Set *2D Behavior* to **Plane Stress** and *Thickness* to **10 mm**.

## 3. Mesh
1. Global element size: **4 mm**, element order **Quadratic**.
2. Add **Edge Sizing** on the fillet edge: **0.5 mm**. Optionally add a Sphere of Influence of radius 15 mm around the corner.
3. Run a convergence check. Solve with fillet sizes of 2, 1, 0.5 and 0.25 mm and note the peak stress and solve time.

## 4. Supports and load
1. **Fixed Support** on the left edge (x = 0).
2. **Force** on the right edge (x = 120): Components, Y = −600 N.

## 5. Results
1. Insert **Equivalent (von Mises) Stress**, **Directional Deformation (Y)** and a **Stress Tool** (safety factor, max equivalent stress, yield).
2. Probe the peak stress at the fillet. Ignore the small singular spots at the corners of the fixed edge.

## 6. Design study
1. Make the fillet radius and thickness **parameters** (tick the P box).
2. Use the **Parameter Set** table to run r = 2, 5, 10, 15 mm, and t = 8 mm with r = 10 mm.
3. Compare your results with `results/03_design_study.csv`. Differences of a few percent are normal between solvers and meshes.

## Expected results (open-source model)
| Variant | Peak von Mises | Safety factor |
|---|---|---|
| r = 5 mm | ~149 MPa | ~1.58 |
| r = 10 mm | ~117 MPa | ~2.0 |
| r = 10 mm, t = 8 mm | ~146 MPa | ~1.6 |
