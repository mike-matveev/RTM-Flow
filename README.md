# MATCVFEM 3D – RTM filling on shell geometries

A MATLAB code for simulating resin flow in resin transfer moulding (RTM) with
the control volume finite element method (CV-FEM). The part is modelled as a
thin shell: a triangulated surface in 3D with a constant thickness. Darcy flow
is solved in the plane of each element, so the same code handles flat plates
in any orientation, curved panels, and multi-panel parts such as T-joints.

This is a 2.5D extension of
[RTM-flow v1.0](https://github.com/ChaffyHurdle/RTM-flow) (Wells, Park, Matveev
& Tretyakov, 2024), which in turn builds on
[MATCVFEM](https://github.com/parkmh/MATCVFEM) by M. Park and M. V. Tretyakov.
The class structure and workflow are the same as in the RTM-flow user manual;
the differences are listed [below](#differences-from-the-2d-rtm-flow-code).

## Requirements

- MATLAB R2022b or later (tested on R2025b).
- No toolboxes. Plotting uses `trisurf`, so the PDE Toolbox is not needed.

## Model

At each time level the pressure on the filled region $D_t$ satisfies

$$-\nabla\cdot(K\nabla p) = 0 \ \text{in } D_t, \qquad p = p_I \ \text{on the inlet}, \qquad n\cdot\nabla p = 0 \ \text{on sealed edges},$$

with $p = 0$ (vent pressure) on the flow front. The front moves with the Darcy
velocity $v = -K\nabla p / (\phi\mu)$. Each node has a control volume (CV);
the solver finds the next CV to fill, advances time to that point, updates
the filled region and repeats until every vent node is full.

For a 1D flow over length $L$ the fill time is
$t_f = \phi\mu L^2 / (2 K p_I)$, which the test scripts use as a check.

## Code structure

| Class | Role |
|---|---|
| `ShellMesh` | Nodes `[N x 3]`, triangles, boundary nodes; per-element local frame, areas, centroids |
| `VoronoiMesh` | Control volumes, outflow vectors, connectivity, panel-junction nodes |
| `Physics` | Viscosity, porosity, thickness, permeability function |
| `Pressure` | FEM pressure solve, inlet/vent/sealed boundary classification |
| `Velocity` | Darcy velocity per element (3D vector in the element plane) |
| `RTMFlow` | Time stepping: flow rates, time step, fill factors, front update |
| `Visualisation` | Plot switches |

Methods live in the `@ClassName` folders. As in the manual, methods that
change an object must return it: `obj = obj.method();`.

## Setting up a simulation

1. **Mesh.** Provide `nodes3d` `[N x 3]`, `elements` `[M x 3]` (1-based) and
   `boundary_nodes`, the nodes on free edges of the surface. Do not include
   nodes on internal junction lines. Mesh generators are included for a
   vertical plate, curved panel, T-shape, U-shape and cross-rib.
2. **Physics.** `Physics(viscosity, porosity, thickness, @permeability)`.
   `permeability(x)` receives the 1x3 element centroid and returns a 2x2
   tensor in the element's local frame.
3. **Boundary conditions.** `Pressure(mesh, physics, @is_inlet, @is_vent, @p_D)`.
   `is_inlet(x)` and `is_vent(x)` receive a 1x3 node position and return true
   on the inlet / vent. Use a tolerance (`abs(x(1)) < 1e-10`) rather than `==`
   for coordinates that come from arithmetic. `p_D(pressure_class)` returns the
   nodal pressure vector with the inlet value set.
4. **Run.** `RTMFlow(mesh, physics, pressure)`, set any plot options, then
   `obj = obj.run();`. The fill time is in `obj.time`.

## Examples
Examples provided in the repo.


## Visualisation

All options are off by default.

| Property | Figure / file |
|---|---|
| `is_plotting_volume` | Fill factor of each CV, figure 1 |
| `is_plotting_pressure` | Pressure, figure 2 |
| `is_plotting_volume_flow` | CV flow rates, figure 3 |
| `is_plotting_velocity` | Velocity arrows, figure 4 |
| `is_animate_volume`, `is_animate_pressure`, `is_animate_volume_flow`, `is_animate_velocity` | Append each frame to `volume.gif`, `pressure.gif`, ... |
| `is_plotting_mesh` | Draw every element edge. The geometry outline (free edges, junctions, folds) is always drawn. |

Plotting every step is slow; turn it off for production runs.

## Test cases

| Script | What it checks |
|---|---|
| `vertical_plate_test.m` | Flat plate in the x–z plane; front should stay straight |
| `curved_plate_test.m` | Curved panel, flow along and around the curve, against the 1D fill time |
| `t_shape_test.m` | Flow along a T-joint (1D check) and across it (front splits at the junction) |
| `u_shape_test.m` | U-channel (three panels, two right-angle folds); plug flow along the channel |
| `cross_rib_test.m` | Plate with two crossing ribs: plate–rib junction lines and a rib–rib column |

## Differences from the 2D RTM-flow code

- `ShellMesh(nodes3d, elements, boundary_nodes)` replaces
  `DelaunayMesh(p, e, t)`. Nodes have three coordinates.
- Each element has a local frame (`local_e1`, `local_e2`, `local_normal`).
  Shape-function gradients and CV normals are computed in this frame, and
  velocities are mapped back to 3D vectors.
- Inlet/vent functions and `permeability` receive 1x3 positions.
- Nodes where panels meet at an angle (`VoronoiMesh.is_junction_node`) are
  handled so that flow passes from one panel to the next.

## Limitations

- **Anisotropic permeability.** The local `e1` axis follows each element's
  first edge, so it changes direction from element to element. A tensor such
  as `diag(k1, k2)` is therefore not aligned with a fixed material direction.
  Use isotropic `k * eye(2)` until a fibre-direction frame is added.
- Viscosity, porosity and thickness are constants.
- Accuracy is first order in mesh size, as for the 2D code.

## References

1. H. Wells, M. Park, M. Matveev, M. Tretyakov. *RTM-flow users manual, v1.0*, 2024.
2. M. Park, M. Tretyakov. Stochastic resin transfer molding process.
   *SIAM/ASA J. Uncertainty Quantification*, 5, 2016.
3. M. Park. MATCVFEM: MATLAB toolbox for the control volume finite element
   method. https://github.com/parkmh/MATCVFEM
4. M. Iglesias, M. Park, M. V. Tretyakov. Bayesian inversion in resin transfer
   molding. *Inverse Problems*, 34(10):105002, 2018.
