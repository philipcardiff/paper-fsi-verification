# 3dTube verification summary

## Mesh study (mesh and time step refined together; robin, backward, quick smoke test)

- Result: PASS
- Analytical wave-speed estimates (m/s): mk_inner_radius 5.477, mk_mean_radius 5.222, mk_inner_radius_tethered 5.742, thick_wall 5.066, thick_wall_inertia 4.811

| Level | Fluid cells | Solid cells | Δt (s) | Cores | u_r,max(A) (mm) | t(u_r,max) (ms) | u_z,min(A) (mm) | u_r,min(A), t > 14 ms (mm) | t_arr(A) (ms) | c_wall (m/s) | c_p (m/s) | Mean FSI iter. | Clock (s) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 16000 | 6400 | 2.5e-05 | 1 | 0.15976 | 7.202 | -0.08819 | n/a | 5.872 | 5.066 | 4.663 | 5.22 | 589 |

Checks:

- PASS: all 1 run(s) completed, and every time step met the robin convergence criteria
- u_r(A) history vs Tukovic et al. (2018): largest difference 8.2% of the published peak (information)
- u_r(A) history vs Lozovskiy et al. (2019): largest difference 35.3% of the published peak (information)
- u_r(A) history vs Eken (2016): largest difference 39.5% of the published peak (information)

- Data: `mesh_robin_backward_quick.csv`

