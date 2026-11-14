# Kratos 3dTube QoIs (definitions of the solids4foam driver)

| run | level factor | scheme | dt (s) | to t (ms) | fluid tets (hex equiv.) | wall hexes | coupling solves mean/max | steps over tol | max residual (m) | wall (min) | threads |
|---|---|---|---|---|---|---|---|---|---|---|---|
| f1_dt1e-4_euler_bfs | 1 | euler | 0.0001 | 20.1 | 96000 (16000) | 6400 | 4.74/15 | 0 | 9.9e-10 | 60.0 | 8 |
| f1_dt1e-4_euler_nobfs | 1 | euler | 0.0001 | 10.3 | 96000 (16000) | 6400 | 8.77/50 | 2 | 7.3e-08 | running | None |

| quantity | f1_dt1e-4_euler_bfs | f1_dt1e-4_euler_nobfs | Lozovskiy et al. (2019), Euler 1e-4 | Eken (2016), Euler 1e-4 | Tuković et al. (2018), backward 2.5e-5 | solids4foam Euler L1, dt 1e-4 (after the fixes; only u_r,max re-checked) | solids4foam Euler L1, dt 1e-4 (before the fixes, Section 6) | solids4foam Euler L2, dt 1e-4 (before the fixes) | solids4foam Backward L1, dt 2.5e-5 (after the fixes) | solids4foam Backward L3, dt 6.25e-6 (after the fixes) |
|---|---|---|---|---|---|---|---|---|---|---|
| u_r,max(A) (mm) | 0.1116 | 0.1145 | 0.1194 | 0.1242 | 0.15694 | 0.12241 | 0.12233 | 0.12107 | 0.15978 | 0.15851 |
| t(u_r,max) (ms) | 7.0255 | 7.0415 | 7 | 6.69 | 7.23 | – | 7.005 | 7.037 | – | 7.252 |
| u_z,min(A) (mm) | -0.068005 | -0.069657 | -0.074 | -0.0794 | – | – | -0.07601 | -0.07515 | -0.08819 | -0.08571 |
| u_r,min(A), t > 14 ms (mm) | -0.068093 | – | -0.0755 | -0.0763 | – | – | -0.07756 | -0.07611 | – | – |
| t_arr(A) (ms) | 5.1589 | 5.1753 | – | – | – | – | – | – | 5.873 | 5.825 |
| c_wall (m/s) | 5.3374 | 5.3359 | – | – | – | – | – | – | – | – |
| c_p (m/s) | 4.3702 | 4.413 | – | – | – | – | – | – | 4.654 | 4.594 |

Differences of the first run from the references:

| quantity | Lozovskiy et al. (2019), Euler 1e-4 | Eken (2016), Euler 1e-4 | Tuković et al. (2018), backward 2.5e-5 | solids4foam Euler L1 |
|---|---|---|---|---|
| u_r,max(A) | -6.5% | -10.1% | -28.9% | -8.8% |
| t(u_r,max) | +0.4% | +5.0% | -2.8% | +0.3% |
| u_z,min(A) | +8.1% | +14.4% | – | +10.5% |
| u_r,min(A), t > 14 ms | +9.8% | +10.8% | – | +12.2% |
| t_arr(A) | – | – | – | – |
| c_wall | – | – | – | – |
| c_p | – | – | – | – |
