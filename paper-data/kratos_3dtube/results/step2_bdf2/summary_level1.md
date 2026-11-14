# Kratos 3dTube QoIs (definitions of the solids4foam driver)

| run | level factor | scheme | dt (s) | to t (ms) | fluid tets (hex equiv.) | wall hexes | coupling solves mean/max | steps over tol | max residual (m) | wall (min) | threads |
|---|---|---|---|---|---|---|---|---|---|---|---|
| f1_dt2.5e-5_bdf2_amg | 1 | bdf2 | 2.5e-05 | 20 | 96000 (16000) | 6400 | 4.62/14 | 0 | 1.0e-09 | 155.9 | 6 |
| f1_dt1e-4_euler_dtau1_nobfs | 1 | euler | 0.0001 | 20.1 | 96000 (16000) | 6400 | 4.48/13 | 0 | 9.9e-10 | 56.2 | 4 |

| quantity | f1_dt2.5e-5_bdf2_amg | f1_dt1e-4_euler_dtau1_nobfs | Lozovskiy et al. (2019), Euler 1e-4 | Eken (2016), Euler 1e-4 | Tuković et al. (2018), backward 2.5e-5 | solids4foam Euler L1, dt 1e-4 (after the fixes; only u_r,max re-checked) | solids4foam Euler L1, dt 1e-4 (before the fixes, Section 6) | solids4foam Euler L2, dt 1e-4 (before the fixes) | solids4foam Backward L1, dt 2.5e-5 (after the fixes) | solids4foam Backward L3, dt 6.25e-6 (after the fixes) |
|---|---|---|---|---|---|---|---|---|---|---|
| u_r,max(A) (mm) | 0.16035 | 0.12334 | 0.1194 | 0.1242 | 0.15694 | 0.12241 | 0.12233 | 0.12107 | 0.15978 | 0.15851 |
| t(u_r,max) (ms) | 7.1753 | 7.1236 | 7 | 6.69 | 7.23 | – | 7.005 | 7.037 | – | 7.252 |
| u_z,min(A) (mm) | -0.081983 | -0.071536 | -0.074 | -0.0794 | – | – | -0.07601 | -0.07515 | -0.08819 | -0.08571 |
| u_r,min(A), t > 14 ms (mm) | -0.10234 | -0.076886 | -0.0755 | -0.0763 | – | – | -0.07756 | -0.07611 | – | – |
| t_arr(A) (ms) | 5.9536 | 5.3338 | – | – | – | – | – | – | 5.873 | 5.825 |
| c_wall (m/s) | 5.1261 | 5.2081 | – | – | – | – | – | – | – | – |
| c_p (m/s) | 4.5783 | 4.4971 | – | – | – | – | – | – | 4.654 | 4.594 |

Differences of the first run from the references:

| quantity | Lozovskiy et al. (2019), Euler 1e-4 | Eken (2016), Euler 1e-4 | Tuković et al. (2018), backward 2.5e-5 | solids4foam Euler L1 |
|---|---|---|---|---|
| u_r,max(A) | +34.3% | +29.1% | +2.2% | +31.0% |
| t(u_r,max) | +2.5% | +7.3% | -0.8% | +2.4% |
| u_z,min(A) | -10.8% | -3.3% | – | -7.9% |
| u_r,min(A), t > 14 ms | -35.5% | -34.1% | – | -31.9% |
| t_arr(A) | – | – | – | – |
| c_wall | – | – | – | – |
| c_p | – | – | – | – |
