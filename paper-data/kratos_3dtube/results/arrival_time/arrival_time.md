# 3dTube arrival times, Kratos against solids4foam

Front rules of the solids4foam driver: wall arrival is the first upward crossing of 50% of the station's incident peak (the first local maximum above 50% of the history maximum); pressure arrival is the first crossing of 50% of the inlet pressure.

## Wall arrival per station (ms), matched levels

| z (mm) | Kratos L1 | s4f L1 | diff | Kratos L2 | s4f L2 | diff | s4f L3 |
|---|---|---|---|---|---|---|---|
| 5 | 1.195 | 1.181 | +0.013 | 1.206 | 1.194 | +0.012 | 1.201 |
| 10 | 2.432 | 2.401 | +0.030 | 2.456 | 2.435 | +0.021 | 2.448 |
| 15 | 3.656 | 3.612 | +0.044 | 3.685 | 3.641 | +0.044 | 3.659 |
| 20 | 4.843 | 4.786 | +0.057 | 4.876 | 4.827 | +0.049 | 4.839 |
| 25 | 5.954 | 5.873 | +0.081 | 5.925 | 5.843 | +0.083 | 5.825 |
| 30 | 6.698 | 6.463 | +0.235 | 6.529 | 6.415 | +0.115 | 6.426 |
| 35 | 7.072 | 7.480 | -0.408 | 7.622 | 7.251 | +0.371 | 7.483 |
| 40 | 8.476 | 8.406 | +0.071 | 8.532 | 8.492 | +0.039 | 8.521 |
| 45 | 9.415 | 8.182 | +1.233 | 9.576 | 8.341 | +1.235 | 8.380 |

Beyond 30 mm the incident front merges with the reflection from the clamped outlet and the first-peak rule can pick different features in the two codes.

## Wall-front speed between z = 5 and 25 mm (m/s)

| | level 1 | level 2 | level 3 |
|---|---|---|---|
| Kratos | 4.203 | 4.238 | – |
| solids4foam | 4.263 | 4.302 | 4.326 |
| Kratos / solids4foam − 1 | -1.42% | -1.50% | – |

## Level 1 against the solids4foam level-1 rerun (full histories to 10 ms)

Crossing-time difference Kratos − solids4foam (ms) at fractions of each code's own incident peak; the least-squares time shift that aligns the two incident fronts over the 3 ms before the solids4foam peak; and the pressure arrival (the solids4foam axis probes take the value of the containing cell, so they carry up to half a cell of position error).

| z (mm) | f = 0.1 | f = 0.25 | f = 0.5 | f = 0.75 | f = 0.9 | peak K / s4f (µm) | front shift (ms) | p arrival K / s4f (ms) |
|---|---|---|---|---|---|---|---|---|
| 5 | +0.013 | -0.034 | +0.013 | +0.032 | +0.070 | 150.99 / 148.24 | -0.025 | 1.021 / 0.953 |
| 10 | +0.311 | -0.010 | +0.030 | +0.048 | +0.062 | 157.98 / 157.13 | +0.025 | 2.121 / 2.027 |
| 15 | -0.102 | +0.040 | +0.044 | +0.081 | +0.106 | 163.28 / 160.28 | +0.043 | 3.207 / 3.077 |
| 20 | -0.073 | +0.265 | +0.057 | +0.096 | +0.087 | 163.78 / 162.11 | +0.058 | 4.261 / 4.110 |
| 25 | +0.090 | +0.012 | +0.082 | +0.100 | +0.034 | 160.35 / 159.76 | +0.068 | 5.295 / 5.121 |
| 30 | +0.522 | +0.008 | +0.235 | +0.129 | -0.012 | 150.86 / 150.41 | +0.068 | 6.342 / 6.161 |

Tuković et al. (2018), Fig. 25(a), same rule: t_arr(A) = 5.889 ms.

