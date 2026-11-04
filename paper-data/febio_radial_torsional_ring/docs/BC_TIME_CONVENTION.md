# The dt/2 time-label offset of load-curve-driven fluid boundary conditions (stock v4.13)

## Finding (PROVEN: source and data)
In the stock FEBio v4.13 **fluid-module** solver a prescribed value driven by a load controller (load curve) is evaluated at the generalized-alpha level time t - (1 - alpha_f) dt, with alpha_f = 1/(1 + rho_inf). At rho_inf = 1 this is t - dt/2, so the state logged at t_n answers a problem whose wall history is delayed by dt/2. The **fluid-FSI** solver evaluates load controllers at the end of the step (no offset), and `math` expressions in boundary-condition values read the end-of-step time in both solvers.

## Source evidence (FEBio tag v4.13, commit 32ae206ff4881dfb54f62296cd1558e58ed9fcc6, clean tree; the run logs report version 4.13.0.32ae206ff)
- `FEBioFluid/FEFluidSolver.cpp`, `FEFluidSolver::InitStep` (lines 519-535): `double t = tp.currentTime; double dt = tp.timeIncrement; double ta = (t > 0) ? t - (1-m_alphaf)*dt : m_alphaf*dt; return FESolver::InitStep(ta);` with the comment "evaluate load curve values at current (or intermediate) time". `FESolver::InitStep(ta)` evaluates all load controllers at `ta`.
- Same file, lines 130-138: `m_alphaf = 1.0/(1+m_rhoi)` for 0 <= rhoi <= 1, and `m_alphaf = 1` for rhoi = -1; so rho_inf = 1 gives ta = t - dt/2 and rho_inf = 0 gives alpha_f = 1, ta = t (no offset). The fluid solver default is rhoi = 0 (line 90).
- `FEBioFluid/FEFluidFSISolver.cpp`, `FEFluidFSISolver::InitStep` (lines 746-764): the alpha-level time version is commented out (`// double ta = ...; // return FESolver::InitStep(ta);`) and the function ends with `return FESolver::InitStep(t);`.
- `FECore/FESolver.cpp`, `FESolver::InitStep(time)` (lines 346-365): `fem.EvaluateLoadControllers(time)`, then data generators and load parameters at the same time. The fluid-module solvers for solutes (`FEFluidSolutesSolver.cpp` line 577) and thermofluids (`FEThermoFluidSolver.cpp` line 653) use the same alpha-level time.
- `FECore/FEScalarValuator.cpp`, `FEMathExpression::value` (line 124): the expression variable t is `fem->GetTime().currentTime`, the end-of-step time.

## Data evidence
| Check | Result |
|---|---|
| RT-B jobs 14209 (M1) and 14210 (M2), `<bc ... lc>` wall velocity = omega(t) (-y, x) | logged wall tangential velocity equals the table value omega(t_n - 0.005) to 2.5e-13 rad/s (max over all 200 steps), whereas |logged - omega(t_n)| reaches 5.9e-3 rad/s; example t = 1.0: logged 0.036760524, table omega(1.0) = 0.031602839, table omega(0.995) = 0.036760524 |
| Diagnostic job 14211 (M1, same input with the table advanced by dt/2) | logged wall velocity equals omega(t_n) to 2.2e-13; errors against the as-run exact solution: velocity 0.186 %, T2 11.42 %, H_f 0.169 %, T3 0.26 %, identical to the post-hoc aligned analysis of job 14209 to 3-4 digits |
| Static rotcheck (job 14187, `math` expression) | evaluated at the end-of-step time exactly (7e-13 m) |

Figure `../figures/fig4_bc_time_offset.png`.

## Consequences
- As run, the offset produces a time-offset floor of about 1.7 % in the velocity error and makes the velocity-convergence gate fail (ratio 1.02); the aligned errors are 0.186 % and 0.044 % (order 2.07). All aligned numbers in this package are labelled post hoc.
- The coupled rotation BC (a `math` expression in the fluid-FSI solver) carries no offset.
- Earlier work (source check only, no new runs): cavity, Hron-Turek CFD1/FSI1 and Womersley used rho_inf = 0 (alpha_f = 1, no offset) for their fluid-module ramps; the ring runs used the fluid-FSI solver with initial-velocity drives; so none of the earlier packages is affected. Any rho_inf = 1 fluid-module run with a load curve would be offset by dt/2 (for a steady state only the ramp end time shifts).
- No change to stock FEBio is proposed; the offset is documented here so that inputs can compensate (shift the table by dt/2 or use rho_inf = 0 in the fluid module). Whether the offset is intended is a question for the FEBio developers and is not asserted here.
