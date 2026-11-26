"""ALE fluid forms in the REFERENCE (initial) configuration.

x = X + d(X,t);  F = I + grad_X d;  J = det F;  mesh velocity w = (d1-d0)/dt.
Spatial gradient of a field f:  grad_x f = grad_X f F^{-1}.

Time discretisation (theta-scheme, theta=1 BE, 1/2 CN) of the conservative ALE form
    d/dt( J u ) + J div_x( (u-w) (x) u ) ,   using  J div_x w := (J1-J0)/dt  (discrete GCL).
So  time + mesh-motion terms = rho [ J1 u1 - J0 u0 - (J1-J0) u_theta ] / dt . phi
which preserves a uniform flow exactly on a moving mesh (GCL satisfied by construction).
Convection/viscous terms are evaluated at level theta using each level's own (F,J).
Pressure and continuity at t_{n+1}.
"""
import ufl


def kin(d):
    I = ufl.Identity(len(d))
    F = I + ufl.grad(d)
    return F, ufl.det(F)


def gradx(f, F):
    return ufl.dot(ufl.grad(f), ufl.inv(F))


def divx(f, F):
    return ufl.tr(gradx(f, F))


def _spatial(u, d, w, rho, mu, phi, F, J):
    gu = gradx(u, F)
    gp = gradx(phi, F)
    conv = rho * J * ufl.dot(ufl.dot(gu, u - w), phi)
    visc = 2.0 * mu * J * ufl.inner(ufl.sym(gu), ufl.sym(gp))
    return conv + visc


def fluid_ale_residual(u1, p1, d1, u0, d0, dt, theta, rho, mu, phi, q, dx, f_body=None, w_levels=None, midpoint=False, p_geom_theta=False):
    """Returns UFL residual (sum of integrals over `dx`). f_body: body force per unit volume
    on the CURRENT configuration (e.g. for manufactured solutions), may be None."""
    F1, J1 = kin(d1)
    F0, J0 = kin(d0)
    w = (d1 - d0) / dt
    w1, w0l = (w, w) if w_levels is None else w_levels   # mesh velocities at levels n+1 and n
    u_th = theta * u1 + (1 - theta) * u0
    time = rho * ufl.dot(J1 * u1 - J0 * u0 - (J1 - J0) * u_th, phi) / dt
    res = time * dx
    if midpoint:
        # all ALE geometry terms at the theta-level state: d_th, u_th, F(d_th), J(d_th); w = (d1-d0)/dt
        d_th = theta * d1 + (1 - theta) * d0
        Fm, Jm_ = kin(d_th)
        res += _spatial(u_th, d_th, w, rho, mu, phi, Fm, Jm_) * dx
    else:
        res += theta * _spatial(u1, d1, w1, rho, mu, phi, F1, J1) * dx
        if theta != 1.0:
            res += (1 - theta) * _spatial(u0, d0, w0l, rho, mu, phi, F0, J0) * dx
    if p_geom_theta:
        res += -p1 * (theta * J1 * divx(phi, F1) + (1 - theta) * J0 * divx(phi, F0)) * dx
    else:
        res += -p1 * J1 * divx(phi, F1) * dx
    res += q * J1 * divx(u1, F1) * dx
    if f_body is not None:
        res += -rho * J1 * ufl.dot(f_body, phi) * dx
    return res


def harmonic_residual(d, psi, dx, k=1.0):
    """Harmonic (Laplace) mesh motion with stiffness weight k (UFL scalar, e.g. DG0 1/vol)."""
    return k * ufl.inner(ufl.grad(d), ufl.grad(psi)) * dx


def fluid_ale_residual_bdf2(u1, p1, d1, u0, d0, um, dm, dt, a0, a1, a2, rho, mu, phi, q, dx, f_body=None):
    """BDF2 variable-by-flag version: coefficients (a0,a1,a2) = (1,-1,0) BE start, (3/2,-2,1/2) BDF2 (Constants).
    Conservative ALE time term with discrete GCL: d_t(Ju) - u d_t J, d_t z := (a0 z1 + a1 z0 + a2 zm)/dt,
    mesh velocity w = d_t d. Spatial terms (convection, viscous, pressure, continuity) at t_{n+1}."""
    F1, J1 = kin(d1)
    _, J0 = kin(d0)
    _, Jm = kin(dm)
    w = (a0 * d1 + a1 * d0 + a2 * dm) / dt
    dJu = (a0 * J1 * u1 + a1 * J0 * u0 + a2 * Jm * um) / dt
    dJ = (a0 * J1 + a1 * J0 + a2 * Jm) / dt
    res = rho * ufl.dot(dJu - u1 * dJ, phi) * dx
    res += _spatial(u1, d1, w, rho, mu, phi, F1, J1) * dx
    res += -p1 * J1 * divx(phi, F1) * dx + q * J1 * divx(u1, F1) * dx
    if f_body is not None:
        res += -rho * J1 * ufl.dot(f_body, phi) * dx
    return res
