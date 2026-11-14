"""3-D compressible neo-Hookean, total Lagrangian (full 3x3, no plane-strain split).

Law (equal to solids4foam neoHookeanElasticMechanicalConstitutiveLaw at 517a44921, l.213-222):
    W(F)  = mu/2 (Ibar1 - 3) + K/4 (J^2 - 1 - 2 ln J),  Ibar1 = J^(-2/3) tr(F F^T)
    sigma = (mu/J) dev(bbar) + K/2 (J^2-1)/J I,          bbar = J^(-2/3) F F^T
    P = dW/dF = mu J^(-2/3) (F - Ibar1' /3 F^-T) + K/2 (J^2-1) F^-T,  Ibar1' = tr(F F^T)   (= J sigma F^-T)
K = 2 mu (1+nu)/(3 (1-2 nu)).  Hessenthaler: mu = 64.2 kPa, nu = 0.45 -> K = 620.6 kPa.
"""
import numpy as np
import ufl

MU, NU = 64.2e3, 0.45
def bulk(mu, nu): return 2.0 * mu * (1.0 + nu) / (3.0 * (1.0 - 2.0 * nu))
K_BULK = bulk(MU, NU)

def energy(F, mu=MU, K=K_BULK):
    J = ufl.det(F)
    return 0.5 * mu * (J ** (-2.0 / 3.0) * ufl.tr(F * F.T) - 3.0) + 0.25 * K * (J ** 2 - 1.0 - 2.0 * ufl.ln(J))

def P_ufl(F, mu=MU, K=K_BULK):
    Fv = ufl.variable(F)
    return ufl.diff(energy(Fv, mu, K), Fv)

def P_analytic(F, mu=MU, K=K_BULK):
    J = ufl.det(F); FiT = ufl.inv(F).T
    return mu * J ** (-2.0 / 3.0) * (F - ufl.tr(F * F.T) / 3.0 * FiT) + 0.5 * K * (J ** 2 - 1.0) * FiT

# numpy reference
def np_W(F, mu=MU, K=K_BULK):
    J = np.linalg.det(F)
    return 0.5 * mu * (J ** (-2 / 3) * np.trace(F @ F.T) - 3) + 0.25 * K * (J ** 2 - 1 - 2 * np.log(J))
def np_P(F, mu=MU, K=K_BULK):
    J = np.linalg.det(F); FiT = np.linalg.inv(F).T
    return mu * J ** (-2 / 3) * (F - np.trace(F @ F.T) / 3 * FiT) + 0.5 * K * (J ** 2 - 1) * FiT
def np_cauchy(F, mu=MU, K=K_BULK):
    J = np.linalg.det(F); b = J ** (-2 / 3) * F @ F.T
    return mu / J * (b - np.trace(b) / 3 * np.eye(3)) + 0.5 * K * (J ** 2 - 1) / J * np.eye(3)
def np_tangent_fd(F, h=1e-6, mu=MU, K=K_BULK):
    A = np.zeros((9, 9))
    for k in range(9):
        E = np.zeros(9); E[k] = h; E = E.reshape(3, 3)
        A[:, k] = ((np_P(F + E, mu, K) - np_P(F - E, mu, K)) / (2 * h)).ravel()
    return A
