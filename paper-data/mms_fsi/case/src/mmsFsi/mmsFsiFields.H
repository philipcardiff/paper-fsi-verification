/*---------------------------------------------------------------------------*\
    Helpers that evaluate the generated manufactured FSI solution
    (derive/generated/mmsFsiExact.H) on OpenFOAM point lists.
\*---------------------------------------------------------------------------*/
#ifndef mmsFsiFields_H
#define mmsFsiFields_H

#include "vectorField.H"
#include "tensorField.H"
#include "symmTensorField.H"
#include "mmsFsiExact.H"

namespace Foam
{
namespace mmsFsiFields
{

struct FluidState
{
    vector U; scalar p; vector f; tensor gradU; vector gradp;
};

inline FluidState fluid(const vector& c, const scalar t)
{
    double U[2], p[1], f[2], G[4], gp[2];
    mmsFsi::fluid(c.x(), c.y(), t, U, p, f, G, gp);
    FluidState s;
    s.U = vector(U[0], U[1], 0);
    s.p = p[0];
    s.f = vector(f[0], f[1], 0);
    // gradU(i,j) = dU_j/dx_i (OpenFOAM convention); G[2*i+j] = dU_i/dx_j
    s.gradU = tensor(G[0], G[2], 0, G[1], G[3], 0, 0, 0, 0);
    s.gradp = vector(gp[0], gp[1], 0);
    return s;
}

struct SolidState
{
    vector D; vector U; vector b0; symmTensor sigma;
};

inline SolidState solid(const vector& C, const scalar t)
{
    double D[2], U[2], b0[2], s[6];
    mmsFsi::solid(C.x(), C.y(), t, D, U, b0, s);
    SolidState st;
    st.D = vector(D[0], D[1], 0);
    st.U = vector(U[0], U[1], 0);
    st.b0 = vector(b0[0], b0[1], 0);
    st.sigma = symmTensor(s[0], s[1], s[2], s[3], s[4], s[5]);
    return st;
}

//- Interface quantities at reference abscissa X: deformed position, Cauchy
//  traction (force per deformed area, from the solid on the fluid side
//  normal) and the fluid outward unit normal
inline void interface(const scalar X, const scalar t, vector& x, vector& tr, vector& n)
{
    double xx[2], tt[2], nn[2];
    mmsFsi::interface(X, t, xx, tt, nn);
    x = vector(xx[0], xx[1], 0);
    tr = vector(tt[0], tt[1], 0);
    n = vector(nn[0], nn[1], 0);
}

} // namespace mmsFsiFields
} // namespace Foam

#endif
