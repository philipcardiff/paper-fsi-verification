/*---------------------------------------------------------------------------*\
Class
    Foam::fv::mmsFsiSource

Description
    Manufactured-solution momentum source for the 2-D ALE FSI test
    (derive/mms_fsi.py).

    region fluid: f(x, t) per unit mass (kinematic momentum equation),
                  evaluated at the current (moving) cell centres;
    region solid: b0(X, t) per unit reference volume (total Lagrangian),
                  evaluated at the cell centres of the undeformed solid mesh.

    The source is explicit and is re-evaluated when the time or the cell
    centres change. It is added to both the plain and the density-weighted
    equations, with the dimensions the equation requires.

Usage
    \verbatim
    mmsSource
    {
        type        mmsFsiSource;
        region      fluid;    // or solid
        fields      (U);      // (D) for the displacement solids
    }
    \endverbatim
\*---------------------------------------------------------------------------*/
#ifndef mmsFsiSource_H
#define mmsFsiSource_H

#include "fvOption.H"
#include "volFields.H"

namespace Foam
{
namespace fv
{

class mmsFsiSource
:
    public fv::option
{
    word region_;
    mutable scalar cachedTime_;
    mutable vectorField cachedC_;
    mutable vectorField cachedS_;

    const vectorField& source() const;

public:

    TypeName("mmsFsiSource");

    mmsFsiSource
    (
        const word& name,
        const word& modelType,
        const dictionary& dict,
        const fvMesh& mesh
    );

    virtual ~mmsFsiSource() = default;

    virtual void addSup(fvMatrix<vector>& eqn, const label fieldi);

    virtual void addSup
    (
        const volScalarField& rho,
        fvMatrix<vector>& eqn,
        const label fieldi
    );

    virtual bool read(const dictionary& dict);
};

} // namespace fv
} // namespace Foam

#endif
