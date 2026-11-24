/*---------------------------------------------------------------------------*\
Class
    Foam::mmsFsiErrorsFunctionObject

Description
    Error norms of the manufactured FSI solution at every time step, written
    to postProcessing/mmsFsiErrors/errors.dat (one line per step):

      fluid U, p     : cells of the moving fluid mesh, volume-weighted
                       L1, L2 and LInf of |error|; p is kinematic. pMean is
                       the L2 error of p after removing the mean difference
      solid D, U     : cells of the (undeformed) solid mesh, volume-weighted
      intD           : solid interface faces, D vs exact D(X), area-weighted
      intTrac        : solid interface faces, applied traction
                       (solidTraction: traction - pressure*n) vs the exact
                       Cauchy traction sigma_s . n_solid at the face
      intP           : fluid interface faces, p vs exact p at the current
                       face centres

Usage
    \verbatim
    mmsErrors
    {
        type            mmsFsiErrors;
        libs            (mmsFsi);
        fluidRegion     fluid;
        solidRegion     solid;
        interfacePatch  interface;
        writeFields     false;
    }
    \endverbatim
\*---------------------------------------------------------------------------*/
#ifndef mmsFsiErrorsFunctionObject_H
#define mmsFsiErrorsFunctionObject_H

#include "functionObject.H"
#include "OFstream.H"
#include "scalarField.H"
#include "autoPtr.H"

namespace Foam
{
class Time;

class mmsFsiErrorsFunctionObject
:
    public functionObject
{
    const Time& time_;
    word fluidRegion_;
    word solidRegion_;
    word interfacePatch_;
    bool writeFields_;
    bool accelerationErrors_;
    autoPtr<OFstream> filePtr_;
        scalarField eTrN_, eTrT_, trMagSf_;

    bool compute();

public:

    TypeName("mmsFsiErrors");

    mmsFsiErrorsFunctionObject
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    virtual ~mmsFsiErrorsFunctionObject() = default;

    virtual bool read(const dictionary& dict);
    virtual bool execute();
    virtual bool write();
};

} // namespace Foam

#endif
