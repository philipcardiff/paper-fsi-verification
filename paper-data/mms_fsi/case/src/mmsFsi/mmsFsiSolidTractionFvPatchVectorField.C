#include "mmsFsiSolidTractionFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "mmsFsiFields.H"

namespace Foam
{
void mmsFsiSolidTractionFvPatchVectorField::updateCoeffs()
{
    if (updated()) return;
    const scalar t = db().time().value();
    const vectorField& Cf = patch().Cf();
    vectorField tr(Cf.size());
    forAll(Cf, i)
    {
        vector x, trf, n;
        mmsFsiFields::interface(Cf[i].x(), t, x, trf, n);
        tr[i] = -trf;
        if (nominal_)
        {
            // Deformed length per unit reference length of the interface
            const scalar eps = 1e-6;
            vector xp, xm, dum1, dum2;
            mmsFsiFields::interface(Cf[i].x() + eps, t, xp, dum1, dum2);
            mmsFsiFields::interface(Cf[i].x() - eps, t, xm, dum1, dum2);
            tr[i] *= mag(xp - xm)/(2*eps);
        }
    }
    traction() = tr;
    pressure() = 0;
    solidTractionFvPatchVectorField::updateCoeffs();
}
makePatchTypeField(fvPatchVectorField, mmsFsiSolidTractionFvPatchVectorField);
}
