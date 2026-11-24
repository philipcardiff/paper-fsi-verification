#include "mmsFsiPressureOracleFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "mmsFsiFields.H"

namespace Foam
{

void mmsFsiPressureOracleFvPatchScalarField::evaluate(const Pstream::commsTypes ct)
{
    fixedGradientFvPatchScalarField::evaluate(ct);

    const scalar t = db().time().value();
    const vectorField& Cf = patch().Cf();
    scalarField& v = *this;

    if (mode_ == "exact")
    {
        forAll(Cf, i)
        {
            v[i] = mmsFsiFields::fluid(Cf[i], t).p;
        }
    }
    else
    {
        // exact non-orthogonal correction k & grad(p) added to the
        // uncorrected value p_P + g (n.d)
        const vectorField n(patch().nf());
        const vectorField d(patch().Cf() - patch().Cn());
        forAll(Cf, i)
        {
            const vector k = d[i] - n[i]*(n[i] & d[i]);
            v[i] += k & mmsFsiFields::fluid(Cf[i], t).gradp;
        }
    }
}

makePatchTypeField(fvPatchScalarField, mmsFsiPressureOracleFvPatchScalarField);

} // namespace Foam
