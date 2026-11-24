#include "mmsFsiFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "mmsFsiFields.H"

namespace Foam
{

// ---- constructors from dictionary: evaluate immediately -------------------

mmsFsiVelocityFvPatchVectorField::mmsFsiVelocityFvPatchVectorField
(
    const fvPatch& p, const DimensionedField<vector, volMesh>& iF, const dictionary& dict
)
: fixedValueFvPatchVectorField(p, iF)
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}

mmsFsiVelocityGradientFvPatchVectorField::mmsFsiVelocityGradientFvPatchVectorField
(
    const fvPatch& p, const DimensionedField<vector, volMesh>& iF, const dictionary& dict
)
: fixedGradientFvPatchVectorField(p, iF)
{
    gradient() = Zero;
    fvPatchVectorField::operator=(patchInternalField());
}

mmsFsiPressureFvPatchScalarField::mmsFsiPressureFvPatchScalarField
(
    const fvPatch& p, const DimensionedField<scalar, volMesh>& iF, const dictionary& dict
)
: fixedValueFvPatchScalarField(p, iF)
{
    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));
}

mmsFsiPressureGradientFvPatchScalarField::mmsFsiPressureGradientFvPatchScalarField
(
    const fvPatch& p, const DimensionedField<scalar, volMesh>& iF, const dictionary& dict
)
: fixedGradientFvPatchScalarField(p, iF)
{
    gradient() = Zero;
    fvPatchScalarField::operator=(patchInternalField());
}

// ---- updateCoeffs -----------------------------------------------------------

void mmsFsiVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated()) return;
    const scalar t = db().time().value();
    const vectorField& Cf = patch().Cf();
    vectorField v(Cf.size());
    forAll(Cf, i) { v[i] = mmsFsiFields::fluid(Cf[i], t).U; }
    fvPatchVectorField::operator==(v);
    fixedValueFvPatchVectorField::updateCoeffs();
}

void mmsFsiVelocityGradientFvPatchVectorField::updateCoeffs()
{
    if (updated()) return;
    const scalar t = db().time().value();
    const vectorField& Cf = patch().Cf();
    const vectorField n(patch().nf());
    forAll(Cf, i) { gradient()[i] = n[i] & mmsFsiFields::fluid(Cf[i], t).gradU; }
    fixedGradientFvPatchVectorField::updateCoeffs();
}

void mmsFsiPressureFvPatchScalarField::updateCoeffs()
{
    if (updated()) return;
    const scalar t = db().time().value();
    const vectorField& Cf = patch().Cf();
    scalarField v(Cf.size());
    forAll(Cf, i) { v[i] = mmsFsiFields::fluid(Cf[i], t).p; }
    fvPatchScalarField::operator==(v);
    fixedValueFvPatchScalarField::updateCoeffs();
}

void mmsFsiPressureGradientFvPatchScalarField::updateCoeffs()
{
    if (updated()) return;
    const scalar t = db().time().value();
    const vectorField& Cf = patch().Cf();
    const vectorField n(patch().nf());
    forAll(Cf, i) { gradient()[i] = n[i] & mmsFsiFields::fluid(Cf[i], t).gradp; }
    fixedGradientFvPatchScalarField::updateCoeffs();
}

makePatchTypeField(fvPatchVectorField, mmsFsiVelocityFvPatchVectorField);
makePatchTypeField(fvPatchVectorField, mmsFsiVelocityGradientFvPatchVectorField);
makePatchTypeField(fvPatchScalarField, mmsFsiPressureFvPatchScalarField);
makePatchTypeField(fvPatchScalarField, mmsFsiPressureGradientFvPatchScalarField);

} // namespace Foam
