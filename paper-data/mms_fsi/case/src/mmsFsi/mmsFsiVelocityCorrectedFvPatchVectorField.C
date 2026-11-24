#include "mmsFsiVelocityCorrectedFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "mmsFsiFields.H"
namespace Foam
{
void mmsFsiVelocityCorrectedFvPatchVectorField::updateCoeffs()
{
    if (updated()) return;
    const scalar t = db().time().value();
    const vectorField& Cf = patch().Cf();
    vectorField v(Cf.size());
    forAll(Cf, i) { v[i] = mmsFsiFields::fluid(Cf[i], t).U; }
    fvPatchVectorField::operator==(v);
    fixedValueCorrectedFvPatchVectorField::updateCoeffs();
}
makePatchTypeField(fvPatchVectorField, mmsFsiVelocityCorrectedFvPatchVectorField);
}
