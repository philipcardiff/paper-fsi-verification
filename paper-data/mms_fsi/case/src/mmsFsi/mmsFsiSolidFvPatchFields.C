#include "mmsFsiSolidFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "mmsFsiFields.H"

namespace Foam
{

void mmsFsiSolidDisplacementFvPatchVectorField::updateCoeffs()
{
    if (updated()) return;
    const scalar t = db().time().value();
    const vectorField& Cf = patch().Cf();   // reference (the solid mesh is not moved)
    vectorField disp(Cf.size());
    forAll(Cf, i) { disp[i] = mmsFsiFields::solid(Cf[i], t).D; }
    totalDisp() = disp;
    fixedDisplacementFvPatchVectorField::updateCoeffs();
}

void mmsFsiSolidVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated()) return;
    const scalar t = db().time().value();
    const vectorField& Cf = patch().Cf();
    vectorField v(Cf.size());
    forAll(Cf, i) { v[i] = mmsFsiFields::solid(Cf[i], t).U; }
    fvPatchVectorField::operator==(v);
    fixedValueFvPatchVectorField::updateCoeffs();
}

makePatchTypeField(fvPatchVectorField, mmsFsiSolidDisplacementFvPatchVectorField);
makePatchTypeField(fvPatchVectorField, mmsFsiSolidVelocityFvPatchVectorField);

} // namespace Foam
