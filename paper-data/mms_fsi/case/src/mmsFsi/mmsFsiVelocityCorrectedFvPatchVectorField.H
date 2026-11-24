/*---------------------------------------------------------------------------*\
    mmsFsiVelocityCorrected: exact manufactured velocity at the face centres
    on a solids4foam fixedValueCorrected condition (non-orthogonal correction
    of snGrad and of the Laplacian boundary coefficients with grad(U)).
\*---------------------------------------------------------------------------*/
#ifndef mmsFsiVelocityCorrectedFvPatchVectorField_H
#define mmsFsiVelocityCorrectedFvPatchVectorField_H
#include "fixedValueCorrectedFvPatchVectorField.H"
namespace Foam
{
class mmsFsiVelocityCorrectedFvPatchVectorField
:
    public fixedValueCorrectedFvPatchVectorField
{
public:
    TypeName("mmsFsiVelocityCorrected");
    mmsFsiVelocityCorrectedFvPatchVectorField
    (const fvPatch& p, const DimensionedField<vector, volMesh>& iF)
    : fixedValueCorrectedFvPatchVectorField(p, iF) {}
    mmsFsiVelocityCorrectedFvPatchVectorField
    (const fvPatch& p, const DimensionedField<vector, volMesh>& iF, const dictionary& dict)
    : fixedValueCorrectedFvPatchVectorField(p, iF, dict) {}
    mmsFsiVelocityCorrectedFvPatchVectorField
    (
        const mmsFsiVelocityCorrectedFvPatchVectorField& ptf,
        const fvPatch& p, const DimensionedField<vector, volMesh>& iF,
        const fvPatchFieldMapper& m
    )
    : fixedValueCorrectedFvPatchVectorField(ptf, p, iF, m) {}
    mmsFsiVelocityCorrectedFvPatchVectorField(const mmsFsiVelocityCorrectedFvPatchVectorField& ptf)
    : fixedValueCorrectedFvPatchVectorField(ptf) {}
    mmsFsiVelocityCorrectedFvPatchVectorField
    (const mmsFsiVelocityCorrectedFvPatchVectorField& ptf, const DimensionedField<vector, volMesh>& iF)
    : fixedValueCorrectedFvPatchVectorField(ptf, iF) {}
    virtual tmp<fvPatchField<vector>> clone() const
    { return tmp<fvPatchField<vector>>(new mmsFsiVelocityCorrectedFvPatchVectorField(*this)); }
    virtual tmp<fvPatchField<vector>> clone(const DimensionedField<vector, volMesh>& iF) const
    { return tmp<fvPatchField<vector>>(new mmsFsiVelocityCorrectedFvPatchVectorField(*this, iF)); }
    virtual void updateCoeffs();
};
}
#endif
