/*---------------------------------------------------------------------------*\
    Exact manufactured solid data on solid patches (reference face centres):

        mmsFsiSolidDisplacement  fixedDisplacement D (with the solids4foam
                                 boundary non-orthogonal correction)
        mmsFsiSolidVelocity      fixedValue U (velocity-form solid)
\*---------------------------------------------------------------------------*/
#ifndef mmsFsiSolidFvPatchFields_H
#define mmsFsiSolidFvPatchFields_H

#include "fixedDisplacementFvPatchVectorField.H"
#include "fixedValueFvPatchFields.H"

namespace Foam
{

class mmsFsiSolidDisplacementFvPatchVectorField
:
    public fixedDisplacementFvPatchVectorField
{
public:
    TypeName("mmsFsiSolidDisplacement");

    mmsFsiSolidDisplacementFvPatchVectorField
    (const fvPatch& p, const DimensionedField<vector, volMesh>& iF)
    : fixedDisplacementFvPatchVectorField(p, iF) {}

    mmsFsiSolidDisplacementFvPatchVectorField
    (const fvPatch& p, const DimensionedField<vector, volMesh>& iF, const dictionary& dict)
    : fixedDisplacementFvPatchVectorField(p, iF, dict) {}

    mmsFsiSolidDisplacementFvPatchVectorField
    (
        const mmsFsiSolidDisplacementFvPatchVectorField& ptf,
        const fvPatch& p, const DimensionedField<vector, volMesh>& iF,
        const fvPatchFieldMapper& m
    )
    : fixedDisplacementFvPatchVectorField(ptf, p, iF, m) {}

    mmsFsiSolidDisplacementFvPatchVectorField(const mmsFsiSolidDisplacementFvPatchVectorField& ptf)
    : fixedDisplacementFvPatchVectorField(ptf) {}

    mmsFsiSolidDisplacementFvPatchVectorField
    (const mmsFsiSolidDisplacementFvPatchVectorField& ptf, const DimensionedField<vector, volMesh>& iF)
    : fixedDisplacementFvPatchVectorField(ptf, iF) {}

    virtual tmp<fvPatchField<vector>> clone() const
    { return tmp<fvPatchField<vector>>(new mmsFsiSolidDisplacementFvPatchVectorField(*this)); }

    virtual tmp<fvPatchField<vector>> clone(const DimensionedField<vector, volMesh>& iF) const
    { return tmp<fvPatchField<vector>>(new mmsFsiSolidDisplacementFvPatchVectorField(*this, iF)); }

    virtual void updateCoeffs();
};


class mmsFsiSolidVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
public:
    TypeName("mmsFsiSolidVelocity");

    mmsFsiSolidVelocityFvPatchVectorField
    (const fvPatch& p, const DimensionedField<vector, volMesh>& iF)
    : fixedValueFvPatchVectorField(p, iF) {}

    mmsFsiSolidVelocityFvPatchVectorField
    (const fvPatch& p, const DimensionedField<vector, volMesh>& iF, const dictionary& dict)
    : fixedValueFvPatchVectorField(p, iF, dict) {}

    mmsFsiSolidVelocityFvPatchVectorField
    (
        const mmsFsiSolidVelocityFvPatchVectorField& ptf,
        const fvPatch& p, const DimensionedField<vector, volMesh>& iF,
        const fvPatchFieldMapper& m
    )
    : fixedValueFvPatchVectorField(ptf, p, iF, m) {}

    mmsFsiSolidVelocityFvPatchVectorField(const mmsFsiSolidVelocityFvPatchVectorField& ptf)
    : fixedValueFvPatchVectorField(ptf) {}

    mmsFsiSolidVelocityFvPatchVectorField
    (const mmsFsiSolidVelocityFvPatchVectorField& ptf, const DimensionedField<vector, volMesh>& iF)
    : fixedValueFvPatchVectorField(ptf, iF) {}

    virtual tmp<fvPatchField<vector>> clone() const
    { return tmp<fvPatchField<vector>>(new mmsFsiSolidVelocityFvPatchVectorField(*this)); }

    virtual tmp<fvPatchField<vector>> clone(const DimensionedField<vector, volMesh>& iF) const
    { return tmp<fvPatchField<vector>>(new mmsFsiSolidVelocityFvPatchVectorField(*this, iF)); }

    virtual void updateCoeffs();
};

} // namespace Foam

#endif
