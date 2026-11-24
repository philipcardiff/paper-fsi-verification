/*---------------------------------------------------------------------------*\
    Boundary conditions with the exact manufactured fluid solution
    (derive/generated/mmsFsiExact.H), evaluated at the current face centres
    and the current time:

        mmsFsiVelocity          fixedValue U
        mmsFsiVelocityGradient  fixedGradient U, gradient = n & grad(U)
        mmsFsiPressure          fixedValue p (kinematic)
        mmsFsiPressureGradient  fixedGradient p, gradient = n & grad(p)
\*---------------------------------------------------------------------------*/
#ifndef mmsFsiFvPatchFields_H
#define mmsFsiFvPatchFields_H

#include "fixedValueFvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

#define MMS_FSI_BC_DECLARE(ClassName, Base, Type, TypeNameStr)                 \
class ClassName : public Base                                                 \
{                                                                             \
public:                                                                       \
    TypeName(TypeNameStr);                                                    \
    ClassName(const fvPatch& p, const DimensionedField<Type, volMesh>& iF)    \
    : Base(p, iF) {}                                                          \
    ClassName                                                                 \
    (                                                                         \
        const fvPatch& p,                                                     \
        const DimensionedField<Type, volMesh>& iF,                            \
        const dictionary& dict                                                \
    );                                                                        \
    ClassName                                                                 \
    (                                                                         \
        const ClassName& ptf,                                                 \
        const fvPatch& p,                                                     \
        const DimensionedField<Type, volMesh>& iF,                            \
        const fvPatchFieldMapper& m                                           \
    )                                                                         \
    : Base(ptf, p, iF, m) {}                                                  \
    ClassName(const ClassName& ptf) : Base(ptf) {}                            \
    ClassName(const ClassName& ptf, const DimensionedField<Type, volMesh>& iF)\
    : Base(ptf, iF) {}                                                        \
    virtual tmp<fvPatchField<Type>> clone() const                             \
    {                                                                         \
        return tmp<fvPatchField<Type>>(new ClassName(*this));                 \
    }                                                                         \
    virtual tmp<fvPatchField<Type>> clone                                     \
    (                                                                         \
        const DimensionedField<Type, volMesh>& iF                             \
    ) const                                                                   \
    {                                                                         \
        return tmp<fvPatchField<Type>>(new ClassName(*this, iF));             \
    }                                                                         \
    virtual void updateCoeffs();                                              \
};

namespace Foam
{
MMS_FSI_BC_DECLARE(mmsFsiVelocityFvPatchVectorField, fixedValueFvPatchVectorField, vector, "mmsFsiVelocity")
MMS_FSI_BC_DECLARE(mmsFsiVelocityGradientFvPatchVectorField, fixedGradientFvPatchVectorField, vector, "mmsFsiVelocityGradient")
MMS_FSI_BC_DECLARE(mmsFsiPressureFvPatchScalarField, fixedValueFvPatchScalarField, scalar, "mmsFsiPressure")
MMS_FSI_BC_DECLARE(mmsFsiPressureGradientFvPatchScalarField, fixedGradientFvPatchScalarField, scalar, "mmsFsiPressureGradient")
}

#endif
