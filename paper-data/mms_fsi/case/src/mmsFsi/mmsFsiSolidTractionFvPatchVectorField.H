/*---------------------------------------------------------------------------*\
    mmsFsiSolidTraction: solidTraction with the exact interface traction that
    the fluid exerts on the solid (as the QM coupling would pass it on with an
    exact fluid traction): traction = -t_f(X, t), pressure = 0, at the
    reference face centres X.
\*---------------------------------------------------------------------------*/
#ifndef mmsFsiSolidTractionFvPatchVectorField_H
#define mmsFsiSolidTractionFvPatchVectorField_H
#include "solidTractionFvPatchVectorField.H"
namespace Foam
{
class mmsFsiSolidTractionFvPatchVectorField
:
    public solidTractionFvPatchVectorField
{
    //- Apply the exact nominal traction (per unit reference area):
    //  -t_f |dx_G/dX|; use with useUndeformedArea yes
    bool nominal_ = false;
public:
    TypeName("mmsFsiSolidTraction");
    mmsFsiSolidTractionFvPatchVectorField
    (const fvPatch& p, const DimensionedField<vector, volMesh>& iF)
    : solidTractionFvPatchVectorField(p, iF) {}
    mmsFsiSolidTractionFvPatchVectorField
    (const fvPatch& p, const DimensionedField<vector, volMesh>& iF, const dictionary& dict)
    : solidTractionFvPatchVectorField(p, iF, dict),
      nominal_(dict.lookupOrDefault<bool>("nominalTraction", false)) {}
    mmsFsiSolidTractionFvPatchVectorField
    (
        const mmsFsiSolidTractionFvPatchVectorField& ptf,
        const fvPatch& p, const DimensionedField<vector, volMesh>& iF,
        const fvPatchFieldMapper& m
    )
    : solidTractionFvPatchVectorField(ptf, p, iF, m), nominal_(ptf.nominal_) {}
    mmsFsiSolidTractionFvPatchVectorField(const mmsFsiSolidTractionFvPatchVectorField& ptf)
    : solidTractionFvPatchVectorField(ptf), nominal_(ptf.nominal_) {}
    mmsFsiSolidTractionFvPatchVectorField
    (const mmsFsiSolidTractionFvPatchVectorField& ptf, const DimensionedField<vector, volMesh>& iF)
    : solidTractionFvPatchVectorField(ptf, iF), nominal_(ptf.nominal_) {}
    virtual tmp<fvPatchField<vector>> clone() const
    { return tmp<fvPatchField<vector>>(new mmsFsiSolidTractionFvPatchVectorField(*this)); }
    virtual tmp<fvPatchField<vector>> clone(const DimensionedField<vector, volMesh>& iF) const
    { return tmp<fvPatchField<vector>>(new mmsFsiSolidTractionFvPatchVectorField(*this, iF)); }
    virtual void updateCoeffs();
    virtual void write(Ostream& os) const
    {
        solidTractionFvPatchVectorField::write(os);
        os.writeEntry("nominalTraction", nominal_);
    }
};
}
#endif
