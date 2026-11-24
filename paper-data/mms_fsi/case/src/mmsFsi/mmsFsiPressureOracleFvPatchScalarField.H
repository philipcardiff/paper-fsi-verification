/*---------------------------------------------------------------------------*\
    mmsFsiPressureOracle: a fixedGradient pressure whose gradient is set by
    newtonIcoFluid's wall-normal momentum balance (as for a plain
    fixedGradient), but whose boundary VALUE is overwritten with the exact
    manufactured pressure after every evaluation. Diagnostic only: it shows
    what an exact interface pressure value does to the rest of the solution.
    Optional: "mode exact" (default) or "mode correction", which instead adds
    k & grad(p_exact) at the face (the exact non-orthogonal correction) to the
    uncorrected value.
\*---------------------------------------------------------------------------*/
#ifndef mmsFsiPressureOracleFvPatchScalarField_H
#define mmsFsiPressureOracleFvPatchScalarField_H

#include "fixedGradientFvPatchFields.H"

namespace Foam
{

class mmsFsiPressureOracleFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    word mode_;

public:
    TypeName("mmsFsiPressureOracle");

    mmsFsiPressureOracleFvPatchScalarField
    (const fvPatch& p, const DimensionedField<scalar, volMesh>& iF)
    : fixedGradientFvPatchScalarField(p, iF), mode_("exact") {}

    mmsFsiPressureOracleFvPatchScalarField
    (const fvPatch& p, const DimensionedField<scalar, volMesh>& iF, const dictionary& dict)
    : fixedGradientFvPatchScalarField(p, iF),
      mode_(dict.getOrDefault<word>("mode", "exact"))
    {
        gradient() = Zero;
        fvPatchScalarField::operator=(patchInternalField());
    }

    mmsFsiPressureOracleFvPatchScalarField
    (
        const mmsFsiPressureOracleFvPatchScalarField& ptf,
        const fvPatch& p, const DimensionedField<scalar, volMesh>& iF,
        const fvPatchFieldMapper& m
    )
    : fixedGradientFvPatchScalarField(ptf, p, iF, m), mode_(ptf.mode_) {}

    mmsFsiPressureOracleFvPatchScalarField(const mmsFsiPressureOracleFvPatchScalarField& ptf)
    : fixedGradientFvPatchScalarField(ptf), mode_(ptf.mode_) {}

    mmsFsiPressureOracleFvPatchScalarField
    (const mmsFsiPressureOracleFvPatchScalarField& ptf, const DimensionedField<scalar, volMesh>& iF)
    : fixedGradientFvPatchScalarField(ptf, iF), mode_(ptf.mode_) {}

    virtual tmp<fvPatchField<scalar>> clone() const
    { return tmp<fvPatchField<scalar>>(new mmsFsiPressureOracleFvPatchScalarField(*this)); }

    virtual tmp<fvPatchField<scalar>> clone(const DimensionedField<scalar, volMesh>& iF) const
    { return tmp<fvPatchField<scalar>>(new mmsFsiPressureOracleFvPatchScalarField(*this, iF)); }

    virtual void evaluate(const Pstream::commsTypes commsType = Pstream::commsTypes::buffered);

    virtual void write(Ostream& os) const
    {
        fixedGradientFvPatchScalarField::write(os);
        os.writeEntry("mode", mode_);
    }
};

} // namespace Foam

#endif
