#include "mmsFsiSource.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "mmsFsiFields.H"

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(mmsFsiSource, 0);
    addToRunTimeSelectionTable(option, mmsFsiSource, dictionary);
}
}


Foam::fv::mmsFsiSource::mmsFsiSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::option(name, modelType, dict, mesh),
    region_(coeffs_.getOrDefault<word>("region", dict.get<word>("region"))),
    cachedTime_(-GREAT),
    cachedC_(),
    cachedS_()
{
    coeffs_.readEntry("fields", fieldNames_);
    fv::option::resetApplied();

    if (region_ != "fluid" && region_ != "solid")
    {
        FatalIOErrorInFunction(dict)
            << "region must be fluid or solid, not " << region_
            << exit(FatalIOError);
    }
    Info<< "    mmsFsiSource: " << region_ << " source on " << fieldNames_
        << endl;
}


const Foam::vectorField& Foam::fv::mmsFsiSource::source() const
{
    const scalar t = mesh_.time().value();
    const vectorField& C = mesh_.C().primitiveField();

    bool update = (t != cachedTime_) || (C.size() != cachedC_.size());
    if (!update)
    {
        forAll(C, i)
        {
            if (C[i] != cachedC_[i])
            {
                update = true;
                break;
            }
        }
    }

    if (update)
    {
        cachedTime_ = t;
        cachedC_ = C;
        cachedS_.setSize(C.size());
        if (region_ == "fluid")
        {
            forAll(C, i)
            {
                cachedS_[i] = mmsFsiFields::fluid(C[i], t).f;
            }
        }
        else
        {
            forAll(C, i)
            {
                cachedS_[i] = mmsFsiFields::solid(C[i], t).b0;
            }
        }
    }

    return cachedS_;
}


void Foam::fv::mmsFsiSource::addSup(fvMatrix<vector>& eqn, const label)
{
    // Body force on the right-hand side: eqn += S, with the dimensions
    // the equation needs (per unit volume of the equation)
    DimensionedField<vector, volMesh> S
    (
        IOobject("mmsFsiSource", mesh_.time().timeName(), mesh_),
        mesh_,
        dimensioned<vector>(eqn.dimensions()/dimVolume, Zero)
    );
    S.field() = source();
    eqn += S;
}


void Foam::fv::mmsFsiSource::addSup
(
    const volScalarField&,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    addSup(eqn, fieldi);
}


bool Foam::fv::mmsFsiSource::read(const dictionary& dict)
{
    return fv::option::read(dict);
}
