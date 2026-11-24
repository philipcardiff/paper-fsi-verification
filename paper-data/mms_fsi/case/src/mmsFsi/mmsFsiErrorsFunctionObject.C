#include "mmsFsiErrorsFunctionObject.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "volFields.H"
#include "fvcDdt.H"
#include "solidTractionFvPatchVectorField.H"
#include "mmsFsiFields.H"
#include "fluidModel.H"

namespace Foam
{
    defineTypeNameAndDebug(mmsFsiErrorsFunctionObject, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        mmsFsiErrorsFunctionObject,
        dictionary
    );
}


namespace
{
using namespace Foam;

// L1, L2, LInf of |e| with weights w
FixedList<scalar, 3> norms(const scalarField& e, const scalarField& w)
{
    const scalar sw = gSum(w);
    FixedList<scalar, 3> n;
    n[0] = gSum(w*e)/sw;
    n[1] = Foam::sqrt(gSum(w*sqr(e))/sw);
    n[2] = gMax(e);
    return n;
}
}


Foam::mmsFsiErrorsFunctionObject::mmsFsiErrorsFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name),
    time_(runTime),
    fluidRegion_("fluid"),
    solidRegion_("solid"),
    interfacePatch_("interface"),
    writeFields_(false)
{
    read(dict);
}


bool Foam::mmsFsiErrorsFunctionObject::read(const dictionary& dict)
{
    fluidRegion_ = dict.getOrDefault<word>("fluidRegion", "fluid");
    solidRegion_ = dict.getOrDefault<word>("solidRegion", "solid");
    interfacePatch_ = dict.getOrDefault<word>("interfacePatch", "interface");
    writeFields_ = dict.getOrDefault<bool>("writeFields", false);
    accelerationErrors_ = dict.getOrDefault<bool>("accelerationErrors", false);
    return true;
}


bool Foam::mmsFsiErrorsFunctionObject::compute()
{
    const scalar t = time_.value();

    if (!filePtr_ && Pstream::master())
    {
        const fileName dir(time_.globalPath()/"postProcessing"/"mmsFsiErrors");
        mkDir(dir);
        filePtr_.reset(new OFstream(dir/"errors.dat"));
        filePtr_() << "# time"
            << " U_L1 U_L2 U_Linf p_L1 p_L2 p_Linf pMean_L2"
            << " D_L1 D_L2 D_Linf Us_L1 Us_L2 Us_Linf"
            << " intD_L1 intD_L2 intD_Linf intTrac_L1 intTrac_L2 intTrac_Linf"
            << " intP_L1 intP_L2 intP_Linf"
            << " UWall_L1 UWall_L2 UWall_Linf pWall_L1 pWall_L2 pWall_Linf"
            << " tracN_L1 tracN_L2 tracN_Linf tracT_L1 tracT_L2 tracT_Linf"
            << " intU_L1 intU_L2 intU_Linf intA_L1 intA_L2 intA_Linf"
            << " intG_L1 intG_L2 intG_Linf" << nl;
        filePtr_().precision(10);
    }

    List<scalar> out;

    // ---- fluid ------------------------------------------------------------
    if (fluidRegion_ == "none")
    {
        for (label k = 0; k < 7; ++k) out.append(0);
    }
    else
    {
        const fvMesh& mesh = time_.lookupObject<fvMesh>(fluidRegion_);
        const volVectorField& U = mesh.lookupObject<volVectorField>("U");
        const volScalarField& p = mesh.lookupObject<volScalarField>("p");
        const vectorField& C = mesh.C().primitiveField();
        const scalarField& V = mesh.V().field();
        scalarField eU(C.size()), ep(C.size());
        forAll(C, i)
        {
            const mmsFsiFields::FluidState s = mmsFsiFields::fluid(C[i], t);
            eU[i] = mag(U[i] - s.U);
            ep[i] = p[i] - s.p;
        }
        const scalar meanDp = gSum(V*ep)/gSum(V);
        for (const scalar v : norms(eU, V)) out.append(v);
        for (const scalar v : norms(mag(ep), V)) out.append(v);
        out.append(Foam::sqrt(gSum(V*sqr(ep - meanDp))/gSum(V)));

        if (writeFields_ && time_.writeTime())
        {
            volScalarField UErr
            (
                IOobject("UError", time_.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
                mesh, dimensionedScalar(dimVelocity, Zero)
            );
            volScalarField pErr
            (
                IOobject("pError", time_.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
                mesh, dimensionedScalar(p.dimensions(), Zero)
            );
            UErr.primitiveFieldRef() = eU;
            pErr.primitiveFieldRef() = ep;
            UErr.write();
            pErr.write();
        }
    }

    // ---- solid ------------------------------------------------------------
    if (solidRegion_ == "none")
    {
        // Fluid-only case: no solid. The traction columns compare the
        // fluid's own interface traction (patchPressureForce and
        // patchViscousForce, as the coupling would pass them on) with the
        // exact fluid traction
        for (label i = 0; i < 9; ++i) out.append(0);

        const fvMesh& mesh = time_.lookupObject<fvMesh>(fluidRegion_);
        const fluidModel& fluid =
            dynamic_cast<const fluidModel&>
            (
                time_.lookupObject<IOdictionary>("fluidProperties")
            );
        const label patchI = mesh.boundaryMesh().findPatchID(interfacePatch_);
        const fvPatch& fp = mesh.boundary()[patchI];
        const scalar rho = mmsFsi::param_rho_f;
        const vectorField nf(fp.nf());
        const vectorField tr
        (
            (-fluid.patchPressureForce(patchI)*nf + fluid.patchViscousForce(patchI))/rho
        );
        const scalar nu = mmsFsi::param_mu_f/mmsFsi::param_rho_f;
        scalarField eTr(nf.size());
        eTrN_.setSize(nf.size());
        eTrT_.setSize(nf.size());
        trMagSf_ = fp.magSf();
        forAll(nf, i)
        {
            const mmsFsiFields::FluidState s = mmsFsiFields::fluid(fp.Cf()[i], t);
            const vector& n = nf[i];
            const vector ex(-s.p*n + nu*(n & (s.gradU + s.gradU.T())));
            const vector e(tr[i] - ex);
            eTr[i] = mag(e);
            eTrN_[i] = mag(e & n);
            eTrT_[i] = mag(e - (e & n)*n);
        }
        for (const scalar v : norms(eTr, trMagSf_)) out.append(v);
    }
    else
    {
        const fvMesh& mesh = time_.lookupObject<fvMesh>(solidRegion_);
        const volVectorField& D = mesh.lookupObject<volVectorField>("D");
        const vectorField& C = mesh.C().primitiveField();
        const scalarField& V = mesh.V().field();
        scalarField eD(C.size()), eUs(C.size(), Zero);
        const volVectorField* UsPtr =
            mesh.foundObject<volVectorField>("U")
          ? &mesh.lookupObject<volVectorField>("U") : nullptr;
        forAll(C, i)
        {
            const mmsFsiFields::SolidState s = mmsFsiFields::solid(C[i], t);
            eD[i] = mag(D[i] - s.D);
            if (UsPtr)
            {
                eUs[i] = mag((*UsPtr)[i] - s.U);
            }
        }
        if (!UsPtr)
        {
            // Displacement solids: velocity from the BDF2 ddt of D is not
            // stored; use ddt(D) with the case's ddt scheme
            const volVectorField Us(fvc::ddt(D));
            forAll(C, i)
            {
                eUs[i] = mag(Us[i] - mmsFsiFields::solid(C[i], t).U);
            }
        }
        for (const scalar v : norms(eD, V)) out.append(v);
        for (const scalar v : norms(eUs, V)) out.append(v);

        // interface
        const label patchI = mesh.boundaryMesh().findPatchID(interfacePatch_);
        const fvPatch& patch = mesh.boundary()[patchI];
        const vectorField& Cf = patch.Cf();
        const scalarField& magSf = patch.magSf();
        const fvPatchVectorField& Dp = D.boundaryField()[patchI];
        scalarField eDi(Cf.size()), eTr(Cf.size());
        eTrN_.setSize(Cf.size());
        eTrT_.setSize(Cf.size());
        trMagSf_ = magSf;
        const solidTractionFvPatchVectorField* stPtr =
            isA<solidTractionFvPatchVectorField>(Dp)
          ? &refCast<const solidTractionFvPatchVectorField>(Dp) : nullptr;
        forAll(Cf, i)
        {
            const mmsFsiFields::SolidState s = mmsFsiFields::solid(Cf[i], t);
            eDi[i] = mag(Dp[i] - s.D);
            vector xG, tr, nf;
            mmsFsiFields::interface(Cf[i].x(), t, xG, tr, nf);
            if (stPtr)
            {
                // applied traction on the solid: sigma_s . n_solid = -tr
                // (nf is the fluid outward normal = - solid outward normal)
                const vector applied =
                    stPtr->traction()[i] + stPtr->pressure()[i]*nf;
                eTr[i] = mag(applied + tr);
                const vector e(applied + tr);
                eTrN_[i] = mag(e & nf);
                eTrT_[i] = mag(e - (e & nf)*nf);
            }
            else
            {
                eTr[i] = 0;
                eTrN_[i] = 0;
                eTrT_[i] = 0;
            }
        }
        for (const scalar v : norms(eDi, magSf)) out.append(v);
        for (const scalar v : norms(eTr, magSf)) out.append(v);
    }

    // ---- fluid interface pressure -------------------------------------------
    if (fluidRegion_ == "none")
    {
        // Solid only: the traction columns, the fluid ones zero
        for (label k = 0; k < 9; ++k) out.append(0);
        for (const scalar v : norms(eTrN_, trMagSf_)) out.append(v);
        for (const scalar v : norms(eTrT_, trMagSf_)) out.append(v);
        for (label k = 0; k < 9; ++k) out.append(0);
    }
    else
    {
        const fvMesh& mesh = time_.lookupObject<fvMesh>(fluidRegion_);
        const volScalarField& p = mesh.lookupObject<volScalarField>("p");
        const label patchI = mesh.boundaryMesh().findPatchID(interfacePatch_);
        const vectorField& Cf = mesh.boundary()[patchI].Cf();
        const scalarField& magSf = mesh.boundary()[patchI].magSf();
        scalarField ep(Cf.size());
        forAll(Cf, i)
        {
            ep[i] = mag(p.boundaryField()[patchI][i] - mmsFsiFields::fluid(Cf[i], t).p);
        }
        for (const scalar v : norms(ep, magSf)) out.append(v);

        // errors in the cells next to the interface (the wall layer)
        const volVectorField& U = mesh.lookupObject<volVectorField>("U");
        const labelUList& fc = mesh.boundary()[patchI].faceCells();
        const scalarField& V = mesh.V().field();
        scalarField eUw(fc.size()), epw(fc.size()), Vw(fc.size());
        const vectorField& C = mesh.C().primitiveField();
        forAll(fc, i)
        {
            const mmsFsiFields::FluidState s = mmsFsiFields::fluid(C[fc[i]], t);
            eUw[i] = mag(U[fc[i]] - s.U);
            epw[i] = mag(p[fc[i]] - s.p);
            Vw[i] = V[fc[i]];
        }
        for (const scalar v : norms(eUw, Vw)) out.append(v);
        for (const scalar v : norms(epw, Vw)) out.append(v);
        if (trMagSf_.size() || solidRegion_ != "none")
        {
            for (const scalar v : norms(eTrN_, trMagSf_)) out.append(v);
            for (const scalar v : norms(eTrT_, trMagSf_)) out.append(v);
        }
        else
        {
            for (label i = 0; i < 6; ++i) out.append(0);
        }

        // the fluid interface velocity against the exact one at the face
        // centres
        scalarField eUi(Cf.size());
        forAll(Cf, i)
        {
            eUi[i] = mag(U.boundaryField()[patchI][i] - mmsFsiFields::fluid(Cf[i], t).U);
        }
        for (const scalar v : norms(eUi, magSf)) out.append(v);

        // the wall acceleration (boundary ddt(U)) against the exact material
        // acceleration, and the interface pressure gradient against the
        // exact normal pressure gradient, at the face centres
        // (optional: fvc::ddt creates old-time mesh data, which changes the
        // solution at round-off level)
        const bool accel = accelerationErrors_;
        const volVectorField ddtU(accel ? volVectorField(fvc::ddt(U)) : volVectorField("zeroDdtU", 0*U));
        const vectorField nfi(mesh.boundary()[patchI].nf());
        const fvPatchScalarField& pb = p.boundaryField()[patchI];
        scalarField eA(Cf.size()), eG(Cf.size());
        forAll(Cf, i)
        {
            const scalar ht = 1e-4;
            const mmsFsiFields::FluidState s = mmsFsiFields::fluid(Cf[i], t);
            const vector a
            (
                (mmsFsiFields::fluid(Cf[i], t + ht).U - mmsFsiFields::fluid(Cf[i], t - ht).U)/(2*ht)
              + (s.U & s.gradU)
            );
            eA[i] = mag(ddtU.boundaryField()[patchI][i] - a);
            eG[i] = mag(pb.snGrad()()[i] - (s.gradp & nfi[i]));
        }
        for (const scalar v : norms(eA, magSf)) out.append(v);
        for (const scalar v : norms(eG, magSf)) out.append(v);
    }

    if (Pstream::master())
    {
        filePtr_() << t;
        for (const scalar v : out) filePtr_() << ' ' << v;
        filePtr_() << endl;
    }
    Info<< "mmsFsiErrors: t = " << t << " U_L2 " << out[1] << " p_L2 " << out[4]
        << " D_L2 " << out[8] << " intTrac_L2 " << out[17] << endl;

    return true;
}


bool Foam::mmsFsiErrorsFunctionObject::execute()
{
    return compute();
}


bool Foam::mmsFsiErrorsFunctionObject::write()
{
    return true;
}
