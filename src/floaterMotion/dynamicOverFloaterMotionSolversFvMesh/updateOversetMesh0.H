/*---------------------------------------------------------------------------*\
    Module Name:       FloatStepper
    Description:       OpenFOAM extension module for fluid-rigid body coupling
    License:           GNU General Public License (GPL) version 3
    Copyright:         2023 Johan Roenby, STROMNING APS
\*---------------------------------------------------------------------------*/
const immiscibleIncompressibleTwoPhaseMixture& mixture
    = lookupObject<immiscibleIncompressibleTwoPhaseMixture>("transportProperties");

const volScalarField& alpha1(mixture.alpha1());
const volScalarField& alpha2(mixture.alpha2());
//const volScalarField& alpha1 = lookupObject<volScalarField>("alpha.water");
//const volScalarField& alpha2 = lookupObject<volScalarField>("alpha.air");
const volScalarField& rho = lookupObject<volScalarField>("rho");
const surfaceScalarField& rhoPhi = lookupObject<surfaceScalarField>("rhoPhi");
const volVectorField& U = lookupObject<volVectorField>("U");
const surfaceScalarField& phif = lookupObject<surfaceScalarField>("phi");
const volScalarField& p_rgh = lookupObject<volScalarField>("p_rgh");
const volScalarField& p = lookupObject<volScalarField>("p");
const surfaceVectorField& Uf = lookupObject<surfaceVectorField>("Uf");
const volScalarField& rAU = lookupObject<volScalarField>("rAU");

// find overset fields and scalars
const volScalarField& cellMask = lookupObject<volScalarField>("cellMask");
const surfaceScalarField& faceMask = lookupObject<surfaceScalarField>("faceMask");
const volScalarField& interpolatedCells = lookupObject<volScalarField>("interpolatedCells");
const tmp<surfaceScalarField>& talphaPhi1Corr0 = lookupObject<tmp<surfaceScalarField>>("talphaPhi1Corr");
// don't like this todo: see if I can change it
const dynamicFvMesh& mesh = lookupObject<dynamicFvMesh>("mesh");

volScalarField cellMask0(cellMask);
surfaceScalarField faceMask0(faceMask);
volScalarField interpolatedCells0(interpolatedCells);
tmp<surfaceScalarField> talphaPhi1Corr00(talphaPhi1Corr0);

volScalarField alpha0(alpha1);
volScalarField alpha20(alpha2);
volVectorField U0("U0", U);
surfaceScalarField phi0(phif);
isoAdvection advector0(alpha0, phi0, U0);
surfaceScalarField rhoPhi0(rhoPhi);
volScalarField rho0(rho);
volScalarField p_rgh0(p_rgh);
volScalarField p0("p0", p);
surfaceVectorField Uf0(Uf);//(Uf());
volScalarField rAU0(rAU);//(rAU());

bodySolver.solve();
disp += bodySolver.pointDisplacement().primitiveField();
fvMesh::movePoints(bodySolver.points0() + disp);
disp -= bodySolver.pointDisplacement().primitiveField();
//movePoints(bodySolver.newPoints());

U0.correctBoundaryConditions();

oversetFvMeshBase::update();
if (polyMesh::changing())
{
    // Do not apply previous time-step mesh compression flux
    // if the mesh topology changed
    if (polyMesh::topoChanging())
    {
        talphaPhi1Corr00.clear();
    }

    // Update cellMask field for blocking out hole cells
    #include "setCellMask0.H"
    #include "setInterpolatedCells0.H"
    #include "correctPhiFaceMask0.H"
}

// At this stage, the mesh point have been moved but all the fields are from the
// old mesh

uniformDimensionedScalarField hRef
(
    IOobject
    (
        "hRef",
        time().constant(),
        *this,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    dimensionedScalar(dimLength, Zero)
);

dimensionedVector g("g", dimAcceleration, Zero);
g = rho.time().lookupObject<uniformDimensionedVectorField>("g");

dimensionedScalar ghRef
(
    mag(g.value()) > SMALL
    ? g & (cmptMag(g.value())/mag(g.value()))*hRef
    : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0)
);
volScalarField gh("gh", (g & C()) - ghRef);
surfaceScalarField ghf("ghf", (g & Cf()) - ghRef);

//MRF.update();

pimpleControl& pimple =
    const_cast<pimpleControl&>
    (
        lookupObject<pimpleControl>("solutionControl")
    );

bool correctPhi
(
    pimple.dict().getOrDefault("correctPhi", dynamic())
);

if (correctPhi)
{
    // Calculate absolute flux
    // from the mapped surface velocity
    phi0 = Sf() & Uf0;

    #include "correctPhi0.H"

    // Make the flux relative to the mesh motion
    fvc::makeRelative(phi0, U0);

//    mixture.correct();
}
/*
if (checkMeshCourantNo)
{
    #include "meshCourantNo.H"
}
*/