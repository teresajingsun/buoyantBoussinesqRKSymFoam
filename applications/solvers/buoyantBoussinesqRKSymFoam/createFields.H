Info<< "Reading thermophysical properties\n" << endl;

Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
/*
Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
*/
Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
Info<< "createPhi\n" << endl;
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    reverseLinear<vector>(mesh).dotInterpolate(mesh.Sf(), U)
);

Info<< "Reading transportProperties\n" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar nu
(
    "nu",
    dimViscosity,
    transportProperties
);

surfaceScalarField nuf
(
    IOobject
    (
        "nuf",
        runTime.constant(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    nu
);

singlePhaseTransportModel laminarTransport(U, phi);
// Thermal expansion coefficient [1/K]
dimensionedScalar beta
(
    "beta",
    dimless/dimTemperature,
    laminarTransport
);

// Reference temperature [K]
dimensionedScalar TRef("TRef", dimTemperature, laminarTransport);

// Laminar Prandtl number
dimensionedScalar Pr("Pr", dimless, laminarTransport);

// Turbulent Prandtl number
dimensionedScalar Prt("Prt", dimless, laminarTransport);

// Turbulence model is required, to run DNS, set:
// 1. <transportModel Newtonian;> in constant/transportProperties, and
// 2. <simulationType laminar;> in constant/turbulenceProperties
Info<< "Creating turbulence model\n" << endl;
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

// Kinematic density for buoyancy force
volScalarField rhok
(
    IOobject
    (
        "rhok",
        runTime.timeName(),
        mesh
    ),
    1.0 - beta*(T - TRef)
);

// kinematic turbulent thermal thermal conductivity m2/s
Info<< "Reading field alphat\n" << endl;
volScalarField alphat
(
    IOobject
    (
        "alphat",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    //p_rgh + rhok*gh //cause issue "cannot be called for a calculatedFvPatchField"
    mesh
);

// Poisson solver solves only for the pressure correction
volScalarField pCorr
(
    IOobject
    (
        "pCorr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("pCorr", p.dimensions(), 0),
    p.boundaryField().types()
);

/*
volScalarField p_rghCorr
(
    IOobject
    (
        "p_rghCorr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("p_rghCorr", p.dimensions(), 0),
    p.boundaryField().types()
);
*/

// Extra pressure field to store pressure at the previous time step
volScalarField pn
(
    IOobject
    (
        "pn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("pn", p.dimensions(), 0),
    p.boundaryField().types()
);

autoPtr<RungeKuttaScheme> schemePtr
(
    RungeKuttaScheme::New(mesh, U, phi, pCorr)
);

RungeKuttaScheme& scheme = schemePtr();

// Read the coefficient to determine the weight of the previous time step
// in setting the pressure predictor
// pnPredCoef = 0 equals Chorin method, see references at the top of RKSymFoam.C
// pnPredCoef = 1 equals Van Kan method, see references at the top of RKSymFoam.C
const scalar pnPredCoef(readScalar(scheme.dict().lookup("pnPredCoef")));

Info<< "Running with predictor pressure = " << pnPredCoef << " * previous time step pressure" << endl;

label pRefCell = 0;
scalar pRefValue = 0.0;

setRefCell(p, scheme.dict(), pRefCell, pRefValue);

if (pCorr.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
}

mesh.setFluxRequired(pCorr.name());

#include "createIncompressibleRadiationModel.H"


