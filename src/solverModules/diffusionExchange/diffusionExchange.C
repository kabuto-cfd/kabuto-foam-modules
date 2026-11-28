/*---------------------------------------------------------------------------*\
  Kabuto CFD modules for OpenFOAM-14
  Copyright (C) 2026 Kabuto CFD
-------------------------------------------------------------------------------
License
    This file is part of kabuto-foam-modules, add-on libraries for OpenFOAM.

    It is free software: you can redistribute it and/or modify it under the
    terms of the GNU General Public License as published by the Free Software
    Foundation, either version 3 of the License, or (at your option) any later
    version.

    It is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
    details.

    You should have received a copy of the GNU General Public License along
    with it.  If not, see <http://www.gnu.org/licenses/>.

    The structure follows the solver module solvers::solid of OpenFOAM
    (Copyright (C) 2022-2026 OpenFOAM Foundation).

\*---------------------------------------------------------------------------*/

#include "diffusionExchange.H"
#include "fvmDdt.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "fvmDiv.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(diffusionExchange, 0);
    addToRunTimeSelectionTable(solver, diffusionExchange, fvMesh);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::solvers::diffusionExchange::readOrUniform
(
    const word& name,
    const dimensionSet& dims,
    const bool required
) const
{
    typeIOobject<volScalarField> io
    (
        name,
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.headerOk())
    {
        Info<< "Reading field " << name << endl;

        return tmp<volScalarField>(new volScalarField(io, mesh));
    }

    const dimensionedScalar value
    (
        name,
        dims,
        required
      ? physicalProperties_.lookup<scalar>(name)
      : physicalProperties_.lookupOrDefault<scalar>(name, 0)
    );

    Info<< "No field " << name << ", using the uniform value "
        << value.value() << endl;

    return volScalarField::New(name, mesh, value);
}


bool Foam::solvers::diffusionExchange::dependenciesModified() const
{
    return runTime.controlDict().modified();
}


bool Foam::solvers::diffusionExchange::read()
{
    solver::read();

    maxDeltaT_ =
        runTime.controlDict().found("maxDeltaT")
      ? runTime.controlDict().lookup<scalar>("maxDeltaT", runTime.userUnits())
      : vGreat;

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::diffusionExchange::diffusionExchange(fvMesh& mesh)
:
    solver(mesh),

    maxDeltaT_(vGreat),

    physicalProperties_
    (
        IOobject
        (
            "physicalProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    C_
    (
        IOobject
        (
            physicalProperties_.lookupOrDefault<word>("field", "C"),
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    D_
    (
        IOobject
        (
            "D",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        readOrUniform("D", dimArea/dimTime, true)
    ),

    kExchange_
    (
        IOobject
        (
            "kExchange",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        readOrUniform("kExchange", dimless/dimTime, false)
    ),

    source_
    (
        IOobject
        (
            "source",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        readOrUniform("source", C_.dimensions()/dimTime, false)
    ),
    phi_(nullptr),
    Cref_
    (
        Function1<scalar>::New
        (
            "Cref",
            runTime.userUnits(),
            C_.dimensions(),
            physicalProperties_
        )
    ),

    C(C_),
    D(D_),
    kExchange(kExchange_)
{
    if (LTS)
    {
        FatalError
            << type()
            << " solver does not support LTS, use 'steadyState' ddtScheme"
            << exit(FatalError);
    }

    // Frozen flux of a flow computed beforehand, if there is one
    typeIOobject<surfaceScalarField> phiIO
    (
        "phi",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (phiIO.headerOk())
    {
        Info<< "Reading the frozen flux phi" << endl;
        phi_.set(new surfaceScalarField(phiIO, mesh));
    }
    else
    {
        Info<< "No field phi: no convection" << endl;
    }

    // Floor of the diffusivity, for the harmonic interpolation to the faces
    const dimensionedScalar Dmin
    (
        "Dmin",
        D_.dimensions(),
        physicalProperties_.lookupOrDefault<scalar>("Dmin", 1e-30)
    );

    label nClosed = 0;
    forAll(D_, celli)
    {
        if (D_[celli] <= Dmin.value())
        {
            nClosed++;
        }
    }
    reduce(nClosed, sumOp<label>());

    D_ = max(D_, Dmin);

    Info<< type() << ": " << nClosed << " of "
        << returnReduce(mesh.nCells(), sumOp<label>())
        << " cells are closed (D <= Dmin = " << Dmin.value() << ")" << nl
        << endl;

    // Read the controls
    read();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::diffusionExchange::~diffusionExchange()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::solvers::diffusionExchange::maxDeltaT() const
{
    return min(fvModels().maxDeltaT(), maxDeltaT_);
}


void Foam::solvers::diffusionExchange::preSolve()
{
    fvModels().preUpdateMesh();

    // Update the mesh for topology change, mesh to mesh mapping
    mesh_.update();
}


void Foam::solvers::diffusionExchange::moveMesh()
{}


void Foam::solvers::diffusionExchange::motionCorrector()
{}


void Foam::solvers::diffusionExchange::prePredictor()
{}


void Foam::solvers::diffusionExchange::momentumTransportPredictor()
{}


void Foam::solvers::diffusionExchange::thermophysicalTransportPredictor()
{}


void Foam::solvers::diffusionExchange::momentumPredictor()
{}


void Foam::solvers::diffusionExchange::thermophysicalPredictor()
{
    // Value of the reservoir at the middle of the time step
    const dimensionedScalar Cref
    (
        "Cref",
        C_.dimensions(),
        Cref_->value(runTime.value() - 0.5*runTime.deltaTValue())
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix CEqn
        (
            fvm::ddt(C_)
          - fvm::laplacian(D_, C_)
          + fvm::Sp(kExchange_, C_)
         ==
            kExchange_*Cref
          + source_
          + fvModels().source(C_)
        );

        if (phi_.valid())
        {
            CEqn += fvm::div(phi_(), C_);
        }

        CEqn.relax();

        fvConstraints().constrain(CEqn);

        CEqn.solve();

        fvConstraints().constrain(C_);
    }
}


void Foam::solvers::diffusionExchange::pressureCorrector()
{}


void Foam::solvers::diffusionExchange::momentumTransportCorrector()
{}


void Foam::solvers::diffusionExchange::thermophysicalTransportCorrector()
{}


void Foam::solvers::diffusionExchange::postSolve()
{}


// ************************************************************************* //
