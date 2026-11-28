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

Class
    Foam::solvers::diffusionExchange

Description
    Solver module for the diffusion of a scalar through a medium whose
    diffusivity varies from cell to cell, with a linear exchange between the
    medium and a reservoir of prescribed, time-dependent value:

        ddt(C) - div(D grad(C)) = kExchange (Cref(t) - C) + fvModels

    Typical use: transport of a solute or of particles through a porous or
    fibrous medium fed from a network that the mesh does not resolve (the
    exchange term), where the effective diffusivity of each cell comes from a
    sub-grid model such as a pore-network model and is zero in the cells the
    particles cannot cross.

    D and kExchange are read as fields from the start time directory; where a
    field file is absent the uniform value of constant/physicalProperties is
    used. Cells with D = 0 are impermeable when D is interpolated to the faces
    with the harmonic scheme (a face is closed as soon as one of its two cells
    is closed, and a jump of D between two cells is handled as two resistances
    in series):

        laplacianSchemes { laplacian(D,C) Gauss harmonic corrected; }

    For that scheme D is floored at Dmin (default 1e-30 m^2/s), which keeps the
    reciprocal finite and leaves the closed cells closed for every practical
    purpose. With kExchange = 0 as well, such a cell keeps its initial value.

    Cref is evaluated at the middle of the time step.

Usage
    system/controlDict:
    \verbatim
    solver          diffusionExchange;
    \endverbatim
    (foamRun loads libdiffusionExchangeSolver.so by that name; no libs entry.)

    constant/physicalProperties:
    \table
        Property  | Description                                  | Req'd | Default
        Cref      | Function1 of time, value of the reservoir    | yes   |
        D         | uniform diffusivity [m^2/s], if no field D   | if no field |
        kExchange | uniform exchange rate [1/s], if no field     | no    | 0
        Dmin      | floor of the diffusivity [m^2/s]             | no    | 1e-30
        field     | name of the solved field                     | no    | C
    \endtable

    Example:
    \verbatim
    Cref
    {
        type    table;
        values  ((0 1) (21600 0.5) (43200 0.25) (86400 0.0625));
    }
    \endverbatim

    Optional maxDeltaT in system/controlDict limits an adjustable time step.

SourceFiles
    diffusionExchange.C

\*---------------------------------------------------------------------------*/

#ifndef diffusionExchange_H
#define diffusionExchange_H

#include "solver.H"
#include "volFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                      Class diffusionExchange Declaration
\*---------------------------------------------------------------------------*/

class diffusionExchange
:
    public solver
{

protected:

    // Control parameters

        scalar maxDeltaT_;


    // Properties

        //- constant/physicalProperties
        IOdictionary physicalProperties_;

        //- Solved field
        volScalarField C_;

        //- Diffusivity [m^2/s]
        volScalarField D_;

        //- Exchange rate with the reservoir [1/s]
        volScalarField kExchange_;

        //- Value of the reservoir as a function of time
        autoPtr<Function1<scalar>> Cref_;


    // Protected Member Functions

        //- Read the field of the start time directory if it is there,
        //  otherwise the uniform value (or the default) from the dictionary
        tmp<volScalarField> readOrUniform
        (
            const word& name,
            const dimensionSet& dims,
            const bool required
        ) const;

        //- Return true if the solver's dependencies have been modified
        virtual bool dependenciesModified() const;

        //- Read controls
        virtual bool read();


public:

    // Public Data

        //- Reference to the solved field
        const volScalarField& C;

        //- Reference to the diffusivity
        const volScalarField& D;

        //- Reference to the exchange rate
        const volScalarField& kExchange;


    //- Runtime type information
    TypeName("diffusionExchange");


    // Constructors

        //- Construct from region mesh
        diffusionExchange(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        diffusionExchange(const diffusionExchange&) = delete;


    //- Destructor
    virtual ~diffusionExchange();


    // Member Functions

        //- Return the current maximum time-step for stable solution
        virtual scalar maxDeltaT() const;

        //- Called at the start of the time-step, before the PIMPLE loop
        virtual void preSolve();

        //- Called at the start of the PIMPLE loop to move the mesh
        virtual void moveMesh();

        //- Corrections that follow mesh motion
        virtual void motionCorrector();

        //- Called at the beginning of the PIMPLE loop
        virtual void prePredictor();

        //- Predict the momentum transport
        virtual void momentumTransportPredictor();

        //- Predict thermophysical transport
        virtual void thermophysicalTransportPredictor();

        //- Construct and optionally solve the momentum equation
        virtual void momentumPredictor();

        //- Construct and solve the diffusion-exchange equation
        virtual void thermophysicalPredictor();

        //- Construct and solve the pressure equation in the PISO loop
        virtual void pressureCorrector();

        //- Correct the momentum transport
        virtual void momentumTransportCorrector();

        //- Correct the thermophysical transport
        virtual void thermophysicalTransportCorrector();

        //- Called after the PIMPLE loop at the end of the time-step
        virtual void postSolve();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const diffusionExchange&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
