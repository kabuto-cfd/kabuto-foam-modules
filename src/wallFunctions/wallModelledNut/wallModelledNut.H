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

    The structure follows twoLayerEpsilon of this library.

Class
    Foam::fv::wallModelledNut

Description
    Applies the eddy-viscosity floor of nutSpaldingLESWallFunction patches
    with nutFloor on: after the subgrid model has set nut (and its boundary
    values), the cells of each such patch's wall-modelled layer get at least
    the eddy viscosity of Spalding's law with the face's friction velocity
    (see nutSpaldingLESWallFunctionFvPatchScalarField.H). The boundary values
    are then re-evaluated so that processor patches see the raised cells.

    The subgrid models of OpenFOAM and of this repository call the
    constraints of nut right after computing it (LESeddyViscosity models'
    correctNut).

Usage
    system/fvConstraints:
    \verbatim
    wallModelledNut
    {
        type        wallModelledNut;
        libs        ("libkabutoWallFunctions.so");
        // Optional, default shown
        nut         nut;
    }
    \endverbatim

SourceFiles
    wallModelledNut.C

\*---------------------------------------------------------------------------*/

#ifndef wallModelledNut_H
#define wallModelledNut_H

#include "fvConstraint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class wallModelledNut Declaration
\*---------------------------------------------------------------------------*/

class wallModelledNut
:
    public fvConstraint
{
    // Private Data

        //- Name of the eddy-viscosity field
        word nutName_;


    // Private Member Functions

        //- Read the coefficients
        void readCoeffs(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("wallModelledNut");


    // Constructors

        //- Construct from components
        wallModelledNut
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        wallModelledNut(const wallModelledNut&) = delete;


    //- Destructor
    virtual ~wallModelledNut()
    {}


    // Member Functions

        //- Return the list of constrained fields (nut)
        virtual wordList constrainedFields() const;

        using fvConstraint::constrain;

        //- Raise nut in the wall-modelled layers
        virtual bool constrain(volScalarField& field) const;

        //- Update for mesh motion
        virtual bool movePoints()
        {
            return true;
        }

        //- Update topology using the given map
        virtual void topoChange(const polyTopoChangeMap&)
        {}

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&)
        {}

        //- Redistribute or update using the given distribution map
        virtual void distribute(const polyDistributionMap&)
        {}

        //- Read dictionary
        virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const wallModelledNut&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
