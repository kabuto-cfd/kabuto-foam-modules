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

    The structure follows fv::limitMag of OpenFOAM
    (Copyright (C) 2016-2026 OpenFOAM Foundation).

Class
    Foam::fv::twoLayerEpsilon

Description
    Wolfstein's two-layer treatment for k-epsilon models in the viscosity-
    affected near-wall region, as in Fluent's enhanced wall treatment, so
    that a k-epsilon model runs on meshes from y+ ~ 1 to the log layer.
    Together with nutNonEquilibriumSpaldingWallFunction and
    epsilonNonEquilibriumSpaldingWallFunction this makes the near-wall
    treatment of epsilon models y+-insensitive.

    With the wall-distance Reynolds number Re_y = y k^1/2/nu, in the cells
    near walls:

    \verbatim
        epsilon_2l = k^3/2/l_eps,   l_eps = y Cl (1 - exp(-Re_y/(2 Cl)))
        nut_2l     = Cmu l_mu k^1/2, l_mu  = y Cl (1 - exp(-Re_y/Amu))
        Cl         = kappa Cmu^-3/4
        lambda     = (1 + tanh((Re_y - ReyStar)/A))/2,
        A          = deltaReyStar/atanh(0.98)
    \endverbatim

    After the epsilon equation is solved, epsilon is blended to
    lambda epsilon + (1 - lambda) epsilon_2l (the two-layer value where
    Re_y << ReyStar, the transported one where Re_y >> ReyStar), and nut to
    lambda nut + (1 - lambda) nut_2l after the model computes it. The cells
    next to walls keep the epsilon the wall functions give them.

    Reference:
    \verbatim
        Wolfshtein, M. (1969). The velocity and temperature distribution of
        one-dimensional flow with turbulence augmentation and pressure
        gradient. Int. J. Heat Mass Transfer, 12, 301-318.
    \endverbatim

Usage
    system/fvConstraints:
    \verbatim
    twoLayerEpsilon
    {
        type        twoLayerEpsilon;
        libs        ("libkabutoWallFunctions.so");

        // Optional, defaults shown
        epsilon         epsilon;
        nut             nut;
        k               k;
        Cmu             0.09;
        kappa           0.41;
        Amu             70;
        ReyStar         200;
        deltaReyStar    20;
    }
    \endverbatim
    The wall distance needs a wallDist entry in system/fvSchemes
    (e.g. wallDist { method meshWave; }).

SourceFiles
    twoLayerEpsilon.C

\*---------------------------------------------------------------------------*/

#ifndef twoLayerEpsilon_H
#define twoLayerEpsilon_H

#include "fvConstraint.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class twoLayerEpsilon Declaration
\*---------------------------------------------------------------------------*/

class twoLayerEpsilon
:
    public fvConstraint
{
    // Private Data

        word epsilonName_;
        word nutName_;
        word kName_;
        scalar Cmu_;
        scalar kappa_;
        scalar Amu_;
        scalar ReyStar_;
        scalar deltaReyStar_;


    // Private Member Functions

        //- Read the coefficients
        void readCoeffs(const dictionary& dict);

        //- Wall-distance Reynolds number and blending factor per cell
        void blending(scalarField& Rey, scalarField& lambda) const;


public:

    //- Runtime type information
    TypeName("twoLayerEpsilon");


    // Constructors

        //- Construct from components
        twoLayerEpsilon
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        twoLayerEpsilon(const twoLayerEpsilon&) = delete;


    //- Destructor
    virtual ~twoLayerEpsilon()
    {}


    // Member Functions

        //- Return the list of constrained fields (epsilon and nut)
        virtual wordList constrainedFields() const;

        using fvConstraint::constrain;

        //- Blend epsilon (after its equation is solved) or nut (after the
        //  model computes it) with the two-layer value near walls
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
        void operator=(const twoLayerEpsilon&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
