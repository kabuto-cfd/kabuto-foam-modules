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

Class
    Foam::epsilonNonEquilibriumSpaldingWallFunctionFvPatchScalarField

Description
    y+-insensitive non-equilibrium wall function for epsilon and G, to go with
    nutNonEquilibriumSpaldingWallFunction. With uTau = ((nu + nut_w)|dU/dn|)^1/2
    from the nut wall function, y+ = uTau y/nu and Spalding's profile U+(y+):

    \verbatim
        gamma   = dU+/dy+                       (viscous share of the stress)
        tau+    = 1 + p+ y+,  p+ = nu (dp/ds)/uTau^3   (total stress / tau_w)
        G       = uTau^4/nu (tau+ - gamma) gamma
        epsilon = 2 nu k/y^2 exp(-Gamma) + Cmu^3/4 k^3/2/(kappa y) exp(-1/Gamma)
        Gamma   = 0.01 y*^4/(1 + 5 y*),  y* = Cmu^1/4 k^1/2 y/nu
    \endverbatim

    G is the production of the Spalding profile under the local total shear
    stress: zero at the wall, 0.25 uTau^4/nu at its peak in the buffer layer,
    tau_w uTau/(kappa y) in the log layer. epsilon blends the viscous and the
    logarithmic values with Kader's (1981) weights, as in the compound wall
    treatment of Popovac and Hanjalic (2007).

    Use the same wall function family on all walls of a case.

Usage
    \verbatim
    <patchName>
    {
        type            epsilonNonEquilibriumSpaldingWallFunction;
        value           $internalField;
    }
    \endverbatim
    Cmu, kappa and E are taken from the nut wall function.

SourceFiles
    epsilonNonEquilibriumSpaldingWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef epsilonNonEquilibriumSpaldingWallFunctionFvPatchScalarField_H
#define epsilonNonEquilibriumSpaldingWallFunctionFvPatchScalarField_H

#include "epsilonNonEquilibriumWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class epsilonNonEquilibriumSpaldingWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class epsilonNonEquilibriumSpaldingWallFunctionFvPatchScalarField
:
    public epsilonNonEquilibriumWallFunctionFvPatchScalarField
{
protected:

    // Protected Member Functions

        //- Calculate and return G and the wall-cell value
        virtual Pair<tmp<scalarField>> calculate
        (
            const momentumTransportModel&
        ) const;


public:

    //- Runtime type information
    TypeName("epsilonNonEquilibriumSpaldingWallFunction");


    // Constructors

        //- Construct from patch, internal field and dictionary
        epsilonNonEquilibriumSpaldingWallFunctionFvPatchScalarField
        (
            const fvPatch& p,
            const DimensionedField<scalar, fvMesh>& iF,
            const dictionary& dict
        )
        :
            epsilonNonEquilibriumWallFunctionFvPatchScalarField(p, iF, dict)
        {}

        //- Construct by mapping given onto a new patch
        epsilonNonEquilibriumSpaldingWallFunctionFvPatchScalarField
        (
            const epsilonNonEquilibriumSpaldingWallFunctionFvPatchScalarField& ptf,
            const fvPatch& p,
            const DimensionedField<scalar, fvMesh>& iF,
            const fieldMapper& mapper
        )
        :
            epsilonNonEquilibriumWallFunctionFvPatchScalarField(ptf, p, iF, mapper)
        {}

        //- Disallow copy without setting internal field reference
        epsilonNonEquilibriumSpaldingWallFunctionFvPatchScalarField
        (
            const epsilonNonEquilibriumSpaldingWallFunctionFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        epsilonNonEquilibriumSpaldingWallFunctionFvPatchScalarField
        (
            const epsilonNonEquilibriumSpaldingWallFunctionFvPatchScalarField& wfpsf,
            const DimensionedField<scalar, fvMesh>& iF
        )
        :
            epsilonNonEquilibriumWallFunctionFvPatchScalarField(wfpsf, iF)
        {}

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, fvMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new epsilonNonEquilibriumSpaldingWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member Operators

        using epsilonNonEquilibriumWallFunctionFvPatchScalarField::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
