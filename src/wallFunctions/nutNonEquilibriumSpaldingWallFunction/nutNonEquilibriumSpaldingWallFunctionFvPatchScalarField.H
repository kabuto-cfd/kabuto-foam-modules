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

    The structure follows nutUSpaldingWallFunctionFvPatchScalarField of OpenFOAM
    (Copyright (C) 2011-2026 OpenFOAM Foundation).

Class
    Foam::nutNonEquilibriumSpaldingWallFunctionFvPatchScalarField

Description
    y+-insensitive non-equilibrium wall function for the turbulent viscosity:
    Spalding's law, valid from the viscous sublayer through the log layer,
    applied to Kim and Choudhury's pressure-gradient-sensitised velocity
    UTilde (wallLaw.H). The friction velocity solves

    \verbatim
        y uTau/nu = y+(UTilde/uTau)                (Spalding's law)
        nut_w     = uTau^2 y/U - nu
    \endverbatim

    so the first cell may sit anywhere from y+ ~ 1 to the log layer. Without
    a pressure gradient this is nutUSpaldingWallFunction. The wall shear never
    reverses: nut_w >= -0.999 nu (UTilde <= 0 gives uTau = 0).

    Use it with epsilonNonEquilibriumSpaldingWallFunction or
    omegaNonEquilibriumSpaldingWallFunction.

    References:
    \verbatim
        Spalding, D. B. (1961). A single formula for the law of the wall.
        Journal of Applied Mechanics, 28(3), 455-458.

        Kim, S.-E., Choudhury, D. (1995).
        A near-wall treatment using wall functions sensitized to pressure
        gradient. ASME FED Vol. 217, Separated and Complex Flows.
    \endverbatim

Usage
    As nutNonEquilibriumWallFunction:
    \verbatim
    <patchName>
    {
        type            nutNonEquilibriumSpaldingWallFunction;
        value           uniform 0;
    }
    \endverbatim

SourceFiles
    nutNonEquilibriumSpaldingWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nutNonEquilibriumSpaldingWallFunctionFvPatchScalarField_H
#define nutNonEquilibriumSpaldingWallFunctionFvPatchScalarField_H

#include "nutNonEquilibriumWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
    Class nutNonEquilibriumSpaldingWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutNonEquilibriumSpaldingWallFunctionFvPatchScalarField
:
    public nutNonEquilibriumWallFunctionFvPatchScalarField
{
protected:

    // Protected Member Functions

        //- Calculate the turbulent viscosity
        virtual tmp<scalarField> nut() const;


public:

    //- Runtime type information
    TypeName("nutNonEquilibriumSpaldingWallFunction");


    // Constructors

        //- Construct from patch, internal field and dictionary
        nutNonEquilibriumSpaldingWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, fvMesh>&,
            const dictionary&
        );

        //- Construct by mapping given onto a new patch
        nutNonEquilibriumSpaldingWallFunctionFvPatchScalarField
        (
            const nutNonEquilibriumSpaldingWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, fvMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        nutNonEquilibriumSpaldingWallFunctionFvPatchScalarField
        (
            const nutNonEquilibriumSpaldingWallFunctionFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        nutNonEquilibriumSpaldingWallFunctionFvPatchScalarField
        (
            const nutNonEquilibriumSpaldingWallFunctionFvPatchScalarField&,
            const DimensionedField<scalar, fvMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, fvMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new nutNonEquilibriumSpaldingWallFunctionFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        //- Friction velocity from Spalding's law with UTilde
        tmp<scalarField> uTau() const;

        //- Calculate and return the y+ at the boundary (uTau-based)
        virtual tmp<scalarField> yPlus() const;


    // Member Operators

        using nutNonEquilibriumWallFunctionFvPatchScalarField::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
