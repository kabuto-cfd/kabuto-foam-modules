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

    The structure follows nutkWallFunctionFvPatchScalarField of OpenFOAM
    (Copyright (C) 2011-2026 OpenFOAM Foundation).

Class
    Foam::nutNonEquilibriumWallFunctionFvPatchScalarField

Description
    Kim and Choudhury's (1995) non-equilibrium wall function for the turbulent
    viscosity: the log law of Launder and Spalding, with the velocity scale
    from k, sensitised to the pressure gradient along the wall.

    \verbatim
        y*     = Cmu^1/4 k^1/2 y/nu
        tau_w  = Cmu^1/4 k^1/2 kappa UTilde/ln(E y*)      y* > yPlusLam
               = nu U/y                                   y* <= yPlusLam
        nut_w  = tau_w y/U - nu
    \endverbatim

    with UTilde the pressure-gradient-sensitised velocity of wallLaw.H (the
    tangential velocity at the cell centre corrected by the pressure gradient
    along it, limited to 0 <= UTilde <= 2 U). Without a pressure gradient
    this is nutkWallFunction. The wall shear never reverses:
    nut_w >= -0.999 nu.

    Use it with epsilonNonEquilibriumWallFunction or
    omegaNonEquilibriumWallFunction, which take the production and
    dissipation in the wall cells from the same two-layer picture.
    Wilcox, Fluent and others recommend the first cell in the log layer
    (y* > 30); for a treatment that does not depend on y+ use
    nutNonEquilibriumSpaldingWallFunction.

    Reference:
    \verbatim
        Kim, S.-E., Choudhury, D. (1995).
        A near-wall treatment using wall functions sensitized to pressure
        gradient. ASME FED Vol. 217, Separated and Complex Flows.
    \endverbatim

Usage
    \table
        Property | Description                       | Required | Default
        Cmu      | model coefficient                 | no       | 0.09
        kappa    | von Karman constant               | no       | 0.41
        E        | log-law constant                  | no       | 9.8
        yvStar   | viscous sublayer thickness in y*  | no       | 11.225
        p        | pressure field                    | no       | p
        rho      | density field (p in Pa only)      | no       | rho
    \endtable

    A pressure field with the dimensions of pressure is divided by the
    density; a kinematic one (incompressible solvers) is used as it is.
    Without a pressure field the wall function is nutkWallFunction.

    \verbatim
    <patchName>
    {
        type            nutNonEquilibriumWallFunction;
        value           uniform 0;
    }
    \endverbatim

SourceFiles
    nutNonEquilibriumWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nutNonEquilibriumWallFunctionFvPatchScalarField_H
#define nutNonEquilibriumWallFunctionFvPatchScalarField_H

#include "nutWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class nutNonEquilibriumWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutNonEquilibriumWallFunctionFvPatchScalarField
:
    public nutWallFunctionFvPatchScalarField
{
protected:

    // Protected Data

        //- Viscous sublayer thickness in y*
        scalar yvStar_;

        //- Name of the pressure field
        word pName_;

        //- Name of the density field
        word rhoName_;


    // Protected Member Functions

        //- Calculate the turbulent viscosity
        virtual tmp<scalarField> nut() const;

        //- Write local wall function variables
        virtual void writeLocalEntries(Ostream&) const;


public:

    //- Runtime type information
    TypeName("nutNonEquilibriumWallFunction");


    // Constructors

        //- Construct from patch, internal field and dictionary
        nutNonEquilibriumWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, fvMesh>&,
            const dictionary&
        );

        //- Construct by mapping given onto a new patch
        nutNonEquilibriumWallFunctionFvPatchScalarField
        (
            const nutNonEquilibriumWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, fvMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        nutNonEquilibriumWallFunctionFvPatchScalarField
        (
            const nutNonEquilibriumWallFunctionFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        nutNonEquilibriumWallFunctionFvPatchScalarField
        (
            const nutNonEquilibriumWallFunctionFvPatchScalarField&,
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
                new nutNonEquilibriumWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        //- Viscous sublayer thickness in y*
        scalar yvStar() const
        {
            return yvStar_;
        }

        //- Tangential velocity relative to the wall at the wall-cell centres
        tmp<scalarField> magUt() const;

        //- Kinematic pressure gradient along the tangential velocity at the
        //  wall-cell centres (zero without a pressure field)
        tmp<scalarField> dpds() const;

        //- Pressure-gradient-sensitised velocity at the wall-cell centres
        tmp<scalarField> UTilde() const;

        //- Calculate and return the y+ at the boundary (k-based, y*)
        virtual tmp<scalarField> yPlus() const;


    // Member Operators

        using nutWallFunctionFvPatchScalarField::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
