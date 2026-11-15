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
    Foam::omegaNonEquilibriumSpaldingWallFunctionFvPatchScalarField

Description
    y+-insensitive non-equilibrium wall function for omega and G, to go with
    nutNonEquilibriumSpaldingWallFunction. With uTau = ((nu + nut_w)|dU/dn|)^1/2
    from the nut wall function, y+ = uTau y/nu and Spalding's profile U+(y+):

    \verbatim
        gamma = dU+/dy+                         (viscous share of the stress)
        tau+  = 1 + p+ y+,  p+ = nu (dp/ds)/uTau^3     (total stress / tau_w)
        G     = uTau^4/nu (tau+ - gamma) gamma
        omega = (omegaVis^2 + omegaLog^2)^1/2
        omegaVis = 6 nu/(beta1 y^2),  omegaLog = uTau/(Cmu^1/2 kappa y)
    \endverbatim

    G is the production of the Spalding profile under the local total shear
    stress: zero at the wall, 0.25 uTau^4/nu at its peak in the buffer layer,
    tau_w uTau/(kappa y) in the log layer. omega blends the viscous and the
    logarithmic values as Menter's automatic wall treatment.

    Use the same wall function family on all walls of a case.

Usage
    \verbatim
    <patchName>
    {
        type            omegaNonEquilibriumSpaldingWallFunction;
        value           $internalField;
    }
    \endverbatim
    Cmu, kappa and E are taken from the nut wall function; beta1 (0.075) as
    omegaNonEquilibriumWallFunction.

SourceFiles
    omegaNonEquilibriumSpaldingWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef omegaNonEquilibriumSpaldingWallFunctionFvPatchScalarField_H
#define omegaNonEquilibriumSpaldingWallFunctionFvPatchScalarField_H

#include "omegaNonEquilibriumWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class omegaNonEquilibriumSpaldingWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class omegaNonEquilibriumSpaldingWallFunctionFvPatchScalarField
:
    public omegaNonEquilibriumWallFunctionFvPatchScalarField
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
    TypeName("omegaNonEquilibriumSpaldingWallFunction");


    // Constructors

        //- Construct from patch, internal field and dictionary
        omegaNonEquilibriumSpaldingWallFunctionFvPatchScalarField
        (
            const fvPatch& p,
            const DimensionedField<scalar, fvMesh>& iF,
            const dictionary& dict
        )
        :
            omegaNonEquilibriumWallFunctionFvPatchScalarField(p, iF, dict)
        {}

        //- Construct by mapping given onto a new patch
        omegaNonEquilibriumSpaldingWallFunctionFvPatchScalarField
        (
            const omegaNonEquilibriumSpaldingWallFunctionFvPatchScalarField& ptf,
            const fvPatch& p,
            const DimensionedField<scalar, fvMesh>& iF,
            const fieldMapper& mapper
        )
        :
            omegaNonEquilibriumWallFunctionFvPatchScalarField(ptf, p, iF, mapper)
        {}

        //- Disallow copy without setting internal field reference
        omegaNonEquilibriumSpaldingWallFunctionFvPatchScalarField
        (
            const omegaNonEquilibriumSpaldingWallFunctionFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        omegaNonEquilibriumSpaldingWallFunctionFvPatchScalarField
        (
            const omegaNonEquilibriumSpaldingWallFunctionFvPatchScalarField& wfpsf,
            const DimensionedField<scalar, fvMesh>& iF
        )
        :
            omegaNonEquilibriumWallFunctionFvPatchScalarField(wfpsf, iF)
        {}

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, fvMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new omegaNonEquilibriumSpaldingWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member Operators

        using omegaNonEquilibriumWallFunctionFvPatchScalarField::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
