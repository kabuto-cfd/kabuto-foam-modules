/*---------------------------------------------------------------------------*\
  Kabuto CFD modules for OpenFOAM-14
  Copyright (C) 2026 Kabuto CFD
  Derived from epsilonWallFunctionFvPatchScalarField of OpenFOAM
  Copyright (C) 2011-2026 OpenFOAM Foundation
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
    Foam::omegaNonEquilibriumWallFunctionFvPatchScalarField

Description
    Kim and Choudhury's (1995) non-equilibrium wall function for omega and
    G, to go with nutNonEquilibriumWallFunction: the production and the
    dissipation in the wall cells are the averages over the cell height
    y_n = 2 y of a two-layer profile, viscous below y_v and logarithmic
    above it (Fluent theory guide), and omega = epsilon/(Cmu k):

    \verbatim
        y_v     = yvStar nu/(Cmu^1/4 k^1/2)
        G       = tau_w^2/(kappa Cmu^1/4 k^1/2 y_n) ln(y_n/y_v)
        epsilon = k/y_n [ 2 nu/y_v + Cmu^3/4 k^1/2/kappa ln(y_n/y_v) ]
        omega   = epsilon/(Cmu k)
    \endverbatim

    with tau_w = (nu + nut_w)|dU/dn| from the nut wall function, which carries
    the pressure-gradient sensitivity. Below y* = yPlusLam (viscous wall cell)
    omega = 6 nu/(beta1 y^2) and G is left as the model computed it.

    Use the same wall function family on all walls of a case (the wall cells
    of all wall-function patches are handled together by one of them).

Usage
    \table
        Property | Description              | Required | Default
        beta1    | model coefficient        | no       | 0.075
    \endtable
    \verbatim
    <patchName>
    {
        type            omegaNonEquilibriumWallFunction;
        value           $internalField;
    }
    \endverbatim
    Cmu, kappa, E and yvStar are taken from the nut wall function.

SourceFiles
    omegaNonEquilibriumWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef omegaNonEquilibriumWallFunctionFvPatchScalarField_H
#define omegaNonEquilibriumWallFunctionFvPatchScalarField_H

#include "wallCellWallFunctionFvPatchScalarField.H"
#include "Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class momentumTransportModel;

/*---------------------------------------------------------------------------*\
            Class omegaNonEquilibriumWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class omegaNonEquilibriumWallFunctionFvPatchScalarField
:
    public wallCellWallFunctionFvPatchScalarField
{
private:

    // Private Data

        //- Wall cell turbulence generation values. Master patch only.
        autoPtr<scalarField> wallCellGPtr_;

        //- Wall cell turbulence dissipation values. Master patch only.
        autoPtr<scalarField> wallCellOmegaPtr_;


    // Private Member Functions


        //- Calculate the turbulent values in wall cells and apply to the
        //  turbulence model's fields
        void updateCoeffsMaster();

        //- Constrain the matrix such that the wall cells take the value
        //  specified by the wall function
        void manipulateMatrixMaster(fvMatrix<scalar>& matrix);


protected:

    // Protected Data

        //- beta1 coefficient (viscous omega)
        scalar beta1_;


    // Protected Member Functions

        //- Calculate and return G and the wall-cell value
        virtual Pair<tmp<scalarField>> calculate
        (
            const momentumTransportModel&
        ) const;


public:

    //- Runtime type information
    TypeName("omegaNonEquilibriumWallFunction");


    // Constructors

        //- Construct from patch, internal field and dictionary
        omegaNonEquilibriumWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, fvMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  omegaNonEquilibriumWallFunctionFvPatchScalarField
        //  onto a new patch
        omegaNonEquilibriumWallFunctionFvPatchScalarField
        (
            const omegaNonEquilibriumWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, fvMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        omegaNonEquilibriumWallFunctionFvPatchScalarField
        (
            const omegaNonEquilibriumWallFunctionFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        omegaNonEquilibriumWallFunctionFvPatchScalarField
        (
            const omegaNonEquilibriumWallFunctionFvPatchScalarField&,
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
                new omegaNonEquilibriumWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    //- Destructor
    virtual ~omegaNonEquilibriumWallFunctionFvPatchScalarField()
    {}


    // Member Functions

        // Mapping functions

            //- Map the given fvPatchField onto this fvPatchField
            virtual void map(const fvPatchScalarField&, const fieldMapper&);

            //- Reset the fvPatchField to the given fvPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const fvPatchScalarField&);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Manipulate matrix
            virtual void manipulateMatrix(fvMatrix<scalar>& matrix);

        // I-O

            //- Write
            virtual void write(Ostream&) const;


    // Member Operators

        //- Inherit assignment
        using wallCellWallFunctionFvPatchScalarField::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
