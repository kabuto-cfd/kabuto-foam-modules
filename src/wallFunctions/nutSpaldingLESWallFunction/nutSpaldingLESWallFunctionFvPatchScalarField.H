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
    Foam::nutSpaldingLESWallFunctionFvPatchScalarField

Description
    Wall-stress model for wall-modelled LES: Spalding's law (wallLaw.H) fed
    with a velocity that is not the raw velocity of the wall cell, the two
    standard remedies for the log-layer mismatch.

    With the raw wall-cell velocity (nutUSpaldingWallFunction), the modelled
    wall stress follows the instantaneous velocity of the same cell, damps the
    resolved turbulence next to the wall, and the mean velocity of the cells
    above it sits above the log law (Yang et al. 2017); the wall stress is
    then too low for the flow outside, by about 15 % in a channel and more
    with a subgrid model that switches off in shear near the wall (WALE,
    sigma). The input velocity U_in can be decoupled from the wall cell

    - in space: taken in the cell sampleLayer layers away from the wall
      (Kawai and Larsson 2012), found by walking from the wall cell through
      the face whose normal is closest to the wall normal (prism layers:
      the cell straight above); the walk stops early at a processor or
      domain boundary;
    - in time: an exponential filter with time scale Tf (Yang et al. 2017),
      \verbatim
          U_f^n = (1 - eps) U_f^(n-1) + eps U_s^n,   eps = min(deltaT/Tf, 1)
      \endverbatim
      updated once per time step and written with the field for restarts.
      Yang et al. find the mismatch removed for Tf > deltaX/U (the time the
      flow takes to cross a wall cell) and unchanged up to much longer Tf.

    Both may be combined. U_s is the velocity of the sample cell relative to
    the wall, less its component normal to the wall face; y_s the distance of
    that cell centre from the wall plane. The friction velocity solves
    Spalding's law at the sample point and the wall stress is applied
    through the wall cell, in the direction of its velocity:
    \verbatim
        y_s uTau/nu = y+(|U_in|/uTau)             (Spalding's law)
        nut_w       = max(uTau^2/|dU/dn|_w - nu, 0)
    \endverbatim
    With sampleLayer 1 and Tf 0 this is nutUSpaldingWallFunction (solved to
    round-off).

    References:
    \verbatim
        Spalding, D. B. (1961). A single formula for the law of the wall.
        Journal of Applied Mechanics, 28(3), 455-458.

        Kawai, S., Larsson, J. (2012). Wall-modeling in large eddy
        simulation: length scales, grid resolution, and accuracy.
        Physics of Fluids, 24(1), 015105.

        Yang, X. I. A., Park, G. I., Moin, P. (2017). Log-layer mismatch and
        modeling of the fluctuating wall stress in wall-modeled large-eddy
        simulations. Physical Review Fluids, 2(10), 104601.
    \endverbatim

Usage
    \table
        Property    | Description                                | Required | Default
        sampleLayer | cell layer of the input velocity (1 = wall cell) | no | 1
        Tf          | time scale of the input filter [s], 0: none | no      | 0
        kappa       | von Karman constant                        | no       | 0.41
        E           | log-law constant                           | no       | 9.8
        Uf          | filtered input velocity (written, restart) | no       |
    \endtable

    In 0/nut, with libs ("libkabutoWallFunctions.so"); in system/controlDict:
    \verbatim
    <patchName>
    {
        type            nutSpaldingLESWallFunction;
        sampleLayer     2;
        Tf              1e-3;
        value           uniform 0;
    }
    \endverbatim

SourceFiles
    nutSpaldingLESWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nutSpaldingLESWallFunctionFvPatchScalarField_H
#define nutSpaldingLESWallFunctionFvPatchScalarField_H

#include "nutWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class nutSpaldingLESWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutSpaldingLESWallFunctionFvPatchScalarField
:
    public nutWallFunctionFvPatchScalarField
{
    // Private Data

        //- Cell layer of the input velocity (1 = the wall cell)
        label sampleLayer_;

        //- Time scale of the exponential input filter (0: no filter)
        scalar Tf_;

        //- Filtered input velocity (wall-parallel, relative to the wall)
        vectorField Uf_;

        //- Time index of the last filter update
        label filterTimeIndex_;

        //- Whether the sample cells have been found. Not the size of
        //  sampleCells_: a processor may hold no face of the patch, and every
        //  processor must take part in the reductions of calcSampling
        mutable bool sampled_;

        //- Sample cell of each face (demand-driven)
        mutable labelList sampleCells_;

        //- Distance of the sample cell centre from the wall plane
        mutable scalarField ySample_;


    // Private Member Functions

        //- Find the sample cells and distances
        void calcSampling() const;

        //- Find them if not done yet (on every processor)
        void checkSampling() const;

        //- Wall-parallel velocity of the sample cells relative to the wall
        tmp<vectorField> sampleVelocity() const;

        //- Input velocity of the wall model: filtered, else sampled
        tmp<vectorField> inputVelocity() const;

        //- Friction velocity from Spalding's law at the sample point
        tmp<scalarField> uTau() const;


protected:

    // Protected Member Functions

        //- Calculate the turbulent viscosity
        virtual tmp<scalarField> nut() const;

        //- Write local wall function variables
        virtual void writeLocalEntries(Ostream&) const;


public:

    //- Runtime type information
    TypeName("nutSpaldingLESWallFunction");


    // Constructors

        //- Construct from patch, internal field and dictionary
        nutSpaldingLESWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, fvMesh>&,
            const dictionary&
        );

        //- Construct by mapping given onto a new patch
        nutSpaldingLESWallFunctionFvPatchScalarField
        (
            const nutSpaldingLESWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, fvMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        nutSpaldingLESWallFunctionFvPatchScalarField
        (
            const nutSpaldingLESWallFunctionFvPatchScalarField&
        ) = delete;

        //- Copy constructor setting internal field reference
        nutSpaldingLESWallFunctionFvPatchScalarField
        (
            const nutSpaldingLESWallFunctionFvPatchScalarField&,
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
                new nutSpaldingLESWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        //- Calculate and return the y+ of the wall cell (model uTau)
        virtual tmp<scalarField> yPlus() const;


        // Mapping functions

            //- Map the given fvPatchField onto this fvPatchField
            virtual void map(const fvPatchScalarField&, const fieldMapper&);

            //- Reset the fvPatchField to the given fvPatchField
            virtual void reset(const fvPatchScalarField&);


        // Evaluation functions

            //- Update the input filter once per time step, then nut
            virtual void updateCoeffs();


    // Member Operators

        using nutWallFunctionFvPatchScalarField::operator=;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
