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
      With Tf auto (the default) the time scale is set per face to that
      crossing time times TfCrossings, deltaX the square root of the wall
      face area and U the larger of |U_s| and |U_f|:
      \verbatim
          eps = min(deltaT max(|U_s|, |U_f|)/(TfCrossings |S_f|^1/2), 1)
      \endverbatim
      so no time scale needs choosing; Tf 0 switches the filter off.

    The defaults (sampleLayer 1, Tf auto, no floor) came out best in a
    channel at Re_tau 2000 and the periodic hill at Re_H 10595 (DNS of
    Krank et al. 2018): in the channel the filter alone brings the wall
    stress to within 4 % of Dean's correlation; a sample further out helps
    the attached channel a little but shortens the separation bubble of the
    hill, where it reads the outer flow instead of the backflow near the
    wall.

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

    The stress set at the wall has to reach the flow outside through the
    cells of the wall-modelled layer, by the resolved and the subgrid stress.
    Subgrid models built for resolved walls (WALE, sigma) fade as y^3 there
    and may carry too little of it, so the velocity jumps between the first
    cells: with the wall-cell input the second cell sits above the log law,
    with a sample further out the wall cell falls below it. With nutFloor
    the eddy viscosity of the cells in layers 1 .. max(sampleLayer, 2) is
    raised (by the fv::wallModelledNut constraint, which must be selected in
    system/fvConstraints) to at least the eddy viscosity of the law itself,
    with this face's friction velocity:
    \verbatim
        nut_floor = nu (dy+/dU+ - 1),  y+ = y uTau/nu   (Spalding's law)
                 ~ kappa y uTau in the log layer
    \endverbatim
    the viscosity that carries the wall stress across a constant-stress layer
    following the law: the wall model's own assumption, with no further
    constant. It is the idea of the two-part and near-wall eddy-viscosity
    models of Schumann (1975), Sullivan et al. (1994) and Templeton et al.
    (2005) restricted to the wall-modelled layer. It is off by default: the
    resolved stress in the layer carries part of the wall stress too, so the
    full law viscosity overcounts it; in the channel the velocity then falls
    below the log law and the wall stress by 3 - 16 %, and in the periodic
    hill it did not improve the separation bubble on the finer grid.

    References:
    \verbatim
        Spalding, D. B. (1961). A single formula for the law of the wall.
        Journal of Applied Mechanics, 28(3), 455-458.

        Schumann, U. (1975). Subgrid scale model for finite difference
        simulations of turbulent flows in plane channels and annuli. Journal
        of Computational Physics, 18(4), 376-404.

        Sullivan, P. P., McWilliams, J. C., Moeng, C.-H. (1994). A subgrid-
        scale model for large-eddy simulation of planetary boundary-layer
        flows. Boundary-Layer Meteorology, 71(3), 247-276.

        Templeton, J. A., Medic, G., Kalitzin, G. (2005). An eddy-viscosity
        based near-wall treatment for coarse grid large-eddy simulation.
        Physics of Fluids, 17(10), 105101.

        Krank, B., Kronbichler, M., Wall, W. A. (2018). Direct numerical
        simulation of flow over periodic hills up to Re_H = 10,595. Flow,
        Turbulence and Combustion, 101(2), 521-551.

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
        Tf          | time scale of the input filter [s], 0: none, auto | no | auto
        TfCrossings | wall-cell crossing times per Tf (Tf auto)  | no       | 4
        nutFloor    | floor nut in the wall-modelled layer       | no       | no
        kappa       | von Karman constant                        | no       | 0.41
        E           | log-law constant                           | no       | 9.8
        Uf          | filtered input velocity (written, restart) | no       |
    \endtable

    In 0/nut, with libs ("libkabutoWallFunctions.so"); in system/controlDict:
    \verbatim
    <patchName>
    {
        type            nutSpaldingLESWallFunction;
        value           uniform 0;
    }
    \endverbatim
    (sampleLayer 1, Tf auto). The floor, if wanted, takes nutFloor yes here
    and the fv::wallModelledNut constraint in system/fvConstraints.
    and with nutFloor, in system/fvConstraints:
    \verbatim
    wallModelledNut
    {
        type            wallModelledNut;
        libs            ("libkabutoWallFunctions.so");
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

        //- Filter time scale set per face from the wall cell's crossing time
        bool TfAuto_;

        //- Crossing times per filter time scale (Tf auto)
        scalar TfCrossings_;

        //- Floor the eddy viscosity of the wall-modelled layer (cells up to
        //  layer max(sampleLayer, 2)) at the eddy viscosity of the law
        bool nutFloor_;

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

        //- Cell of each face in layers 1 .. max(sampleLayer, 2) (-1 where
        //  the walk stopped early)
        mutable List<labelList> layerCells_;

        //- Their distances from the wall plane
        mutable List<scalarField> layerY_;

        //- Friction velocity of the last evaluation
        mutable scalarField uTau_;


    // Private Member Functions

        //- Find the sample cells and distances
        void calcSampling() const;

        //- Find them if not done yet (on every processor)
        void checkSampling() const;

        //- Whether the input velocity is filtered in time
        bool filtered() const
        {
            return Tf_ > 0 || TfAuto_;
        }

        //- Number of cell layers walked: the sample layer, and at least two
        //  (the floor needs both cells of the face between layers 1 and 2)
        label nWalk() const
        {
            return max(sampleLayer_, label(2));
        }

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

        //- Whether the eddy viscosity floor is on
        bool nutFloor() const
        {
            return nutFloor_;
        }

        //- Raise the cell values of nut in the wall-modelled layer to at
        //  least the eddy viscosity of Spalding's law with this patch's
        //  friction velocity; returns the number of cells raised. Called by
        //  fv::wallModelledNut after the subgrid model has set nut
        label limitNut(scalarField& nut) const;


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
