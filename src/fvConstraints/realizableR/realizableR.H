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
    Foam::fv::realizableR

Description
    Realizability limiter for a Reynolds-stress field.

    OpenFOAM's ReynoldsStress::boundNormalStress bounds only the normal
    stresses from below (kMin); the shear stresses are left free, and on
    layered meshes the stress tensor can lose positive definiteness (shear
    components grow while the normal ones sit at kMin) until the solution
    diverges.

    In every cell whose smallest eigenvalue of R falls below

        floor = max(minEigenvalue, relativeMinEigenvalue*tr(R)/3)

    the eigenvalues are clipped to that floor and R is rebuilt from the
    unchanged eigenvectors:

        R = sum_i max(lambda_i, floor) e_i e_i^T

    The eigenvectors of a clipped cell come from Jacobi rotations, which keep
    them orthonormal for repeated eigenvalues (Foam::eigenVectors returns the
    same vector twice there). Other cells and the boundary values are not
    changed. The Reynolds-stress
    models (LRR, SSG) call fvConstraints::constrain(R) right after solving
    for R, which is where this applies.

Usage
    Example usage, in system/fvConstraints:
    \verbatim
    realizableR
    {
        type        realizableR;
        libs        ("libkabutoFvConstraints.so");

        field       R;           // Optional, default R
        minEigenvalue          1e-6;  // Optional [m^2/s^2], default 1e-6
        relativeMinEigenvalue  1e-3;  // Optional, default 1e-3
    }
    \endverbatim

SourceFiles
    realizableR.C

\*---------------------------------------------------------------------------*/

#ifndef realizableR_H
#define realizableR_H

#include "fvConstraint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                         Class realizableR Declaration
\*---------------------------------------------------------------------------*/

class realizableR
:
    public fvConstraint
{
    // Private Data

        //- Name of the Reynolds-stress field
        word fieldName_;

        //- Absolute floor of the eigenvalues [m^2/s^2]
        scalar minEigenvalue_;

        //- Floor of the eigenvalues relative to tr(R)/3
        scalar relativeMinEigenvalue_;


    // Private Member Functions

        //- Read the coefficients
        void readCoeffs(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("realizableR");


    // Constructors

        //- Construct from components
        realizableR
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        realizableR(const realizableR&) = delete;


    //- Destructor
    virtual ~realizableR()
    {}


    // Static Member Functions

        //- Clip the eigenvalues of r to
        //  max(minEigenvalue, relativeMinEigenvalue*tr(r)/3), keeping its
        //  eigenvectors. Returns true if r was changed.
        static bool clip
        (
            symmTensor& r,
            const scalar minEigenvalue,
            const scalar relativeMinEigenvalue
        );


    // Member Functions

        //- Return the list of constrained fields
        virtual wordList constrainedFields() const;

        //- Clip the eigenvalues of R, returning true
        virtual bool constrain(volSymmTensorField& R) const;

        //- Update for mesh motion
        virtual bool movePoints();

        //- Update topology using the given map
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Redistribute or update using the given distribution map
        virtual void distribute(const polyDistributionMap&);

        //- Read dictionary
        virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const realizableR&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
