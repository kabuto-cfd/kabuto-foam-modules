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

\*---------------------------------------------------------------------------*/

#include "realizableR.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(realizableR, 0);
    addToRunTimeSelectionTable
    (
        fvConstraint,
        realizableR,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::realizableR::readCoeffs(const dictionary& dict)
{
    fieldName_ = dict.lookupOrDefault<word>("field", "R");
    minEigenvalue_ = dict.lookupOrDefault<scalar>("minEigenvalue", 1e-6);
    relativeMinEigenvalue_ =
        dict.lookupOrDefault<scalar>("relativeMinEigenvalue", 1e-3);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::realizableR::realizableR
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvConstraint(name, modelType, mesh, dict),
    fieldName_("R"),
    minEigenvalue_(1e-6),
    relativeMinEigenvalue_(1e-3)
{
    readCoeffs(coeffs(dict));
}


// * * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * //

namespace
{

// Eigen-decomposition of a symmetric tensor by cyclic Jacobi rotations.
// Returns the eigenvalues in lambda and the eigenvectors as the rows of V,
// orthonormal also for repeated eigenvalues (Foam::eigenVectors returns the
// same vector twice there).
void jacobiEigen(const Foam::symmTensor& t, Foam::vector& lambda, Foam::tensor& V)
{
    using namespace Foam;

    scalar a[3][3] =
    {
        {t.xx(), t.xy(), t.xz()},
        {t.xy(), t.yy(), t.yz()},
        {t.xz(), t.yz(), t.zz()}
    };
    scalar v[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

    const scalar scale = max(mag(t), vSmall);

    for (label sweep = 0; sweep < 50; sweep++)
    {
        const scalar off =
            sqr(a[0][1]) + sqr(a[0][2]) + sqr(a[1][2]);

        if (off < sqr(1e-15*scale))
        {
            break;
        }

        for (label p = 0; p < 2; p++)
        {
            for (label q = p + 1; q < 3; q++)
            {
                if (mag(a[p][q]) < vSmall)
                {
                    continue;
                }

                const scalar theta = (a[q][q] - a[p][p])/(2*a[p][q]);
                const scalar tn =
                    sign(theta)/(mag(theta) + Foam::sqrt(sqr(theta) + 1));
                const scalar c = 1/Foam::sqrt(sqr(tn) + 1);
                const scalar s = tn*c;

                for (label k = 0; k < 3; k++)
                {
                    const scalar akp = a[k][p], akq = a[k][q];
                    a[k][p] = c*akp - s*akq;
                    a[k][q] = s*akp + c*akq;
                }
                for (label k = 0; k < 3; k++)
                {
                    const scalar apk = a[p][k], aqk = a[q][k];
                    a[p][k] = c*apk - s*aqk;
                    a[q][k] = s*apk + c*aqk;
                }
                for (label k = 0; k < 3; k++)
                {
                    const scalar vkp = v[k][p], vkq = v[k][q];
                    v[k][p] = c*vkp - s*vkq;
                    v[k][q] = s*vkp + c*vkq;
                }
            }
        }
    }

    lambda = vector(a[0][0], a[1][1], a[2][2]);

    // Columns of v are the eigenvectors; return them as rows
    V = tensor
    (
        v[0][0], v[1][0], v[2][0],
        v[0][1], v[1][1], v[2][1],
        v[0][2], v[1][2], v[2][2]
    );
}

} // End anonymous namespace


bool Foam::fv::realizableR::clip
(
    symmTensor& r,
    const scalar minEigenvalue,
    const scalar relativeMinEigenvalue
)
{
    const scalar floor =
        max(minEigenvalue, relativeMinEigenvalue*max(tr(r), 0)/3);

    // Cheap test first: OpenFOAM's analytical eigenvalues
    if (cmptMin(eigenValues(r)) >= floor)
    {
        return false;
    }

    vector lambda;
    tensor ev;
    jacobiEigen(r, lambda, ev);

    r =
        max(lambda.x(), floor)*sqr(ev.x())
      + max(lambda.y(), floor)*sqr(ev.y())
      + max(lambda.z(), floor)*sqr(ev.z());

    return true;
}


// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

Foam::wordList Foam::fv::realizableR::constrainedFields() const
{
    return wordList(1, fieldName_);
}


bool Foam::fv::realizableR::constrain(volSymmTensorField& R) const
{
    symmTensorField& Ri = R.primitiveFieldRef();

    label nClipped = 0;

    forAll(Ri, celli)
    {
        if (clip(Ri[celli], minEigenvalue_, relativeMinEigenvalue_))
        {
            nClipped++;
        }
    }

    R.correctBoundaryConditions();

    Info<< type() << ' ' << name() << ": "
        << returnReduce(nClipped, sumOp<label>()) << " of "
        << returnReduce(Ri.size(), sumOp<label>()) << " cells clipped"
        << endl;

    return true;
}


bool Foam::fv::realizableR::movePoints()
{
    return true;
}


void Foam::fv::realizableR::topoChange(const polyTopoChangeMap&)
{}


void Foam::fv::realizableR::mapMesh(const polyMeshMap&)
{}


void Foam::fv::realizableR::distribute(const polyDistributionMap&)
{}


bool Foam::fv::realizableR::read(const dictionary& dict)
{
    if (fvConstraint::read(dict))
    {
        readCoeffs(coeffs(dict));
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
