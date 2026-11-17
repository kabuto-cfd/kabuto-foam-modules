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
#include <cmath>

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
    // Work on r scaled by its largest component: the cubic of the analytical
    // eigenvalues and the squares in the Jacobi rotations overflow for very
    // large R (a diverging run) otherwise
    const scalar s = max(cmptMax(cmptMag(r)), vSmall);
    const symmTensor rs(r/s);

    const scalar floor =
        max(minEigenvalue, relativeMinEigenvalue*max(tr(r), 0)/3)/s;

    // Cheap test first: OpenFOAM's analytical eigenvalues
    if (cmptMin(eigenValues(rs)) >= floor)
    {
        return false;
    }

    vector lambda;
    tensor ev;
    jacobiEigen(rs, lambda, ev);

    r =
        s
       *(
            max(lambda.x(), floor)*sqr(ev.x())
          + max(lambda.y(), floor)*sqr(ev.y())
          + max(lambda.z(), floor)*sqr(ev.z())
        );

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
    const scalarField& V = mesh().V();
    const vectorField& C = mesh().C();

    // Stop with the location rather than a floating-point exception in the
    // eigenvalue computation when R is not finite
    label nNonFinite = 0;
    label firstNonFinite = -1;

    forAll(Ri, celli)
    {
        const symmTensor& r = Ri[celli];

        for (direction cmpt = 0; cmpt < symmTensor::nComponents; cmpt++)
        {
            if (!std::isfinite(r[cmpt]))
            {
                if (firstNonFinite < 0)
                {
                    firstNonFinite = celli;
                }
                nNonFinite++;
                break;
            }
        }
    }

    if (returnReduce(nNonFinite, sumOp<label>()))
    {
        if (firstNonFinite >= 0)
        {
            Pout<< type() << ' ' << name() << ": " << nNonFinite
                << " cells with non-finite " << fieldName_
                << ", first: cell " << firstNonFinite
                << " at " << C[firstNonFinite]
                << ", " << fieldName_ << " = " << Ri[firstNonFinite] << endl;
        }

        FatalErrorInFunction
            << "Non-finite " << fieldName_ << " in "
            << returnReduce(nNonFinite, sumOp<label>()) << " cells"
            << " (locations above)" << exit(FatalError);
    }

    label nClipped = 0;

    // Turbulent kinetic energy added by the clipping: tr(R)/2 after - before
    scalar dkV = 0;
    scalar kV = 0;
    scalar dkMax = 0;
    label celldkMax = -1;

    forAll(Ri, celli)
    {
        const scalar tr0 = tr(Ri[celli]);
        kV += 0.5*tr0*V[celli];

        if (clip(Ri[celli], minEigenvalue_, relativeMinEigenvalue_))
        {
            nClipped++;

            const scalar dk = 0.5*(tr(Ri[celli]) - tr0);
            dkV += dk*V[celli];

            if (dk > dkMax)
            {
                dkMax = dk;
                celldkMax = celli;
            }
        }
    }

    R.correctBoundaryConditions();

    // Location of the largest added k over all processors
    const scalar dkMaxAll = returnReduce(dkMax, maxOp<scalar>());
    vector CdkMax(-great*vector::one);
    scalar kdkMax = -great;
    if (celldkMax >= 0 && dkMax == dkMaxAll)
    {
        CdkMax = C[celldkMax];
        kdkMax = 0.5*tr(Ri[celldkMax]) - dkMax;
    }
    reduce(CdkMax, maxOp<vector>());
    reduce(kdkMax, maxOp<scalar>());

    Info<< type() << ' ' << name() << ": "
        << returnReduce(nClipped, sumOp<label>()) << " of "
        << returnReduce(Ri.size(), sumOp<label>()) << " cells clipped, "
        << "k added (volume integral) "
        << returnReduce(dkV, sumOp<scalar>())/max(returnReduce(kV, sumOp<scalar>()), vSmall)
        << " of the total";

    if (dkMaxAll > 0)
    {
        Info<< ", max " << dkMaxAll << " (k before " << kdkMax
            << ") at " << CdkMax;
    }

    Info<< endl;

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
