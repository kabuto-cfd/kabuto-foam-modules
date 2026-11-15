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

bool Foam::fv::realizableR::clip
(
    symmTensor& r,
    const scalar minEigenvalue,
    const scalar relativeMinEigenvalue
)
{
    const vector lambda = eigenValues(r);
    const scalar floor =
        max(minEigenvalue, relativeMinEigenvalue*max(tr(r), 0)/3);

    if (cmptMin(lambda) >= floor)
    {
        return false;
    }

    // Rows of ev are the eigenvectors, in the order of lambda
    const tensor ev = eigenVectors(r, lambda);

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
