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

#include "wallModelledNut.H"
#include "nutSpaldingLESWallFunctionFvPatchScalarField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(wallModelledNut, 0);
    addToRunTimeSelectionTable
    (
        fvConstraint,
        wallModelledNut,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::wallModelledNut::readCoeffs(const dictionary& dict)
{
    nutName_ = dict.lookupOrDefault<word>("nut", "nut");
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::wallModelledNut::wallModelledNut
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvConstraint(name, modelType, mesh, dict)
{
    readCoeffs(coeffs(dict));
}


// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

Foam::wordList Foam::fv::wallModelledNut::constrainedFields() const
{
    return wordList(1, nutName_);
}


bool Foam::fv::wallModelledNut::constrain(volScalarField& field) const
{
    if (field.name() != nutName_)
    {
        return false;
    }

    scalarField& nut = field.primitiveFieldRef();

    label n = 0, nPatches = 0;
    forAll(field.boundaryField(), patchi)
    {
        if
        (
            isA<nutSpaldingLESWallFunctionFvPatchScalarField>
            (
                field.boundaryField()[patchi]
            )
        )
        {
            const nutSpaldingLESWallFunctionFvPatchScalarField& wf =
                refCast<const nutSpaldingLESWallFunctionFvPatchScalarField>
                (
                    field.boundaryField()[patchi]
                );

            if (wf.nutFloor())
            {
                n += wf.limitNut(nut);
                nPatches++;
            }
        }
    }

    // Processor (and coupled) patches must see the raised cells; every
    // processor takes part
    field.correctBoundaryConditions();

    if (debug)
    {
        Info<< type() << ": raised nut in " << returnReduce(n, sumOp<label>())
            << " cells" << endl;
    }

    return returnReduce(nPatches, sumOp<label>()) > 0;
}


bool Foam::fv::wallModelledNut::read(const dictionary& dict)
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
