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

#include "twoLayerEpsilon.H"
#include "momentumTransportModel.H"
#include "wallDist.H"
#include "wallFvPatch.H"
#include "fvMatrix.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(twoLayerEpsilon, 0);
    addToRunTimeSelectionTable
    (
        fvConstraint,
        twoLayerEpsilon,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::twoLayerEpsilon::readCoeffs(const dictionary& dict)
{
    epsilonName_ = dict.lookupOrDefault<word>("epsilon", "epsilon");
    nutName_ = dict.lookupOrDefault<word>("nut", "nut");
    kName_ = dict.lookupOrDefault<word>("k", "k");
    Cmu_ = dict.lookupOrDefault<scalar>("Cmu", 0.09);
    kappa_ = dict.lookupOrDefault<scalar>("kappa", 0.41);
    Amu_ = dict.lookupOrDefault<scalar>("Amu", 70);
    ReyStar_ = dict.lookupOrDefault<scalar>("ReyStar", 200);
    deltaReyStar_ = dict.lookupOrDefault<scalar>("deltaReyStar", 20);
}


void Foam::fv::twoLayerEpsilon::blending
(
    scalarField& Rey,
    scalarField& lambda
) const
{
    const volScalarField& k = mesh().lookupObject<volScalarField>(kName_);
    const volScalarField& y = wallDist::New(mesh()).y();
    const momentumTransportModel& mtm =
        mesh().lookupType<momentumTransportModel>();
    const tmp<volScalarField> tnu = mtm.nu();
    const volScalarField& nu = tnu();

    const scalar A = deltaReyStar_/atanh(0.98);

    Rey = y.primitiveField()*sqrt(max(k.primitiveField(), scalar(0)))
       /nu.primitiveField();
    lambda = 0.5*(1 + tanh((Rey - ReyStar_)/A));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::twoLayerEpsilon::twoLayerEpsilon
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

Foam::wordList Foam::fv::twoLayerEpsilon::constrainedFields() const
{
    return wordList({epsilonName_, nutName_});
}


bool Foam::fv::twoLayerEpsilon::constrain(volScalarField& field) const
{
    const bool isEpsilon = field.name() == epsilonName_;
    const bool isNut = field.name() == nutName_;

    if (!isEpsilon && !isNut)
    {
        return false;
    }

    const volScalarField& k = mesh().lookupObject<volScalarField>(kName_);
    const volScalarField& y = wallDist::New(mesh()).y();

    scalarField Rey(mesh().nCells()), lambda(mesh().nCells());
    blending(Rey, lambda);

    // The wall cells keep the epsilon the wall functions give them
    boolList wallCell(mesh().nCells(), false);
    forAll(mesh().boundary(), patchi)
    {
        if (isA<wallFvPatch>(mesh().boundary()[patchi]))
        {
            UIndirectList<bool>
            (
                wallCell,
                mesh().boundary()[patchi].faceCells()
            ) = true;
        }
    }

    const scalar Cl = kappa_*pow(Cmu_, -0.75);
    const scalar A = isEpsilon ? 2*Cl : Amu_;

    scalarField& psi = field.primitiveFieldRef();

    label n = 0;
    forAll(psi, celli)
    {
        if ((isEpsilon && wallCell[celli]) || lambda[celli] > 1 - 1e-3)
        {
            continue;
        }

        const scalar kc = max(k[celli], scalar(0));
        const scalar l =
            max(y[celli]*Cl*(1 - exp(-Rey[celli]/A)), vSmall);

        const scalar psi2l = isEpsilon ? pow(kc, 1.5)/l : Cmu_*l*sqrt(kc);

        psi[celli] = lambda[celli]*psi[celli] + (1 - lambda[celli])*psi2l;
        n++;
    }

    return n;
}


bool Foam::fv::twoLayerEpsilon::read(const dictionary& dict)
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
