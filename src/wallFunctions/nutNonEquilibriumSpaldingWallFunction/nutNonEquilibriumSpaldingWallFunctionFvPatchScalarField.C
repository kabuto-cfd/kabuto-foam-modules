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

#include "nutNonEquilibriumSpaldingWallFunctionFvPatchScalarField.H"
#include "wallLaw.H"
#include "momentumTransportModel.H"
#include "fieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::nutNonEquilibriumSpaldingWallFunctionFvPatchScalarField::nut() const
{
    const label patchi = patch().index();
    const momentumTransportModel& turbModel =
        db().lookupType<momentumTransportModel>(internalField().group());

    const scalarField& y = turbModel.yb()[patchi];
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const scalarField Ut(magUt());
    const scalarField ut(uTau());

    tmp<scalarField> tnutw(new scalarField(patch().size(), 0.0));
    scalarField& nutw = tnutw.ref();

    forAll(nutw, facei)
    {
        if (Ut[facei] > vSmall)
        {
            nutw[facei] = max
            (
                sqr(ut[facei])*y[facei]/Ut[facei] - nuw[facei],
               -0.999*nuw[facei]
            );
        }
    }

    return tnutw;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nutNonEquilibriumSpaldingWallFunctionFvPatchScalarField::
nutNonEquilibriumSpaldingWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, fvMesh>& iF,
    const dictionary& dict
)
:
    nutNonEquilibriumWallFunctionFvPatchScalarField(p, iF, dict)
{}


Foam::nutNonEquilibriumSpaldingWallFunctionFvPatchScalarField::
nutNonEquilibriumSpaldingWallFunctionFvPatchScalarField
(
    const nutNonEquilibriumSpaldingWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, fvMesh>& iF,
    const fieldMapper& mapper
)
:
    nutNonEquilibriumWallFunctionFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::nutNonEquilibriumSpaldingWallFunctionFvPatchScalarField::
nutNonEquilibriumSpaldingWallFunctionFvPatchScalarField
(
    const nutNonEquilibriumSpaldingWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, fvMesh>& iF
)
:
    nutNonEquilibriumWallFunctionFvPatchScalarField(wfpsf, iF)
{}


// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::nutNonEquilibriumSpaldingWallFunctionFvPatchScalarField::uTau() const
{
    const label patchi = patch().index();
    const momentumTransportModel& turbModel =
        db().lookupType<momentumTransportModel>(internalField().group());

    const scalarField& y = turbModel.yb()[patchi];
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const scalarField Utilde(UTilde());

    tmp<scalarField> tuTau(new scalarField(patch().size(), 0.0));
    scalarField& ut = tuTau.ref();

    forAll(ut, facei)
    {
        ut[facei] = wallLaw::spaldingUTau
        (
            Utilde[facei], y[facei], nuw[facei], kappa_, E_
        );
    }

    return tuTau;
}


Foam::tmp<Foam::scalarField>
Foam::nutNonEquilibriumSpaldingWallFunctionFvPatchScalarField::yPlus() const
{
    const label patchi = patch().index();
    const momentumTransportModel& turbModel =
        db().lookupType<momentumTransportModel>(internalField().group());

    const scalarField& y = turbModel.yb()[patchi];
    const tmp<scalarField> tnuw = turbModel.nu(patchi);

    return y*uTau()/tnuw();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        nutNonEquilibriumSpaldingWallFunctionFvPatchScalarField
    );
}

// ************************************************************************* //
