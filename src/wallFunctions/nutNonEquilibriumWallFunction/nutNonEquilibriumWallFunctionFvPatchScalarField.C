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

#include "nutNonEquilibriumWallFunctionFvPatchScalarField.H"
#include "wallLaw.H"
#include "momentumTransportModel.H"
#include "fieldMapper.H"
#include "volFields.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::nutNonEquilibriumWallFunctionFvPatchScalarField::nut() const
{
    const label patchi = patch().index();
    const momentumTransportModel& turbModel =
        db().lookupType<momentumTransportModel>(internalField().group());

    const scalarField& y = turbModel.yb()[patchi];
    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const scalarField Ut(magUt());
    const scalarField Utilde(UTilde());

    const scalar Cmu25 = pow025(Cmu_);

    tmp<scalarField> tnutw(new scalarField(patch().size(), 0.0));
    scalarField& nutw = tnutw.ref();

    forAll(nutw, facei)
    {
        if (Ut[facei] < vSmall)
        {
            continue;
        }

        const label celli = patch().faceCells()[facei];
        const scalar yStar = Cmu25*y[facei]*sqrt(k[celli])/nuw[facei];
        const scalar ratio = Utilde[facei]/Ut[facei];

        if (yStar > yPlusLam_)
        {
            nutw[facei] =
                nuw[facei]*(yStar*kappa_/log(E_*yStar)*ratio - 1);
        }
        // else viscous wall cell: tau_w = nu U/y, nut_w = 0 (as Fluent)

        nutw[facei] = max(nutw[facei], -0.999*nuw[facei]);
    }

    return tnutw;
}


void Foam::nutNonEquilibriumWallFunctionFvPatchScalarField::writeLocalEntries
(
    Ostream& os
) const
{
    nutWallFunctionFvPatchScalarField::writeLocalEntries(os);
    writeEntry(os, "yvStar", yvStar_);
    writeEntryIfDifferent<word>(os, "p", "p", pName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nutNonEquilibriumWallFunctionFvPatchScalarField::
nutNonEquilibriumWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, fvMesh>& iF,
    const dictionary& dict
)
:
    nutWallFunctionFvPatchScalarField(p, iF, dict),
    yvStar_(dict.lookupOrDefault<scalar>("yvStar", 11.225)),
    pName_(dict.lookupOrDefault<word>("p", "p")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho"))
{}


Foam::nutNonEquilibriumWallFunctionFvPatchScalarField::
nutNonEquilibriumWallFunctionFvPatchScalarField
(
    const nutNonEquilibriumWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, fvMesh>& iF,
    const fieldMapper& mapper
)
:
    nutWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    yvStar_(ptf.yvStar_),
    pName_(ptf.pName_),
    rhoName_(ptf.rhoName_)
{}


Foam::nutNonEquilibriumWallFunctionFvPatchScalarField::
nutNonEquilibriumWallFunctionFvPatchScalarField
(
    const nutNonEquilibriumWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, fvMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(wfpsf, iF),
    yvStar_(wfpsf.yvStar_),
    pName_(wfpsf.pName_),
    rhoName_(wfpsf.rhoName_)
{}


// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::nutNonEquilibriumWallFunctionFvPatchScalarField::magUt() const
{
    const label patchi = patch().index();
    const momentumTransportModel& turbModel =
        db().lookupType<momentumTransportModel>(internalField().group());

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const vectorField Ur(Uw.patchInternalField() - Uw);
    const vectorField n(patch().nf());

    return mag(Ur - n*(n & Ur));
}


Foam::tmp<Foam::scalarField>
Foam::nutNonEquilibriumWallFunctionFvPatchScalarField::dpds() const
{
    const label patchi = patch().index();

    tmp<scalarField> tdpds(new scalarField(patch().size(), 0.0));

    if (!db().foundObject<volScalarField>(pName_))
    {
        return tdpds;
    }

    const momentumTransportModel& turbModel =
        db().lookupType<momentumTransportModel>(internalField().group());

    const volScalarField& p = db().lookupObject<volScalarField>(pName_);
    vectorField gradp
    (
        fvc::grad(p)().boundaryField()[patchi].patchInternalField()
    );

    // Pressure in Pa: divide by the density of the wall cells
    if (p.dimensions() == dimPressure)
    {
        if (!db().foundObject<volScalarField>(rhoName_))
        {
            FatalErrorInFunction
                << "Pressure field " << pName_ << " has the dimensions of "
                << "pressure but there is no density field " << rhoName_
                << " to make it kinematic" << exit(FatalError);
        }

        const volScalarField& rho =
            db().lookupObject<volScalarField>(rhoName_);
        gradp /= rho.boundaryField()[patchi].patchInternalField();
    }

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const vectorField Ur(Uw.patchInternalField() - Uw);
    const vectorField n(patch().nf());
    const vectorField Ut(Ur - n*(n & Ur));

    scalarField& dps = tdpds.ref();
    forAll(dps, facei)
    {
        const scalar magUt = mag(Ut[facei]);
        if (magUt > vSmall)
        {
            dps[facei] = gradp[facei] & (Ut[facei]/magUt);
        }
    }

    return tdpds;
}


Foam::tmp<Foam::scalarField>
Foam::nutNonEquilibriumWallFunctionFvPatchScalarField::UTilde() const
{
    const label patchi = patch().index();
    const momentumTransportModel& turbModel =
        db().lookupType<momentumTransportModel>(internalField().group());

    const scalarField& y = turbModel.yb()[patchi];
    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const scalarField Ut(magUt());
    const scalarField dps(dpds());

    tmp<scalarField> tUtilde(new scalarField(patch().size()));
    scalarField& Utilde = tUtilde.ref();

    forAll(Utilde, facei)
    {
        const label celli = patch().faceCells()[facei];
        Utilde[facei] = wallLaw::UTilde
        (
            Ut[facei], dps[facei], y[facei], k[celli], nuw[facei],
            Cmu_, kappa_, yvStar_
        );
    }

    return tUtilde;
}


Foam::tmp<Foam::scalarField>
Foam::nutNonEquilibriumWallFunctionFvPatchScalarField::yPlus() const
{
    const label patchi = patch().index();
    const momentumTransportModel& turbModel =
        db().lookupType<momentumTransportModel>(internalField().group());

    const scalarField& y = turbModel.yb()[patchi];
    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();
    tmp<scalarField> kwc = k.boundaryField()[patchi].patchInternalField();
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    return pow025(Cmu_)*y*sqrt(kwc)/nuw;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        nutNonEquilibriumWallFunctionFvPatchScalarField
    );
}

// ************************************************************************* //
