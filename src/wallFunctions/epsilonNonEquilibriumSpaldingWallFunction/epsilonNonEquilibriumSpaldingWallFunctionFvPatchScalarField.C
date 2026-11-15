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

#include "epsilonNonEquilibriumSpaldingWallFunctionFvPatchScalarField.H"
#include "nutNonEquilibriumWallFunctionFvPatchScalarField.H"
#include "wallLaw.H"
#include "momentumTransportModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::scalarField>>
Foam::epsilonNonEquilibriumSpaldingWallFunctionFvPatchScalarField::calculate
(
    const momentumTransportModel& mtm
) const
{
    const label patchi = patch().index();

    const nutWallFunctionFvPatchScalarField& nutw =
        nutWallFunctionFvPatchScalarField::nutw(mtm, patchi);

    const scalarField& y = mtm.yb()[patchi];

    const tmp<scalarField> tnuw = mtm.nu(patchi);
    const scalarField& nuw = tnuw();

    const tmp<volScalarField> tk = mtm.k();
    const volScalarField& k = tk();

    const fvPatchVectorField& Uw = mtm.U().boundaryField()[patchi];

    const scalarField magGradUw(mag(Uw.snGrad()));

    // Pressure gradient along the flow, if the nut wall function has it
    const scalarField dpds
    (
        isA<nutNonEquilibriumWallFunctionFvPatchScalarField>(nutw)
      ? refCast<const nutNonEquilibriumWallFunctionFvPatchScalarField>
        (nutw).dpds()
      : tmp<scalarField>(new scalarField(patch().size(), 0.0))
    );

    const scalar Cmu = nutw.Cmu();
    const scalar Cmu25 = pow025(Cmu);
    const scalar kappa = nutw.kappa();
    const scalar E = nutw.E();

    Pair<tmp<scalarField>> GandValue
    (
        tmp<scalarField>(new scalarField(patch().size())),
        tmp<scalarField>(new scalarField(patch().size()))
    );
    scalarField& G = GandValue.first().ref();
    scalarField& value = GandValue.second().ref();

    forAll(nutw, facei)
    {
        const label celli = patch().faceCells()[facei];
        const scalar kc = max(k[celli], small*small);

        const scalar tauw =
            max((nutw[facei] + nuw[facei])*magGradUw[facei], 0);
        const scalar ut = sqrt(tauw);

        if (ut > vSmall)
        {
            const scalar yPlus = ut*y[facei]/nuw[facei];
            const scalar Uplus = wallLaw::spaldingUPlus(yPlus, kappa, E);
            const scalar gamma =
                1/wallLaw::spaldingDYPlusDUPlus(Uplus, kappa, E);
            const scalar pPlus = nuw[facei]*dpds[facei]/pow3(ut);
            const scalar tauPlus = max(1 + pPlus*yPlus, 0);

            G[facei] = pow4(ut)/nuw[facei]*max(tauPlus - gamma, 0)*gamma;
        }
        else
        {
            G[facei] = 0;
        }

        {
            const scalar uk = Cmu25*sqrt(kc);
            const scalar yStar = uk*y[facei]/nuw[facei];
            const scalar Gamma = wallLaw::blendingGamma(yStar);
            const scalar wVis = exp(-Gamma);
            const scalar wLog = Gamma > small ? exp(-1/Gamma) : 0;

            value[facei] =
                wVis*2*kc*nuw[facei]/sqr(y[facei])
              + wLog*pow3(uk)/(kappa*y[facei]);
        }
    }

    return GandValue;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        epsilonNonEquilibriumSpaldingWallFunctionFvPatchScalarField
    );
}

// ************************************************************************* //
