/*---------------------------------------------------------------------------*\
  Kabuto CFD modules for OpenFOAM-14
  Copyright (C) 2026 Kabuto CFD
  Derived from epsilonWallFunctionFvPatchScalarField of OpenFOAM
  Copyright (C) 2011-2026 OpenFOAM Foundation
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

#include "epsilonNonEquilibriumWallFunctionFvPatchScalarField.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "nutNonEquilibriumWallFunctionFvPatchScalarField.H"
#include "momentumTransportModel.H"
#include "fvMatrix.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::scalarField>>
Foam::epsilonNonEquilibriumWallFunctionFvPatchScalarField::calculate
(
    const momentumTransportModel& mtm
) const
{
    const label patchi = patch().index();

    const nutWallFunctionFvPatchScalarField& nutw =
        nutWallFunctionFvPatchScalarField::nutw(mtm, patchi);

    const scalar yvStar =
        isA<nutNonEquilibriumWallFunctionFvPatchScalarField>(nutw)
      ? refCast<const nutNonEquilibriumWallFunctionFvPatchScalarField>
        (nutw).yvStar()
      : 11.225;

    const scalarField& y = mtm.yb()[patchi];

    const tmp<scalarField> tnuw = mtm.nu(patchi);
    const scalarField& nuw = tnuw();

    const tmp<volScalarField> tk = mtm.k();
    const volScalarField& k = tk();

    const fvPatchVectorField& Uw = mtm.U().boundaryField()[patchi];

    const scalarField magGradUw(mag(Uw.snGrad()));

    const scalar Cmu25 = pow025(nutw.Cmu());
    const scalar Cmu75 = pow(nutw.Cmu(), 0.75);
    const scalar kappa = nutw.kappa();

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
        const scalar sqrtk = sqrt(kc);
        const scalar uk = Cmu25*sqrtk;
        const scalar yStar = uk*y[facei]/nuw[facei];

        if (yStar > nutw.yPlusLam())
        {
            const scalar tauw = (nutw[facei] + nuw[facei])*magGradUw[facei];
            const scalar yn = 2*y[facei];
            const scalar yv = yvStar*nuw[facei]/uk;
            const scalar lnr = log(max(yn/yv, 1.0));

            G[facei] = sqr(tauw)/(kappa*uk*yn)*lnr;
            value[facei] = kc/yn*(2*nuw[facei]/yv + Cmu75*sqrtk/kappa*lnr);
        }
        else
        {
            G[facei] = 0;
            value[facei] = 2*kc*nuw[facei]/sqr(y[facei]);
        }
    }

    return GandValue;
}


void Foam::epsilonNonEquilibriumWallFunctionFvPatchScalarField::updateCoeffsMaster()
{
    if (patch().index() != masterPatchIndex())
    {
        return;
    }

    // Lookup the momentum transport model
    const momentumTransportModel& mtm =
        db().lookupType<momentumTransportModel>(internalField().group());

    // Get mutable references to the turbulence fields
    VolInternalField<scalar>& G =
        db().lookupObjectRef<VolInternalField<scalar>>(mtm.GName());
    volScalarField& epsilon =
        const_cast<volScalarField&>
        (
            static_cast<const volScalarField&>
            (
                internalField()
            )
        );

    const volScalarField::Boundary& epsilonBf = epsilon.boundaryField();

    // Make all processors build the near wall distances
    mtm.yb();

    // Evaluate all the wall functions
    PtrList<scalarField> Gpfs(epsilonBf.size());
    PtrList<scalarField> epsilonPfs(epsilonBf.size());
    forAll(epsilonBf, patchi)
    {
        if (isA<epsilonNonEquilibriumWallFunctionFvPatchScalarField>(epsilonBf[patchi]))
        {
            const epsilonNonEquilibriumWallFunctionFvPatchScalarField& epsilonPf =
                refCast<const epsilonNonEquilibriumWallFunctionFvPatchScalarField>
                (epsilonBf[patchi]);

            Pair<tmp<scalarField>> GandEpsilon = epsilonPf.calculate(mtm);

            Gpfs.set(patchi, GandEpsilon.first().ptr());
            epsilonPfs.set(patchi, GandEpsilon.second().ptr());
        }
    }

    // Average the values into the wall-adjacent cells and store
    wallCellGPtr_.reset(patchFieldsToWallCellField(Gpfs).ptr());
    wallCellEpsilonPtr_.reset(patchFieldsToWallCellField(epsilonPfs).ptr());

    // Set the fractional proportion of the value in the wall cells
    UIndirectList<scalar>(G, wallCells()) =
        (1 - wallCellFraction())*scalarField(G, wallCells())
      + wallCellFraction()*wallCellGPtr_();
    UIndirectList<scalar>(epsilon, wallCells()) =
        (1 - wallCellFraction())*scalarField(epsilon, wallCells())
      + wallCellFraction()*wallCellEpsilonPtr_();
}


void Foam::epsilonNonEquilibriumWallFunctionFvPatchScalarField::manipulateMatrixMaster
(
    fvMatrix<scalar>& matrix
)
{
    if (patch().index() != masterPatchIndex())
    {
        return;
    }

    matrix.setValues(wallCells(), wallCellEpsilonPtr_(), wallCellFraction());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::epsilonNonEquilibriumWallFunctionFvPatchScalarField::
epsilonNonEquilibriumWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, fvMesh>& iF,
    const dictionary& dict
)
:
    wallCellWallFunctionFvPatchScalarField(p, iF, dict),
    wallCellGPtr_(nullptr),
    wallCellEpsilonPtr_(nullptr)
{}


Foam::epsilonNonEquilibriumWallFunctionFvPatchScalarField::
epsilonNonEquilibriumWallFunctionFvPatchScalarField
(
    const epsilonNonEquilibriumWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, fvMesh>& iF,
    const fieldMapper& mapper
)
:
    wallCellWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    wallCellGPtr_(nullptr),
    wallCellEpsilonPtr_(nullptr)
{}


Foam::epsilonNonEquilibriumWallFunctionFvPatchScalarField::
epsilonNonEquilibriumWallFunctionFvPatchScalarField
(
    const epsilonNonEquilibriumWallFunctionFvPatchScalarField& ewfpsf,
    const DimensionedField<scalar, fvMesh>& iF
)
:
    wallCellWallFunctionFvPatchScalarField(ewfpsf, iF),
    wallCellGPtr_(nullptr),
    wallCellEpsilonPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::epsilonNonEquilibriumWallFunctionFvPatchScalarField::map
(
    const fvPatchScalarField& ptf,
    const fieldMapper& mapper
)
{
    wallCellWallFunctionFvPatchScalarField::map(ptf, mapper);
    wallCellGPtr_.clear();
    wallCellEpsilonPtr_.clear();
}


void Foam::epsilonNonEquilibriumWallFunctionFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    wallCellWallFunctionFvPatchScalarField::reset(ptf);
    wallCellGPtr_.clear();
    wallCellEpsilonPtr_.clear();
}


void Foam::epsilonNonEquilibriumWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    initMaster();

    updateCoeffsMaster();

    operator==(patchInternalField());

    fvPatchScalarField::updateCoeffs();
}


void Foam::epsilonNonEquilibriumWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix
)
{
    if (manipulatedMatrix())
    {
        return;
    }

    if (masterPatchIndex() == -1)
    {
        FatalErrorInFunction
            << "updateCoeffs must be called before manipulateMatrix"
            << exit(FatalError);
    }

    manipulateMatrixMaster(matrix);

    fvPatchScalarField::manipulateMatrix(matrix);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        epsilonNonEquilibriumWallFunctionFvPatchScalarField
    );
}


// ************************************************************************* //
