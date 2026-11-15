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

#include "omegaNonEquilibriumWallFunctionFvPatchScalarField.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "nutNonEquilibriumWallFunctionFvPatchScalarField.H"
#include "momentumTransportModel.H"
#include "fvMatrix.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::scalarField>>
Foam::omegaNonEquilibriumWallFunctionFvPatchScalarField::calculate
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

    const VolInternalField<scalar>& Gvf =
        db().lookupObject<VolInternalField<scalar>>(mtm.GName());

    const scalar Cmu = nutw.Cmu();
    const scalar Cmu25 = pow025(Cmu);
    const scalar Cmu75 = pow(Cmu, 0.75);
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
            const scalar epsilon =
                kc/yn*(2*nuw[facei]/yv + Cmu75*sqrtk/kappa*lnr);
            value[facei] = epsilon/(Cmu*kc);
        }
        else
        {
            G[facei] = Gvf[celli];
            value[facei] = 6*nuw[facei]/(beta1_*sqr(y[facei]));
        }
    }

    return GandValue;
}


void Foam::omegaNonEquilibriumWallFunctionFvPatchScalarField::updateCoeffsMaster()
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
    volScalarField& omega =
        const_cast<volScalarField&>
        (
            static_cast<const volScalarField&>
            (
                internalField()
            )
        );

    const volScalarField::Boundary& omegaBf = omega.boundaryField();

    // Make all processors build the near wall distances
    mtm.yb();

    // Evaluate all the wall functions
    PtrList<scalarField> Gpfs(omegaBf.size());
    PtrList<scalarField> omegaPfs(omegaBf.size());
    forAll(omegaBf, patchi)
    {
        if (isA<omegaNonEquilibriumWallFunctionFvPatchScalarField>(omegaBf[patchi]))
        {
            const omegaNonEquilibriumWallFunctionFvPatchScalarField& omegaPf =
                refCast<const omegaNonEquilibriumWallFunctionFvPatchScalarField>
                (omegaBf[patchi]);

            Pair<tmp<scalarField>> GandOmega = omegaPf.calculate(mtm);

            Gpfs.set(patchi, GandOmega.first().ptr());
            omegaPfs.set(patchi, GandOmega.second().ptr());
        }
    }

    // Average the values into the wall-adjacent cells and store
    wallCellGPtr_.reset(patchFieldsToWallCellField(Gpfs).ptr());
    wallCellOmegaPtr_.reset(patchFieldsToWallCellField(omegaPfs).ptr());

    // Set the fractional proportion of the value in the wall cells
    UIndirectList<scalar>(G, wallCells()) =
        (1 - wallCellFraction())*scalarField(G, wallCells())
      + wallCellFraction()*wallCellGPtr_();
    UIndirectList<scalar>(omega, wallCells()) =
        (1 - wallCellFraction())*scalarField(omega, wallCells())
      + wallCellFraction()*wallCellOmegaPtr_();
}


void Foam::omegaNonEquilibriumWallFunctionFvPatchScalarField::manipulateMatrixMaster
(
    fvMatrix<scalar>& matrix
)
{
    if (patch().index() != masterPatchIndex())
    {
        return;
    }

    matrix.setValues(wallCells(), wallCellOmegaPtr_(), wallCellFraction());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::omegaNonEquilibriumWallFunctionFvPatchScalarField::
omegaNonEquilibriumWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, fvMesh>& iF,
    const dictionary& dict
)
:
    wallCellWallFunctionFvPatchScalarField(p, iF, dict),
    wallCellGPtr_(nullptr),
    wallCellOmegaPtr_(nullptr),
    beta1_(dict.lookupOrDefault<scalar>("beta1", 0.075))
{}


Foam::omegaNonEquilibriumWallFunctionFvPatchScalarField::
omegaNonEquilibriumWallFunctionFvPatchScalarField
(
    const omegaNonEquilibriumWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, fvMesh>& iF,
    const fieldMapper& mapper
)
:
    wallCellWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    wallCellGPtr_(nullptr),
    wallCellOmegaPtr_(nullptr),
    beta1_(ptf.beta1_)
{}


Foam::omegaNonEquilibriumWallFunctionFvPatchScalarField::
omegaNonEquilibriumWallFunctionFvPatchScalarField
(
    const omegaNonEquilibriumWallFunctionFvPatchScalarField& ewfpsf,
    const DimensionedField<scalar, fvMesh>& iF
)
:
    wallCellWallFunctionFvPatchScalarField(ewfpsf, iF),
    wallCellGPtr_(nullptr),
    wallCellOmegaPtr_(nullptr),
    beta1_(ewfpsf.beta1_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::omegaNonEquilibriumWallFunctionFvPatchScalarField::map
(
    const fvPatchScalarField& ptf,
    const fieldMapper& mapper
)
{
    wallCellWallFunctionFvPatchScalarField::map(ptf, mapper);
    wallCellGPtr_.clear();
    wallCellOmegaPtr_.clear();
}


void Foam::omegaNonEquilibriumWallFunctionFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    wallCellWallFunctionFvPatchScalarField::reset(ptf);
    wallCellGPtr_.clear();
    wallCellOmegaPtr_.clear();
}


void Foam::omegaNonEquilibriumWallFunctionFvPatchScalarField::updateCoeffs()
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


void Foam::omegaNonEquilibriumWallFunctionFvPatchScalarField::manipulateMatrix
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
void Foam::omegaNonEquilibriumWallFunctionFvPatchScalarField::write
(
    Ostream& os
) const
{
    wallCellWallFunctionFvPatchScalarField::write(os);
    writeEntry(os, "beta1", beta1_);
}


namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        omegaNonEquilibriumWallFunctionFvPatchScalarField
    );
}


// ************************************************************************* //
