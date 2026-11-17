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

    Parts follow RASModels::LRR and RASModels::kOmega2006 of OpenFOAM
    (Copyright (C) 2011-2026 OpenFOAM Foundation).

\*---------------------------------------------------------------------------*/

#include "stressOmega.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "wallFvPatch.H"
#include "fvm.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
void stressOmega<BasicMomentumTransportModel>::boundOmega()
{
    omega_ = max(omega_, k_/(this->nutMaxCoeff_*this->nu()));
}


template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal> stressOmega<BasicMomentumTransportModel>::beta
(
    const volTensorField& gradU
) const
{
    const volSymmTensorField::Internal S(symm(gradU()));
    const volSymmTensorField::Internal Shat(S - 0.5*tr(S)*I);
    const volTensorField::Internal W(skew(gradU.v()));

    const volScalarField::Internal chiOmega
    (
        typedName("chiOmega"),
        mag((W & W) && Shat)/pow3(betaStar_*omega_.v())
    );

    return beta0_*(1 + 85*chiOmega)/(1 + 100*chiOmega);
}


template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal>
stressOmega<BasicMomentumTransportModel>::CDkOmega() const
{
    return max
    (
        sigmaDo_*(fvc::grad(k_)().v() & fvc::grad(omega_)().v())/omega_(),
        dimensionedScalar(dimless/sqr(dimTime), 0)
    );
}


template<class BasicMomentumTransportModel>
void stressOmega<BasicMomentumTransportModel>::correctNut()
{
    this->nut_ = k_/omega_;
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}


template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix>
stressOmega<BasicMomentumTransportModel>::omegaSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
stressOmega<BasicMomentumTransportModel>::stressOmega
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity,
    const word& type
)
:
    ReynoldsStress<RASModel<BasicMomentumTransportModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    ),

    betaStar_("betaStar", this->typeDict(type), 0.09),
    beta0_("beta0", this->typeDict(type), 0.0708),
    gamma_("gamma", this->typeDict(type), 0.52),
    C1_("C1", this->typeDict(type), 1.8),
    C2_("C2", this->typeDict(type), 10.0/19.0),
    sigma_("sigma", this->typeDict(type), 0.5),
    sigmaStar_("sigmaStar", this->typeDict(type), 0.6),
    sigmaDo_("sigmaDo", this->typeDict(type), 0.125),

    k_
    (
        IOobject
        (
            this->groupName("k"),
            this->runTime_.name(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0.5*tr(this->R_)
    ),
    omega_
    (
        IOobject
        (
            this->groupName("omega"),
            this->runTime_.name(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensions::turbulentOmega
    )
{
    if (type == typeName)
    {
        this->boundNormalStress(this->R_);
        k_ = 0.5*tr(this->R_);
        boundOmega();
    }
}


// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool stressOmega<BasicMomentumTransportModel>::read()
{
    if (ReynoldsStress<RASModel<BasicMomentumTransportModel>>::read())
    {
        betaStar_.readIfPresent(this->typeDict());
        beta0_.readIfPresent(this->typeDict());
        gamma_.readIfPresent(this->typeDict());
        C1_.readIfPresent(this->typeDict());
        C2_.readIfPresent(this->typeDict());
        sigma_.readIfPresent(this->typeDict());
        sigmaStar_.readIfPresent(this->typeDict());
        sigmaDo_.readIfPresent(this->typeDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
void stressOmega<BasicMomentumTransportModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volSymmTensorField& R = this->R_;
    const Foam::fvModels& fvModels(Foam::fvModels::New(this->mesh_));
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(this->mesh_)
    );

    ReynoldsStress<RASModel<BasicMomentumTransportModel>>::correct();

    tmp<volTensorField> tgradU(fvc::grad(U));
    const volTensorField& gradU = tgradU();

    // Production tensor, its counterpart with the transposed velocity
    // gradient (LRR-QI), and the production of k
    volSymmTensorField P(-twoSymm(R & gradU));
    // D = -twoSymm(R & gradU^T) is written as -twoSymm(gradU & R) (the same
    // tensor) so that no operand of the inner product is a temporary:
    // R & gradU.T() computes the product in the storage of the temporary
    // gradU.T() and returns wrong off-diagonal components (OpenFOAM-14).
    volSymmTensorField D(-twoSymm(gradU & R));
    volScalarField::Internal G(this->GName(), 0.5*tr(P()));

    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();

    // Specific dissipation rate equation. The production gamma G omega/k is
    // a source where G > 0 and an implicit sink where G < 0.
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(alpha, rho, omega_)
      + fvm::div(alphaRhoPhi, omega_)
      - fvm::laplacian(alpha*rho*DomegaEff(), omega_)
     ==
      - fvm::SuSp(-gamma_*alpha()*rho()*G/k_(), omega_)
      - fvm::Sp(beta(gradU)*alpha()*rho()*omega_(), omega_)
      + alpha()*rho()*CDkOmega()
      + omegaSource()
      + fvModels.source(alpha, rho, omega_)
    );

    omegaEqn.ref().relax();
    fvConstraints.constrain(omegaEqn.ref());
    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
    solve(omegaEqn);
    fvConstraints.constrain(omega_);
    boundOmega();


    // Correct the trace of the tensorial production to be consistent
    // with the near-wall generation from the wall functions (as LRR)
    const fvPatchList& patches = this->mesh_.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isA<wallFvPatch>(curPatch))
        {
            forAll(curPatch, facei)
            {
                const label celli = curPatch.faceCells()[facei];
                const scalar f = min
                (
                    mag(G[celli])/(0.5*mag(tr(P[celli])) + small),
                    1.0
                );
                P[celli] *= f;
                D[celli] *= f;
            }
        }
    }

    // LRR-QI pressure-strain coefficients
    const dimensionedScalar alphaHat((8 + C2_)/11);
    const dimensionedScalar betaHat((8*C2_ - 2)/11);
    const dimensionedScalar gammaHat((60*C2_ - 4)/55);

    // Reynolds stress equation
    tmp<fvSymmTensorMatrix> REqn
    (
        fvm::ddt(alpha, rho, R)
      + fvm::div(alphaRhoPhi, R)
      - fvm::laplacian(alpha*rho*DREff(), R)
      + fvm::Sp(betaStar_*C1_*alpha*rho*omega_, R)
     ==
        alpha*rho*P
      + ((2.0/3.0)*betaStar_*(C1_ - 1)*I)*alpha*rho*omega_*k_
      - alphaHat*alpha*rho*dev(P)
      - betaHat*alpha*rho*dev(D)
      - (gammaHat*alpha*rho*k_)*dev(symm(gradU))
      + this->RSource()
      + fvModels.source(alpha, rho, R)
    );

    REqn.ref().relax();
    fvConstraints.constrain(REqn.ref());
    solve(REqn);
    fvConstraints.constrain(R);

    this->boundNormalStress(R);

    k_ = 0.5*tr(R);
    boundOmega();

    correctNut();

    // Correct wall shear-stresses when applying wall-functions
    this->correctWallShearStress(R);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
