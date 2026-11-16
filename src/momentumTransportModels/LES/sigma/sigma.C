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

#include "sigma.H"
#include "sigmaOperator.H"
#include "fvConstraints.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
tmp<volScalarField> sigma<BasicMomentumTransportModel>::Dsigma
(
    const volTensorField& gradU
) const
{
    tmp<volScalarField> tD
    (
        volScalarField::New
        (
            this->groupName("Dsigma"),
            this->mesh_,
            dimensionedScalar(dimless/dimTime, 0)
        )
    );
    volScalarField& D = tD.ref();

    scalarField& Di = D.primitiveFieldRef();
    const tensorField& gi = gradU.primitiveField();
    forAll(Di, celli)
    {
        Di[celli] = sigmaOperator::Dsigma(gi[celli]);
    }

    volScalarField::Boundary& Db = D.boundaryFieldRef();
    forAll(Db, patchi)
    {
        scalarField& Dp = Db[patchi];
        const tensorField& gp = gradU.boundaryField()[patchi];
        forAll(Dp, facei)
        {
            Dp[facei] = sigmaOperator::Dsigma(gp[facei]);
        }
    }

    return tD;
}


template<class BasicMomentumTransportModel>
void sigma<BasicMomentumTransportModel>::correctNut()
{
    this->nut_ =
        sqr(Csigma_*this->delta())*Dsigma(fvc::grad(this->U_));
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
sigma<BasicMomentumTransportModel>::sigma
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
    LESeddyViscosity<BasicMomentumTransportModel>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    ),

    Ck_("Ck", this->typeDict(type), 0.094),
    Csigma_("Csigma", this->typeDict(type), 1.35)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool sigma<BasicMomentumTransportModel>::read()
{
    if (LESeddyViscosity<BasicMomentumTransportModel>::read())
    {
        Ck_.readIfPresent(this->typeDict());
        Csigma_.readIfPresent(this->typeDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
tmp<volScalarField> sigma<BasicMomentumTransportModel>::k() const
{
    return volScalarField::New
    (
        this->groupName("k"),
        sqr(this->nut_/(Ck_*this->delta()))
    );
}


template<class BasicMomentumTransportModel>
tmp<volScalarField> sigma<BasicMomentumTransportModel>::epsilon() const
{
    const volScalarField k(this->k());

    return volScalarField::New
    (
        this->groupName("epsilon"),
        this->Ce_*k*sqrt(k)/this->delta()
    );
}


template<class BasicMomentumTransportModel>
void sigma<BasicMomentumTransportModel>::correct()
{
    LESeddyViscosity<BasicMomentumTransportModel>::correct();
    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// ************************************************************************* //
