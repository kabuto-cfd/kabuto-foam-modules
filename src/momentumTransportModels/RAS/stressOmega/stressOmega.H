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

    The structure follows RASModels::LRR and RASModels::kOmega2006 of
    OpenFOAM (Copyright (C) 2011-2026 OpenFOAM Foundation).

Class
    Foam::RASModels::stressOmega

Description
    Wilcox's stress-omega Reynolds-stress model (2006) for incompressible and
    compressible flows: six transport equations for the Reynolds stresses with
    the LRR-QI pressure-strain model, closed with the omega equation of the
    k-omega (2006) model.

    Reference:
    \verbatim
        Wilcox, D. C. (2006).
        Turbulence modeling for CFD, 3rd edition, section 6.3.3.
        La Canada, CA: DCW Industries, Inc.
    \endverbatim

    Written for the Reynolds stress R = <u'u'> (Wilcox's tau_ij = -R_ij),
    with P = -twoSymm(R & grad(U)), D = -twoSymm(R & grad(U)^T),
    G = tr(P)/2 (production of k), S = symm(grad(U)), W = skew(grad(U)),
    k = tr(R)/2 and nut = k/omega:

    \verbatim
        D(R)/Dt = P + Pi - 2/3 betaStar omega k I
                + div((nu + sigmaStar nut) grad(R))

        Pi = betaStar C1 omega (2/3 k I - R)
           - alphaHat dev(P) - betaHat dev(D) - gammaHat k dev(S)

        D(omega)/Dt = gamma G omega/k - beta omega^2
                    + sigmad/omega grad(k) & grad(omega)
                    + div((nu + sigma nut) grad(omega))

        beta   = beta0 fBeta,  fBeta = (1 + 85 chi)/(1 + 100 chi),
        chi    = |(W & W) && (S - 1/2 tr(S) I)|/(betaStar omega)^3,
        sigmad = sigmaDo where grad(k) & grad(omega) > 0, else 0
    \endverbatim

    The pressure-strain coefficients are those of LRR-QI written with C2:
    alphaHat = (8 + C2)/11, betaHat = (8 C2 - 2)/11,
    gammaHat = (60 C2 - 4)/55. Use omegaWallFunction / kqRWallFunction /
    nutkWallFunction for wall functions, or R = 0 and a large omega (Wilcox's
    or Menter's wall value) at walls resolved down to y+ ~ 1.

    The default model coefficients are
    \verbatim
        stressOmegaCoeffs
        {
            betaStar    0.09;
            beta0       0.0708;
            gamma       0.52;       // Wilcox's alpha = 13/25
            C1          1.8;
            C2          0.526316;   // 10/19
            sigma       0.5;
            sigmaStar   0.6;
            sigmaDo     0.125;

            couplingFactor 0.0;
        }
    \endverbatim

    k is computed from R and written; omega must be given (0/omega).

SourceFiles
    stressOmega.C

\*---------------------------------------------------------------------------*/

#ifndef stressOmega_H
#define stressOmega_H

#include "RASModel.H"
#include "ReynoldsStress.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                         Class stressOmega Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class stressOmega
:
    public ReynoldsStress<RASModel<BasicMomentumTransportModel>>
{
protected:

    // Protected data

        // Model coefficients

            dimensionedScalar betaStar_;
            dimensionedScalar beta0_;
            dimensionedScalar gamma_;
            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar sigma_;
            dimensionedScalar sigmaStar_;
            dimensionedScalar sigmaDo_;


        // Fields

            volScalarField k_;
            volScalarField omega_;


    // Protected Member Functions

        //- Bound omega so that nut does not exceed nutMaxCoeff*nu
        void boundOmega();

        //- Vortex-stretching factor times beta0
        tmp<volScalarField::Internal> beta(const volTensorField& gradU) const;

        //- Cross-diffusion term sigmad/omega grad(k) & grad(omega)
        tmp<volScalarField::Internal> CDkOmega() const;

        //- Update the eddy viscosity nut = k/omega
        virtual void correctNut();

        //- Source term for the omega equation
        virtual tmp<fvScalarMatrix> omegaSource() const;


public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;


    //- Runtime type information
    TypeName("stressOmega");


    // Constructors

        //- Construct from components
        stressOmega
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const viscosity& viscosity,
            const word& type = typeName
        );

        //- Disallow default bitwise copy construction
        stressOmega(const stressOmega&) = delete;


    //- Destructor
    virtual ~stressOmega()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return volScalarField::New
            (
                this->groupName("epsilon"),
                betaStar_*k_*omega_
            );
        }

        //- Return the turbulence specific dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }

        //- Return the effective diffusivity for R
        tmp<volScalarField> DREff() const
        {
            return volScalarField::New
            (
                "DREff",
                sigmaStar_*this->nut_ + this->nu()
            );
        }

        //- Return the effective diffusivity for omega
        tmp<volScalarField> DomegaEff() const
        {
            return volScalarField::New
            (
                "DomegaEff",
                sigma_*this->nut_ + this->nu()
            );
        }

        //- Solve the turbulence equations and correct eddy viscosity
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const stressOmega&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "stressOmega.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
