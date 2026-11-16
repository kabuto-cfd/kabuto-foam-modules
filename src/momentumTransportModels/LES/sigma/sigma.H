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

    The structure follows LESModels::WALE of OpenFOAM
    (Copyright (C) 2011-2026 OpenFOAM Foundation).

Class
    Foam::LESModels::sigma

Description
    The sigma subgrid-scale model of Nicoud et al. (2011): an eddy viscosity
    built from the singular values sigma1 >= sigma2 >= sigma3 of the resolved
    velocity gradient,

    \verbatim
        nut     = (Csigma delta)^2 D_sigma
        D_sigma = sigma3 (sigma1 - sigma2)(sigma2 - sigma3)/sigma1^2
        k       = (nut/(Ck delta))^2
        epsilon = Ce k^(3/2)/delta
    \endverbatim

    (sigmaOperator.H). Unlike Smagorinsky's model it vanishes in pure shear
    and decays as y^3 at a wall without a damping function; unlike WALE it
    also vanishes in solid-body rotation, for which WALE returns an eddy
    viscosity proportional to the rotation rate (in a forced vortex, the core
    of a cyclone or a swirl tube). It vanishes too for isotropic and
    axisymmetric expansion or contraction.

    Csigma = 1.35 is the value of the reference, derived from Smagorinsky's
    constant 0.18 for isotropic turbulence. Some implementations write
    nut = C delta^2 D_sigma with C = 1.68 (calibrated on decaying isotropic
    turbulence), which is Csigma = sqrt(1.68) = 1.30 here.

    Reference:
    \verbatim
        Nicoud, F., Baya Toda, H., Cabrit, O., Bose, S., Lee, J. (2011).
        Using singular values to build a subgrid-scale model for large eddy
        simulations.
        Physics of Fluids, 23(8), 085106.
    \endverbatim

Usage
    \table
        Property | Description                          | Required | Default
        Csigma   | model coefficient                    | no       | 1.35
        Ck       | k coefficient (k and epsilon only)   | no       | 0.094
        Ce       | epsilon coefficient (LESModel)       | no       | 1.048
    \endtable

    In constant/momentumTransport, with the library loaded by
    libs ("libkabutoMomentumTransportModels.so"); in system/controlDict:
    \verbatim
    simulationType  LES;

    LES
    {
        model           sigma;
        delta           cubeRootVol;
        cubeRootVol
        {
            deltaCoeff      1;
        }
    }
    \endverbatim

See also
    Foam::LESModels::WALE
    Foam::LESModels::Smagorinsky

SourceFiles
    sigma.C

\*---------------------------------------------------------------------------*/

#ifndef sigma_H
#define sigma_H

#include "LESModel.H"
#include "LESeddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{

/*---------------------------------------------------------------------------*\
                            Class sigma Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMomentumTransportModel>
class sigma
:
    public LESeddyViscosity<BasicMomentumTransportModel>
{
protected:

    // Protected data

        dimensionedScalar Ck_;
        dimensionedScalar Csigma_;


    // Protected Member Functions

        //- The differential operator D_sigma of the given velocity gradient
        tmp<volScalarField> Dsigma(const volTensorField& gradU) const;

        //- Update the SGS eddy-viscosity
        virtual void correctNut();


public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;


    //- Runtime type information
    TypeName("sigma");


    // Constructors

        //- Construct from components
        sigma
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
        sigma(const sigma&) = delete;


    //- Destructor
    virtual ~sigma()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Return SGS kinetic energy
        virtual tmp<volScalarField> k() const;

        //- Return sub-grid dissipation rate
        virtual tmp<volScalarField> epsilon() const;

        //- Correct Eddy-Viscosity and related properties
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const sigma&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sigma.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
