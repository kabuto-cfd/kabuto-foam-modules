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

Namespace
    Foam::wallLaw

Description
    Point-wise laws of the wall used by the non-equilibrium wall functions.
    All quantities are kinematic (pressure gradient divided by density).

    UTilde: Kim and Choudhury's (1995) pressure-gradient-sensitised velocity,
    in the form of the Fluent theory guide,
    \verbatim
        y_v    = yvStar nu/(Cmu^1/4 k^1/2)                 (yvStar = 11.225)
        UTilde = U - 1/2 dp/ds [ y_v/(kappa k^1/2) ln(y/y_v)
                               + (y - y_v)/(kappa k^1/2) + y_v^2/nu ]   y > y_v
        UTilde = U - 1/2 dp/ds y^2/nu                                  y <= y_v
    \endverbatim
    The second line is the exact laminar profile with a pressure gradient
    (tau_w = nu UTilde/y); the two meet at y = y_v.

    The correction is limited to |UTilde - U| <= U (0 <= UTilde <= 2 U): it
    vanishes with the velocity, so a wall cell near a separation or
    reattachment line, where U -> 0 and its direction turns, gets no wall
    shear from the pressure gradient alone (unlimited, the shear along a
    vanishing velocity fed back on itself and the pitzDaily recirculation
    did not close).

    Spalding's law (1961), continuous from the viscous sublayer through the
    log layer, with e^(-kappa B) = 1/E:
    \verbatim
        y+ = U+ + 1/E [ e^(kappa U+) - 1 - kappa U+ - (kappa U+)^2/2
                       - (kappa U+)^3/6 ]
    \endverbatim

SourceFiles
    wallLaw.H (inline only)

\*---------------------------------------------------------------------------*/

#ifndef wallLaw_H
#define wallLaw_H

#include "scalar.H"
#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace wallLaw
{

//- Kim and Choudhury's pressure-gradient-sensitised velocity
inline scalar UTilde
(
    const scalar U,         // tangential velocity at the cell centre
    const scalar dpds,      // kinematic pressure gradient along U
    const scalar y,         // wall distance of the cell centre
    const scalar k,         // turbulent kinetic energy at the cell centre
    const scalar nu,        // kinematic viscosity
    const scalar Cmu,
    const scalar kappa,
    const scalar yvStar     // viscous sublayer thickness in y*, 11.225
)
{
    const scalar sqrtk = Foam::sqrt(max(k, 0));
    const scalar uk = Foam::pow(Cmu, 0.25)*sqrtk;
    const scalar yv = uk > vSmall ? yvStar*nu/uk : great;

    scalar dU = 0;
    if (y <= yv)
    {
        dU = -0.5*dpds*sqr(y)/nu;
    }
    else
    {
        dU =
          - 0.5*dpds
           *(
                yv/(kappa*sqrtk)*Foam::log(y/yv)
              + (y - yv)/(kappa*sqrtk)
              + sqr(yv)/nu
            );
    }

    return U + max(min(dU, U), -U);
}


//- Spalding's y+ for a given U+
inline scalar spaldingYPlus
(
    const scalar Uplus,
    const scalar kappa,
    const scalar E
)
{
    const scalar x = kappa*Uplus;
    return Uplus + (Foam::exp(x) - 1 - x - sqr(x)/2 - pow3(x)/6)/E;
}


//- dy+/dU+ of Spalding's law
inline scalar spaldingDYPlusDUPlus
(
    const scalar Uplus,
    const scalar kappa,
    const scalar E
)
{
    const scalar x = kappa*Uplus;
    return 1 + kappa*(Foam::exp(x) - 1 - x - sqr(x)/2)/E;
}


//- U+ of Spalding's law for a given y+ (Newton, from the smaller of the
//  linear and logarithmic estimates, which never overshoots)
inline scalar spaldingUPlus
(
    const scalar yPlus,
    const scalar kappa,
    const scalar E
)
{
    if (yPlus <= 0)
    {
        return 0;
    }

    scalar Uplus = yPlus;
    if (E*yPlus > 1)
    {
        Uplus = min(yPlus, Foam::log(E*yPlus)/kappa);
    }

    for (label i = 0; i < 100; i++)
    {
        const scalar f = spaldingYPlus(Uplus, kappa, E) - yPlus;
        const scalar dU = f/spaldingDYPlusDUPlus(Uplus, kappa, E);
        Uplus -= dU;

        if (mag(dU) < 1e-12*max(Uplus, 1))
        {
            break;
        }
    }

    return Uplus;
}


//- Eddy viscosity of Spalding's law, nu_t/nu = dy+/dU+ - 1 at y+: the
//  total stress (1 + nu_t/nu) dU+/dy+ = 1 of a constant-stress layer
//  that follows the law. ~ kappa y+ in the log layer, ~ y+^3 at the wall
inline scalar spaldingNutByNu
(
    const scalar yPlus,
    const scalar kappa,
    const scalar E
)
{
    return
        spaldingDYPlusDUPlus(spaldingUPlus(yPlus, kappa, E), kappa, E) - 1;
}


//- Friction velocity from Spalding's law for a velocity U at distance y
//  (bisection-safeguarded Newton; 0 for U <= 0)
inline scalar spaldingUTau
(
    const scalar U,
    const scalar y,
    const scalar nu,
    const scalar kappa,
    const scalar E
)
{
    if (U <= 0 || y <= 0)
    {
        return 0;
    }

    // f(uTau) = y uTau/nu - y+(U/uTau) increases monotonically with uTau
    auto f = [&](const scalar uTau)
    {
        return y*uTau/nu - spaldingYPlus(U/uTau, kappa, E);
    };

    // Laminar estimate: there y+ = U+ <= Spalding's y+(U+), so f <= 0
    scalar lo = Foam::sqrt(nu*U/y);
    scalar hi = 2*lo;
    for (label i = 0; i < 200 && f(hi) < 0; i++)
    {
        lo = hi;
        hi *= 2;
    }

    scalar uTau = 0.5*(lo + hi);
    for (label i = 0; i < 200; i++)
    {
        const scalar fu = f(uTau);
        if (fu < 0)
        {
            lo = uTau;
        }
        else
        {
            hi = uTau;
        }

        const scalar Uplus = U/uTau;
        const scalar df =
            y/nu + spaldingDYPlusDUPlus(Uplus, kappa, E)*Uplus/uTau;
        scalar next = uTau - fu/df;
        if (next <= lo || next >= hi)
        {
            next = 0.5*(lo + hi);
        }

        if (mag(next - uTau) < 1e-12*uTau)
        {
            return next;
        }
        uTau = next;
    }

    return uTau;
}


//- Exponential blending function of Kader (1981), as used by Popovac and
//  Hanjalic (2007): the viscous value is weighted by exp(-Gamma), the
//  logarithmic one by exp(-1/Gamma)
inline scalar blendingGamma(const scalar yStar)
{
    return 0.01*Foam::pow4(yStar)/(1 + 5*yStar);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace wallLaw
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
