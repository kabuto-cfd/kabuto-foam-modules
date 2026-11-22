/*---------------------------------------------------------------------------*\
  Kabuto CFD modules for OpenFOAM-14
  Copyright (C) 2026 Kabuto CFD
-------------------------------------------------------------------------------
License
    This file is part of kabuto-foam-modules, add-on libraries for OpenFOAM,
    distributed under the GNU General Public License version 3 or later.

Application
    Test-wallLaw

Description
    Checks the point-wise laws of the wall of libkabutoWallFunctions
    (src/wallFunctions/wallLaw/wallLaw.H): Spalding's law and its inverses,
    the friction velocity from a velocity sample anywhere from y+ = 0.5 to
    3000, and Kim and Choudhury's pressure-gradient-sensitised velocity.
    Exits with status 1 on a failure.

\*---------------------------------------------------------------------------*/

#include "wallLaw.H"
#include "IOstreams.H"

using namespace Foam;

namespace
{

label nFailed = 0;

void check(const bool ok, const char* what)
{
    Info<< (ok ? "    ok    " : "    FAIL  ") << what << nl;
    if (!ok) nFailed++;
}

} // End anonymous namespace


int main()
{
    const scalar kappa = 0.41, E = 9.8, Cmu = 0.09, yvStar = 11.225;

    Info<< "Spalding's law" << nl;
    {
        scalar worst = 0;
        for (scalar Uplus = 0.1; Uplus < 35; Uplus *= 1.3)
        {
            const scalar yPlus = wallLaw::spaldingYPlus(Uplus, kappa, E);
            worst = max
            (
                worst,
                mag(wallLaw::spaldingUPlus(yPlus, kappa, E) - Uplus)/Uplus
            );
        }
        check(worst < 1e-10, "U+(y+(U+)) = U+");

        check
        (
            mag(wallLaw::spaldingUPlus(1e-3, kappa, E) - 1e-3) < 1e-9,
            "U+ = y+ in the viscous sublayer"
        );
        const scalar yp = 1000;
        check
        (
            mag(wallLaw::spaldingUPlus(yp, kappa, E) - Foam::log(E*yp)/kappa)
          < 0.01,
            "log law at y+ = 1000"
        );
        check
        (
            mag(wallLaw::spaldingDYPlusDUPlus(0, kappa, E) - 1) < 1e-15,
            "dU+/dy+ = 1 at the wall"
        );
    }

    Info<< "eddy viscosity of Spalding's law" << nl;
    {
        // Constant stress: (1 + nu_t/nu) dU+/dy+ = 1, by finite differences
        scalar worst = 0;
        for (scalar yp = 0.5; yp < 3000; yp *= 1.5)
        {
            const scalar h = 1e-4*yp;
            const scalar dUdy =
            (
                wallLaw::spaldingUPlus(yp + h, kappa, E)
              - wallLaw::spaldingUPlus(yp - h, kappa, E)
            )/(2*h);
            worst = max
            (
                worst,
                mag((1 + wallLaw::spaldingNutByNu(yp, kappa, E))*dUdy - 1)
            );
        }
        check(worst < 1e-6, "(1 + nu_t/nu) dU+/dy+ = 1");

        const scalar yp = 3000;
        check
        (
            mag(wallLaw::spaldingNutByNu(yp, kappa, E)/(kappa*yp) - 1) < 0.01,
            "nu_t/nu -> kappa y+ in the log layer"
        );

        const scalar order = Foam::log
        (
            wallLaw::spaldingNutByNu(0.2, kappa, E)
           /wallLaw::spaldingNutByNu(0.1, kappa, E)
        )/Foam::log(scalar(2));
        check(order > 2.9 && order < 3.1, "nu_t/nu ~ y+^3 at the wall");
    }

    Info<< "friction velocity from a sample at y+ = 0.5 .. 3000" << nl;
    {
        const scalar nu = 1.5e-5, uTau = 0.3;
        scalar worst = 0;
        for (scalar yPlus = 0.5; yPlus < 3000; yPlus *= 1.5)
        {
            const scalar y = yPlus*nu/uTau;
            const scalar U = uTau*wallLaw::spaldingUPlus(yPlus, kappa, E);
            worst = max
            (
                worst,
                mag(wallLaw::spaldingUTau(U, y, nu, kappa, E) - uTau)/uTau
            );
        }
        check(worst < 1e-9, "uTau recovered");
        check
        (
            wallLaw::spaldingUTau(-1, 1e-3, nu, kappa, E) == 0,
            "uTau = 0 for U <= 0"
        );
    }

    Info<< "pressure-gradient-sensitised velocity" << nl;
    {
        const scalar nu = 1e-5, k = 0.01, U = 2, y = 5e-3;
        check
        (
            wallLaw::UTilde(U, 0, y, k, nu, Cmu, kappa, yvStar) == U,
            "UTilde = U without a pressure gradient"
        );

        // Laminar profile with a pressure gradient, U = tau y/nu + dpds y^2/2nu,
        // inside the viscous sublayer: UTilde gives tau y/nu exactly
        const scalar uk = Foam::pow(Cmu, 0.25)*Foam::sqrt(k);
        const scalar yv = yvStar*nu/uk;
        const scalar tau = 1e-3, dpds = 0.5, yl = 0.5*yv;
        const scalar Ul = tau*yl/nu + 0.5*dpds*sqr(yl)/nu;
        check
        (
            mag
            (
                wallLaw::UTilde(Ul, dpds, yl, k, nu, Cmu, kappa, yvStar)
              - tau*yl/nu
            ) < 1e-12*Ul,
            "laminar profile: tau_w = nu UTilde/y"
        );

        const scalar below =
            wallLaw::UTilde(U, dpds, yv*(1 - 1e-9), k, nu, Cmu, kappa, yvStar);
        const scalar above =
            wallLaw::UTilde(U, dpds, yv*(1 + 1e-9), k, nu, Cmu, kappa, yvStar);
        check(mag(above - below) < 1e-8*U, "continuous at y = y_v");

        check
        (
            wallLaw::UTilde(U, dpds, y, k, nu, Cmu, kappa, yvStar) < U
         && wallLaw::UTilde(U, -dpds, y, k, nu, Cmu, kappa, yvStar) > U,
            "adverse gradient lowers it, favourable raises it"
        );
        check
        (
            wallLaw::UTilde(1e-3, 1e3, y, k, nu, Cmu, kappa, yvStar) >= 0
         && wallLaw::UTilde(1e-3, -1e3, y, k, nu, Cmu, kappa, yvStar)
         <= 2e-3*(1 + 1e-12),
            "limited to 0 <= UTilde <= 2 U"
        );
    }

    Info<< "blending" << nl;
    {
        check
        (
            Foam::exp(-1/max(wallLaw::blendingGamma(1), small)) < 1e-100
         && Foam::exp(-wallLaw::blendingGamma(100)) < 1e-100,
            "viscous at y* = 1, logarithmic at y* = 100"
        );
    }

    Info<< nl << (nFailed ? "FAILED" : "passed") << nl;
    return nFailed ? 1 : 0;
}


// ************************************************************************* //
