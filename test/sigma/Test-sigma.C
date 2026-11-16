/*---------------------------------------------------------------------------*\
  Kabuto CFD modules for OpenFOAM-14
  Copyright (C) 2026 Kabuto CFD
-------------------------------------------------------------------------------
License
    This file is part of kabuto-foam-modules, add-on libraries for OpenFOAM,
    distributed under the GNU General Public License version 3 or later.

Application
    Test-sigma

Description
    Checks the differential operator of the sigma subgrid-scale model
    (src/momentumTransportModels/LES/sigma/sigmaOperator.H): the singular
    values against the eigenvalues of g^T g, a known value, invariance under
    rotation, and that D_sigma vanishes for pure shear, solid-body rotation,
    isotropic and axisymmetric expansion and two-dimensional gradients, and
    decays as y^3 at a no-slip wall. Exits with status 1 on a failure.

\*---------------------------------------------------------------------------*/

#include "sigmaOperator.H"
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

// Deterministic pseudo-random numbers in [-1, 1]
scalar rnd()
{
    static unsigned long s = 12345;
    s = (1103515245*s + 12345) % 2147483648UL;
    return 2*scalar(s)/2147483648.0 - 1;
}

tensor rndTensor()
{
    return tensor
    (
        rnd(), rnd(), rnd(), rnd(), rnd(), rnd(), rnd(), rnd(), rnd()
    );
}

// Rotation by angle t about the unit axis a
tensor rotation(const vector& a, const scalar t)
{
    const scalar c = Foam::cos(t), s = Foam::sin(t);
    return
        c*tensor::I + (1 - c)*(a*a)
      + s*tensor(0, -a.z(), a.y(), a.z(), 0, -a.x(), -a.y(), a.x(), 0);
}

} // End anonymous namespace


int main()
{
    Info<< "Singular values" << nl;
    {
        scalar worst = 0;
        for (label i = 0; i < 2000; i++)
        {
            const tensor g(rndTensor());
            const vector s(sigmaOperator::singularValues(g));
            const vector l(eigenValues(symm(g.T() & g)));   // ascending
            const scalar scale = max(l.z(), small);
            worst = max(worst, mag(sqr(s.x()) - l.z())/scale);
            worst = max(worst, mag(sqr(s.y()) - l.y())/scale);
            worst = max(worst, mag(sqr(s.z()) - l.x())/scale);
        }
        check(worst < 1e-8, "sigma_i^2 = eigenvalues of g^T g (2000 random)");

        const vector s(sigmaOperator::singularValues(tensor(3, 0, 0, 0, 2, 0, 0, 0, 1)));
        check
        (
            mag(s - vector(3, 2, 1)) < 1e-12,
            "diag(3, 2, 1) -> (3, 2, 1)"
        );
        check
        (
            mag(sigmaOperator::Dsigma(tensor(3, 0, 0, 0, 2, 0, 0, 0, 1)) - 1.0/9)
          < 1e-12,
            "D_sigma(diag(3, 2, 1)) = 1/9"
        );
    }

    Info<< "Invariance and sign" << nl;
    {
        scalar worst = 0, minD = great;
        for (label i = 0; i < 500; i++)
        {
            const tensor g(rndTensor());
            const vector a(vector(rnd(), rnd(), rnd() + 2)/mag(vector(rnd(), rnd(), rnd() + 2)));
            const tensor Q(rotation(a/mag(a), 3*rnd()));
            const scalar D = sigmaOperator::Dsigma(g);
            const scalar DQ = sigmaOperator::Dsigma(Q & g & Q.T());
            worst = max(worst, mag(D - DQ)/max(mag(g), small));
            minD = min(minD, D);
        }
        check(worst < 1e-10, "D_sigma(Q g Q^T) = D_sigma(g)");
        check(minD >= 0, "D_sigma >= 0");
    }

    Info<< "Flows without subgrid dissipation" << nl;
    {
        const scalar gamma = 7.3;
        const scalar tol = 1e-10*gamma;

        // u = gamma y
        check
        (
            sigmaOperator::Dsigma(tensor(0, 0, 0, gamma, 0, 0, 0, 0, 0)) < tol,
            "pure shear"
        );

        // u = Omega x r about z
        check
        (
            sigmaOperator::Dsigma(tensor(0, gamma, 0, -gamma, 0, 0, 0, 0, 0))
          < tol,
            "solid-body rotation"
        );

        // Solid-body rotation about an oblique axis
        {
            const vector w(vector(1, 2, 3)/mag(vector(1, 2, 3)));
            const tensor W
            (
                gamma*tensor(0, -w.z(), w.y(), w.z(), 0, -w.x(), -w.y(), w.x(), 0)
            );
            check(sigmaOperator::Dsigma(W) < tol, "rotation about (1, 2, 3)");
        }

        check
        (
            sigmaOperator::Dsigma(gamma*tensor::I) < tol,
            "isotropic expansion"
        );
        check
        (
            sigmaOperator::Dsigma(tensor(-2*gamma, 0, 0, 0, gamma, 0, 0, 0, gamma))
          < tol,
            "axisymmetric contraction"
        );

        // Any two-dimensional gradient (u = u(x, y), w = 0)
        scalar worst = 0;
        for (label i = 0; i < 200; i++)
        {
            const tensor g(rnd(), rnd(), 0, rnd(), rnd(), 0, 0, 0, 0);
            worst = max(worst, sigmaOperator::Dsigma(g)/max(mag(g), small));
        }
        check(worst < 1e-10, "two-dimensional gradients (200 random)");

        // Contrast: WALE's operator |Sd|^3/(|S|^5 + |Sd|^5/2) is non-zero in
        // solid-body rotation; checked here so the test documents the
        // difference
        const tensor W(0, gamma, 0, -gamma, 0, 0, 0, 0, 0);
        const symmTensor Sd(dev(symm(W & W)));
        check(magSqr(Sd) > small, "(WALE's operator is not zero in rotation)");
    }

    Info<< "Near-wall behaviour" << nl;
    {
        // Taylor expansion at a no-slip wall y = 0 (x = z = 0):
        //   u = y (a1 + b1 x + c1 z) + d1 y^2,  w = y (a3 + b3 x + c3 z) + d3 y^2,
        //   v = e y^2 with e = -(b1 + c3)/2 (continuity);
        // g_ij = du_i/dx_j. D_sigma ~ y^3 (Nicoud et al. 2011)
        const scalar a1 = 1.0, b1 = 0.3, c1 = -0.5, d1 = 0.7;
        const scalar a3 = 0.6, b3 = -0.2, c3 = 0.4, d3 = -0.9;
        const scalar e = -(b1 + c3)/2;
        auto gWall = [&](const scalar y)
        {
            return tensor
            (
                b1*y, a1 + 2*d1*y, c1*y,
                0,    2*e*y,       0,
                b3*y, a3 + 2*d3*y, c3*y
            );
        };
        const scalar D1 = sigmaOperator::Dsigma(gWall(1e-2));
        const scalar D2 = sigmaOperator::Dsigma(gWall(0.5e-2));
        const scalar order = Foam::log(D1/D2)/Foam::log(scalar(2));
        Info<< "    order of D_sigma at the wall: " << order << nl;
        check(order > 2.8 && order < 3.2, "D_sigma ~ y^3");
    }

    Info<< nl << (nFailed ? "FAILED" : "passed") << nl;

    return nFailed ? 1 : 0;
}


// ************************************************************************* //
