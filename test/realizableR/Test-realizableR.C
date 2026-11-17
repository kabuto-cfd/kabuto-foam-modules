/*---------------------------------------------------------------------------*\
  Kabuto CFD modules for OpenFOAM-14
  Copyright (C) 2026 Kabuto CFD
-------------------------------------------------------------------------------
License
    This file is part of kabuto-foam-modules, add-on libraries for OpenFOAM,
    distributed under the GNU General Public License version 3 or later.

Application
    Test-realizableR

Description
    Checks fv::realizableR::clip: a realizable tensor is left alone; a
    non-realizable one gets its eigenvalues raised to the floor while its
    eigenvectors stay; repeated eigenvalues and a zero tensor work.
    Exits with status 1 on a failure.

\*---------------------------------------------------------------------------*/

#include "realizableR.H"
#include "tensor.H"
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

// Rotation about an oblique axis, so no test tensor is aligned with x, y, z
tensor rotation()
{
    const vector a = vector(1, 2, 3)/mag(vector(1, 2, 3));
    const scalar t = 0.7;
    const scalar c = Foam::cos(t), s = Foam::sin(t);
    return
        c*tensor::I + (1 - c)*sqr(a)
      + s*tensor(0, -a.z(), a.y(), a.z(), 0, -a.x(), -a.y(), a.x(), 0);
}

symmTensor rotated(const vector& lambda)
{
    const tensor Q = rotation();
    const tensor D(lambda.x(), 0, 0, 0, lambda.y(), 0, 0, 0, lambda.z());
    return symm(Q & D & Q.T());
}

} // End anonymous namespace


int main()
{
    const scalar minEig = 1e-6, relEig = 1e-3, tol = 1e-9;

    Info<< "realizable tensor" << nl;
    {
        const symmTensor r0 = rotated(vector(1, 2, 3));
        symmTensor r = r0;
        check(!fv::realizableR::clip(r, minEig, relEig), "not clipped");
        check(mag(r - r0) < tol, "unchanged");
    }

    Info<< "one negative eigenvalue" << nl;
    {
        const symmTensor r0 = rotated(vector(-0.5, 2, 3));
        symmTensor r = r0;
        check(fv::realizableR::clip(r, minEig, relEig), "clipped");

        const vector l0 = eigenValues(r0);
        const vector l = eigenValues(r);
        const scalar floor = max(minEig, relEig*tr(r0)/3);
        check(mag(cmptMin(l) - floor) < tol, "smallest eigenvalue = floor");
        check
        (
            mag(l.y() - l0.y()) < tol && mag(l.z() - l0.z()) < tol,
            "other eigenvalues kept"
        );

        const tensor e0 = eigenVectors(r0, l0);
        const vector v[3] = {e0.x(), e0.y(), e0.z()};
        const scalar lNew[3] = {floor, l0.y(), l0.z()};
        bool aligned = true;
        for (label i = 0; i < 3; i++)
        {
            aligned = aligned && mag((r & v[i]) - lNew[i]*v[i]) < 1e-8;
        }
        check(aligned, "eigenvectors kept");
    }

    Info<< "shear-dominated tensor (normal stresses small)" << nl;
    {
        symmTensor r(1e-4, 0.5, 0, 1e-4, 0, 1e-4);
        check(fv::realizableR::clip(r, minEig, relEig), "clipped");
        check(cmptMin(eigenValues(r)) > 0, "positive definite");
        check(r.xx() > 0 && r.yy() > 0 && r.zz() > 0, "normal stresses > 0");
    }

    Info<< "two equal negative eigenvalues" << nl;
    {
        symmTensor r = rotated(vector(-1, -1, 2));
        check(fv::realizableR::clip(r, minEig, relEig), "clipped");
        const vector l = eigenValues(r);
        check(cmptMin(l) > 0 && mag(l.z() - 2) < 1e-8, "floor, floor, 2");
    }

    Info<< "zero tensor" << nl;
    {
        symmTensor r = Zero;
        check(fv::realizableR::clip(r, minEig, relEig), "clipped");
        check(mag(r - minEig*symmTensor::I) < 1e-15, "minEigenvalue*I");
    }

    Info<< "very large tensor (diverging run: the cubic of R overflows)" << nl;
    {
        const scalar big = 1e150;
        const symmTensor r0 = rotated(vector(-0.5, 2, 3));
        symmTensor r1 = r0;
        fv::realizableR::clip(r1, minEig, relEig);

        symmTensor r = big*r0;
        check(fv::realizableR::clip(r, minEig, relEig), "clipped");
        bool finite = true;
        for (direction c = 0; c < symmTensor::nComponents; c++)
        {
            finite = finite && std::isfinite(r[c]);
        }
        check(finite, "finite");
        check(mag(r/big - r1) < tol, "= big*(clipped tensor)");
    }

    Info<< "10000 random tensors" << nl;
    {
        // Small linear congruential generator: the same numbers everywhere
        unsigned long long seed = 12345;
        auto rnd = [&seed]()
        {
            seed = seed*6364136223846793005ULL + 1442695040888963407ULL;
            return scalar(seed >> 11)/scalar(1ULL << 53)*2 - 1;
        };

        label nClipped = 0, nBad = 0;
        for (label i = 0; i < 10000; i++)
        {
            const symmTensor r0(rnd(), rnd(), rnd(), rnd(), rnd(), rnd());
            symmTensor r = r0;
            const bool clipped = fv::realizableR::clip(r, minEig, relEig);
            const scalar floor = max(minEig, relEig*max(tr(r0), 0)/3);
            const vector l = eigenValues(r);
            const vector l0 = eigenValues(r0);

            // Eigenvalues above the floor are kept, the others become it,
            // to the round-off of the analytical eigenvalues (~1e-9 |r|)
            const scalar eps = 1e-8*max(mag(r0), 1);
            bool ok = cmptMin(l) > floor - eps;
            if (clipped)
            {
                nClipped++;
                for (direction d = 0; d < 3; d++)
                {
                    const scalar expect = max(l0[d], floor);
                    ok = ok && mag(l[d] - expect) < eps;
                }
                // The clipped tensor must commute with the original one
                // (same eigenvectors)
                ok = ok && mag((r & r0) - (r0 & r)) < eps;
            }
            else
            {
                ok = ok && r == r0;
            }
            if (!ok && nBad < 3)
            {
                Info<< "    r0 " << r0 << " l0 " << l0 << " floor " << floor
                    << nl << "    r  " << r << " l " << l << " comm "
                    << mag((r & r0) - (r0 & r)) << nl;
            }
            if (!ok) nBad++;
        }
        Info<< "    " << nClipped << " clipped" << nl;
        check(nBad == 0, "eigenvalues max(lambda, floor), eigenvectors kept");
    }

    Info<< nl << (nFailed ? "FAILED" : "passed") << nl;
    return nFailed ? 1 : 0;
}


// ************************************************************************* //
