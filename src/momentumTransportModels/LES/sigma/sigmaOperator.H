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
    Foam::sigmaOperator

Description
    Point-wise differential operator of the sigma subgrid-scale model
    (Nicoud et al. 2011), unit-tested by test/sigma.

    The singular values sigma1 >= sigma2 >= sigma3 >= 0 of the velocity
    gradient g are the square roots of the eigenvalues of G = g^T g, taken in
    closed form from the invariants of G (Nicoud et al. 2011, appendix):
    \verbatim
        I1 = tr(G),  I2 = (I1^2 - tr(G G))/2,  I3 = det(G)
        a1 = I1^2/9 - I2/3
        a2 = I1^3/27 - I1 I2/6 + I3/2
        a3 = acos(a2/a1^(3/2))/3
        sigma1 = sqrt(I1/3 + 2 sqrt(a1) cos(a3))
        sigma2 = sqrt(I1/3 - 2 sqrt(a1) cos(pi/3 + a3))
        sigma3 = |det(g)|/(sigma1 sigma2)

        D_sigma = sigma3 (sigma1 - sigma2)(sigma2 - sigma3)/sigma1^2
    \endverbatim
    The closed form never fails (unlike a general cubic solver on a matrix
    with nearly repeated eigenvalues); the argument of acos is clipped to
    [-1, 1] and the eigenvalues to >= 0 against round-off. sigma3 comes from
    the determinant rather than its own root, which would carry the square
    root of the round-off (~1e-8 sigma1) and spoil the zero of D_sigma for
    two-dimensional gradients and its y^3 decay at a wall.

    D_sigma vanishes whenever the resolved velocity gradient is
    two-dimensional or two-component (sigma3 = 0: pure shear, solid-body
    rotation), for isotropic expansion or contraction (sigma1 = sigma2 =
    sigma3) and for axisymmetric expansion or contraction (two equal singular
    values), and it decays as y^3 at a no-slip wall.

SourceFiles
    sigmaOperator.H (inline only)

\*---------------------------------------------------------------------------*/

#ifndef sigmaOperator_H
#define sigmaOperator_H

#include "tensor.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace sigmaOperator
{

//- Singular values (sigma1, sigma2, sigma3), descending, of a velocity
//  gradient (either orientation: g and g^T have the same singular values)
inline vector singularValues(const tensor& g)
{
    const tensor G(g.T() & g);

    const scalar I1 = tr(G);
    const scalar I2 = 0.5*(sqr(I1) - tr(G & G));
    const scalar I3 = det(G);

    const scalar a1 = sqr(I1)/9 - I2/3;

    if (a1 <= small*sqr(max(I1, vSmall)))
    {
        // Three (nearly) equal eigenvalues
        const scalar s = Foam::sqrt(max(I1/3, 0));
        return vector(s, s, s);
    }

    const scalar a2 = pow3(I1)/27 - I1*I2/6 + I3/2;
    const scalar a3 =
        Foam::acos(max(min(a2/Foam::pow(a1, 1.5), scalar(1)), scalar(-1)))
       /3;

    const scalar piBy3 = constant::mathematical::pi/3;
    const scalar r = 2*Foam::sqrt(a1);

    scalar l1 = I1/3 + r*Foam::cos(a3);
    scalar l2 = I1/3 - r*Foam::cos(piBy3 + a3);
    scalar l3 = I1/3 - r*Foam::cos(piBy3 - a3);

    // a3 in [0, pi/3] orders them already; enforce it against round-off
    if (l2 > l1) Swap(l1, l2);
    if (l3 > l2) Swap(l2, l3);
    if (l2 > l1) Swap(l1, l2);

    const scalar s1 = Foam::sqrt(max(l1, 0));
    const scalar s2 = Foam::sqrt(max(l2, 0));

    // The smallest from sigma1 sigma2 sigma3 = |det(g)|: exactly zero for a
    // two-dimensional gradient, and free of the square root of the round-off
    // of l3 (~sqrt(eps) sigma1) that the closed form leaves
    const scalar s3 =
        s1*s2 > vSmall
      ? min(mag(det(g))/(s1*s2), s2)
      : Foam::sqrt(max(l3, 0));

    return vector(s1, s2, s3);
}


//- The sigma model's differential operator D_sigma (1/s)
inline scalar Dsigma(const tensor& g)
{
    const vector s(singularValues(g));

    if (s.x() <= vSmall)
    {
        return 0;
    }

    return max
    (
        s.z()*(s.x() - s.y())*(s.y() - s.z())/sqr(s.x()),
        scalar(0)
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace sigmaOperator
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
