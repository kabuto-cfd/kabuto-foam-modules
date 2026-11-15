# kabuto-foam-modules

Add-on modules for [OpenFOAM-14](https://openfoam.org) (OpenFOAM Foundation release) by Kabuto CFD.
Every module builds with `wmake` against an unmodified OpenFOAM-14 and is loaded at run time
through `libs (...)`; nothing in OpenFOAM itself is changed.

## Modules

| Module | Library | What it does |
|---|---|---|
| [`realizableR`](src/fvConstraints/realizableR/realizableR.H) | `libkabutoFvConstraints.so` | fvConstraint that keeps a Reynolds-stress field realizable: in cells where the smallest eigenvalue of **R** drops below a floor, the eigenvalues are clipped to that floor and **R** is rebuilt from its unchanged eigenvectors |
| [`stressOmega`](src/momentumTransportModels/RAS/stressOmega/stressOmega.H) | `libkabutoMomentumTransportModels.so` | Wilcox's stress-ω Reynolds-stress model (2006): six stress equations with the LRR-QI pressure-strain model, closed with the ω equation of k-ω (2006); integrates to the wall or runs with wall functions; incompressible and compressible |

Planned: non-equilibrium wall functions, and solver modules for `foamRun`.

### realizableR

OpenFOAM's `ReynoldsStress::boundNormalStress` bounds the normal stresses from below but
leaves the shear stresses free. On meshes with thin prism layers this can let **R** lose
positive definiteness (shear components grow while the normal ones sit at `kMin`) until an
LRR or SSG run diverges. The limiter acts right after **R** is solved for (the models call
`fvConstraints.constrain(R)` there):

```
floor = max(minEigenvalue, relativeMinEigenvalue * tr(R)/3)
R     = sum_i max(lambda_i, floor) e_i e_i^T      (only where min(lambda) < floor)
```

In an unsteady LRR run of a cyclone with layered walls it touched about 0.1 % of the cells per
time step and let the run go on past the point where it had diverged without it.

`system/fvConstraints`:

```
realizableR
{
    type        realizableR;
    libs        ("libkabutoFvConstraints.so");

    field       R;                    // optional, default R
    minEigenvalue          1e-6;      // optional [m^2/s^2], default 1e-6
    relativeMinEigenvalue  1e-3;      // optional, default 1e-3
}
```

The log reports, every time it acts, how many cells were clipped:
`realizableR realizableR: <n> of <N> cells clipped`.

### stressOmega

The ε-based LRR and SSG of OpenFOAM need wall functions or damping near walls. Wilcox's
stress-ω model carries the length scale in ω instead, so it can be integrated down to the wall
(R = 0 there, ω from `omegaWallFunction`) like k-ω, while keeping the full stress anisotropy
that drives secondary flows in ducts and the effects of swirl and curvature.

```
D(R)/Dt     = P + Π - 2/3 β* ω k I + ∇·((ν + σ* k/ω) ∇R)
Π           = β* C1 ω (2/3 k I - R) - α̂ dev(P) - β̂ dev(D) - γ̂ k dev(S)
D(ω)/Dt     = γ G ω/k - β ω² + σd/ω ∇k·∇ω + ∇·((ν + σ k/ω) ∇ω)
```

with P = −twoSymm(R·∇U), D = −twoSymm(R·∇Uᵀ), G = tr(P)/2, β = β0 fβ (vortex stretching),
and the 2006 coefficients (C1 = 1.8, C2 = 10/19, γ = 0.52, β0 = 0.0708, β* = 0.09, σ = 0.5,
σ* = 0.6, σdo = 1/8), all settable in `stressOmegaCoeffs`.

`constant/momentumTransport` and `system/controlDict`:

```
RAS { model stressOmega; turbulence on; }

libs ("libkabutoMomentumTransportModels.so");
```

Fields: `0/R`, `0/omega` (k is computed from R), `0/nut`. Divergence schemes for `div(phi,R)`
and `div(phi,omega)`.

Checks (`tutorials/stressOmega/`):

- `channel395`: plane channel integrated to the wall (first cell at y+ ≈ 0.8) at
  Re_b = 13 750. Re_τ = 400.5 against 399 from Dean's correlation (k-ω 2006 on the same mesh:
  388.6); U⁺ follows y⁺ below y⁺ ≈ 5 and the log law (κ = 0.41, B = 5.0) within 0.4 at
  y⁺ = 30–200; the normal stresses come out anisotropic in the right order
  (⟨u'u'⟩ > ⟨w'w'⟩ > ⟨v'v'⟩).
- `ductSecondaryFlow`: OpenFOAM's duct tutorial (wall functions) with the model switched from
  LRR. Run to convergence, both give the corner vortices that eddy-viscosity models miss, with
  the flow along the corner bisector going into the corner: largest secondary velocity 2.95 %
  of the bulk velocity (LRR 2.59 %), 1.2 % into the corner at 5 mm from it (LRR 1.7 %).
  stressOmega needed 14 931 SIMPLE iterations to the tutorial's residual targets, LRR 5 207.

The compressible variant is built and registered but not yet run on a test case.

## Build

```bash
source /opt/openfoam14/etc/bashrc     # or <your OpenFOAM-14>/etc/bashrc
git clone https://github.com/kabuto-cfd/kabuto-foam-modules.git
cd kabuto-foam-modules
./Allwmake -j
```

The libraries go to `$FOAM_USER_LIBBIN`, which OpenFOAM searches when a case names them in
`libs`. `./Allwclean` removes the build.

## Tests and tutorials

- `Test-realizableR` (built by `Allwmake`) checks the eigenvalue clipping on a few tensors.
- `tutorials/realizableR/ductSecondaryFlow/Allrun [iterations]` copies OpenFOAM-14's
  `ductSecondaryFlow` tutorial (LRR in a square duct), adds the limiter and runs it.
- `tutorials/stressOmega/ductSecondaryFlow/Allrun [iterations]` runs the same duct with
  stressOmega.
- `tutorials/stressOmega/channel395/Allrun` runs the wall-resolved channel and fails when Re_τ
  is more than 3 % off Dean's correlation.

GitHub Actions builds every module against the openfoam.org OpenFOAM-14 package on
Ubuntu 24.04 and runs all of these for every change.

## Licence

GNU General Public License v3 or later (`COPYING`), as OpenFOAM itself.
Copyright (C) 2026 Kabuto CFD.

## Trademarks

This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the
OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM® and OpenCFD® trade marks,
nor by the OpenFOAM Foundation.
