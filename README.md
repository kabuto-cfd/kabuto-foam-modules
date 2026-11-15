# kabuto-foam-modules

Add-on modules for [OpenFOAM-14](https://openfoam.org) (OpenFOAM Foundation release) by Kabuto CFD.
Every module builds with `wmake` against an unmodified OpenFOAM-14 and is loaded at run time
through `libs (...)`; nothing in OpenFOAM itself is changed.

## Modules

| Module | Library | What it does |
|---|---|---|
| [`realizableR`](src/fvConstraints/realizableR/realizableR.H) | `libkabutoFvConstraints.so` | fvConstraint that keeps a Reynolds-stress field realizable: in cells where the smallest eigenvalue of **R** drops below a floor, the eigenvalues are clipped to that floor and **R** is rebuilt from its unchanged eigenvectors |
| [`stressOmega`](src/momentumTransportModels/RAS/stressOmega/stressOmega.H) | `libkabutoMomentumTransportModels.so` | Wilcox's stress-ω Reynolds-stress model (2006): six stress equations with the LRR-QI pressure-strain model, closed with the ω equation of k-ω (2006); integrates to the wall or runs with wall functions; incompressible and compressible |
| [`nut/epsilon/omegaNonEquilibriumWallFunction`](src/wallFunctions/nutNonEquilibriumWallFunction/nutNonEquilibriumWallFunctionFvPatchScalarField.H) | `libkabutoWallFunctions.so` | Kim & Choudhury's non-equilibrium wall functions: log law sensitised to the pressure gradient along the wall, production and dissipation in the wall cell averaged over a two-layer profile; first cell in the log layer |
| [`nut/epsilon/omegaNonEquilibriumSpaldingWallFunction`](src/wallFunctions/nutNonEquilibriumSpaldingWallFunction/nutNonEquilibriumSpaldingWallFunctionFvPatchScalarField.H) | `libkabutoWallFunctions.so` | y⁺-insensitive version: Spalding's law with the same pressure-gradient correction, production from the Spalding profile under the local total stress; first cell anywhere from y⁺ ≈ 1 to the log layer |
| [`twoLayerEpsilon`](src/wallFunctions/twoLayerEpsilon/twoLayerEpsilon.H) | `libkabutoWallFunctions.so` | fvConstraint: Wolfstein's two-layer ε and ν_t in the viscosity-affected cells (Re_y < 200), which ε models need to be y⁺-insensitive |

Planned: solver modules for `foamRun`.

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

### Non-equilibrium wall functions

Two families, each a nut boundary condition and the matching epsilon or omega one (Cmu, κ, E
and `yvStar` are set on the nut condition; `p` / `rho` name the pressure and density fields,
p in Pa is divided by ρ):

| | `…NonEquilibriumWallFunction` | `…NonEquilibriumSpaldingWallFunction` |
|---|---|---|
| velocity law | log law (y* > 11.5), laminar below | Spalding's law, any y⁺ |
| pressure gradient | Ũ = U − ½ (dp/ds) […] (Kim & Choudhury) | the same Ũ |
| G and ε / ω in the wall cell | averages over the cell of the two-layer profile | G = u_τ⁴/ν (τ⁺ − γ) γ with γ = dU⁺/dy⁺ of Spalding's profile and τ⁺ = 1 + p⁺y⁺; ε blended viscous/log (Kader), ω = (ω_vis² + ω_log²)^½ |
| first cell | y⁺ ≳ 30 | y⁺ ≈ 1 … 100+ (ε models with `twoLayerEpsilon`) |

The pressure-gradient correction is limited to 0 ≤ Ũ ≤ 2U, so it fades where the near-wall
velocity vanishes (separation and reattachment lines); unlimited, it fed back on itself there.

```
wall { type nutNonEquilibriumSpaldingWallFunction;  value uniform 0; }        // 0/nut
wall { type omegaNonEquilibriumSpaldingWallFunction; value $internalField; }  // 0/omega
wall { type kqRWallFunction;                          value $internalField; }  // 0/k

libs ("libkabutoWallFunctions.so");                                            // controlDict
```

For an ε model, the y⁺-insensitive family also needs `twoLayerEpsilon` in `system/fvConstraints`
(and `wallDist { method meshWave; }` in `fvSchemes`): a high-Re ε model in the cells between the
wall cell and y⁺ ≈ 30 overpredicts the mixing by itself, which no wall function can correct.
With it the ε equation's initial residual stays near 1e-2 (the near-wall ε is imposed after each
solve); judge convergence from the other fields.

Checks:

- `tutorials/wallFunctions/channelYPlus`: plane channel at Re_τ = 1864 (Dean's correlation),
  first cell centre at y⁺ = 1 … 100. Deviation of Re_τ from 1864 (`./Allrun -all`):

  | model | treatment | y⁺ = 1 | 2 | 5 | 11 | 20 | 40 | 100 |
  |---|---|---|---|---|---|---|---|---|
  | k-ε | OpenFOAM standard | +37.4 % | +36.3 % | +28.7 % | +10.2 % | −3.0 % | −2.9 % | −2.9 % |
  | k-ε | NonEquilibrium | +37.4 % | +36.3 % | +28.7 % | +10.2 % | −1.2 % | −4.5 % | −6.8 % |
  | k-ε | NonEquilibriumSpalding | +37.3 % | +35.9 % | +30.0 % | +7.9 % | +1.5 % | −1.0 % | −1.7 % |
  | k-ε | NonEquilibriumSpalding + twoLayerEpsilon | −2.0 % | −0.9 % | +4.6 % | +4.1 % | +1.9 % | −1.4 % | −1.8 % |
  | k-ω SST | OpenFOAM standard | −3.4 % | −4.6 % | +0.8 % | +11.6 % | −3.1 % | −3.0 % | −2.9 % |
  | k-ω SST | NonEquilibrium | −3.4 % | −4.6 % | +0.8 % | +0.9 % | −2.3 % | −5.0 % | −6.8 % |
  | k-ω SST | NonEquilibriumSpalding | −3.5 % | −4.7 % | +0.1 % | +3.2 % | +2.0 % | −0.1 % | −1.8 % |

  The channel is driven by a body force (meanVelocityForce), so the pressure gradient is zero
  and the table shows the near-wall treatment alone. The two-layer averaging of the
  NonEquilibrium family lowers the friction by up to 7 % on coarse wall cells (y⁺ = 100).
- `tutorials/wallFunctions/pitzDailySteady`: OpenFOAM's backward-facing step (k-ε; first-cell
  y⁺ 0.3 … 27 on the lower wall). Reattachment length: standard 6.76 h, NonEquilibrium 6.72 h,
  NonEquilibriumSpalding 6.95 h, with twoLayerEpsilon 6.91 h. On this mesh the near-wall cells
  sit mostly in the buffer layer, where the NonEquilibrium family falls back to its laminar
  branch, so this shows that the treatments are stable in separated flow rather than which is
  more accurate.
- `Test-wallLaw`: Spalding's law and its inverses, the friction velocity from a sample at
  y⁺ = 0.5 … 3000, and the pressure-gradient correction (no gradient, exact laminar profile,
  continuity at y_v, the limit).

References: Kim & Choudhury (1995), ASME FED 217; the form of the Fluent theory guide.
Spalding (1961), J. Appl. Mech. 28, 455. Kader (1981), Int. J. Heat Mass Transfer 24, 1541.
Popovac & Hanjalić (2007), Flow Turbul. Combust. 78, 177. Wolfshtein (1969), Int. J. Heat
Mass Transfer 12, 301.

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
- `Test-wallLaw` checks the laws of the wall used by the wall functions.
- `tutorials/wallFunctions/channelYPlus/Allrun [-all]` runs the channel y⁺ sweep and fails when
  a y⁺-insensitive treatment is more than 6 % off; `tutorials/wallFunctions/pitzDailySteady/Allrun
  [std|ne|nes|nes2l]` runs the backward-facing step and fails when it does not reattach between
  5 and 9 step heights.

GitHub Actions builds every module against the openfoam.org OpenFOAM-14 package on
Ubuntu 24.04 and runs all of these for every change.

## Licence

GNU General Public License v3 or later (`COPYING`), as OpenFOAM itself.
Copyright (C) 2026 Kabuto CFD.

## Trademarks

This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the
OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM® and OpenCFD® trade marks,
nor by the OpenFOAM Foundation.
