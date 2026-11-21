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
| [`sigma`](src/momentumTransportModels/LES/sigma/sigma.H) | `libkabutoMomentumTransportModels.so` | Nicoud et al.'s σ subgrid-scale model (2011): eddy viscosity from the singular values of the velocity gradient; zero in pure shear and, unlike WALE, in solid-body rotation (vortex cores, swirl); y³ at walls without damping; incompressible and compressible |
| [`nutSpaldingLESWallFunction`](src/wallFunctions/nutSpaldingLESWallFunction/nutSpaldingLESWallFunctionFvPatchScalarField.H) | `libkabutoWallFunctions.so` | Wall-stress model for wall-modelled LES: Spalding's law fed with the velocity of a cell layer away from the wall and / or filtered in time, the remedies for the log-layer mismatch |

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
  Re_b = 13 750. Re_τ = 385.4 against 392 from the DNS of Moser, Kim & Mansour (1999) (k-ω 2006 on
  the same mesh: 388.6); U⁺ follows y⁺ below y⁺ ≈ 5 and the log law (κ = 0.41, B = 5.0) within 1.2 at
  y⁺ = 30–200; the normal stresses come out anisotropic in the right order
  (⟨u'u'⟩ > ⟨w'w'⟩ > ⟨v'v'⟩).
- `ductSecondaryFlow`: OpenFOAM's duct tutorial (wall functions) with the model switched from
  LRR. Run to convergence, both give the corner vortices that eddy-viscosity models miss, with
  the flow along the corner bisector going into the corner: largest secondary velocity 2.95 %
  of the bulk velocity (LRR 2.59 %), 1.2 % into the corner at 5 mm from it (LRR 1.7 %).
  stressOmega needed 14 931 SIMPLE iterations to the tutorial's residual targets, LRR 5 207.

The compressible variant is built and registered but not yet run on a test case.

### sigma

WALE's eddy viscosity vanishes in pure shear and decays as y³ at a wall, but it does not vanish
in solid-body rotation: in a forced vortex (the core of a cyclone or a swirl tube, which
laminarises) it returns a subgrid viscosity proportional to the rotation rate. Smagorinsky's
model has the opposite defect. The σ model of Nicoud et al. (2011) is built from the singular
values σ1 ≥ σ2 ≥ σ3 of the resolved velocity gradient g and vanishes in both:

```
nut     = (Cσ Δ)² Dσ,   Dσ = σ3 (σ1 − σ2)(σ2 − σ3)/σ1²
k       = (nut/(Ck Δ))²,   ε = Ce k^(3/2)/Δ
```

σ1, σ2 come in closed form from the invariants of gᵀg; σ3 = |det g|/(σ1 σ2), which keeps the
zero of Dσ for two-dimensional gradients and its y³ decay exact (its own square root would carry
√ε round-off). Cσ = 1.35 is the value of the reference; implementations that write
nut = C Δ² Dσ with C = 1.68 correspond to Cσ = 1.30.

```
LES { model sigma; delta cubeRootVol; cubeRootVol { deltaCoeff 1; } }   // momentumTransport

libs ("libkabutoMomentumTransportModels.so");                            // controlDict
```

Checks: `Test-sigma` (singular values against the eigenvalues of gᵀg for 2000 random tensors,
Dσ(diag(3, 2, 1)) = 1/9, invariance under rotation, Dσ = 0 for pure shear, rotation about any
axis, isotropic and axisymmetric expansion and two-dimensional gradients, order 2.99 of the
decay at a no-slip wall from a Taylor expansion that satisfies continuity);
`tutorials/sigma/channel395/Allrun [steps]` runs OpenFOAM-14's wall-resolved channel with the
model switched from WALE.

Reference: Nicoud, Baya Toda, Cabrit, Bose & Lee (2011), Phys. Fluids 23, 085106.

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

### nutSpaldingLESWallFunction

In a wall-modelled LES the first cell sits in the log layer and a wall function supplies the wall
stress. Fed with the instantaneous velocity of the wall cell itself (`nutUSpaldingWallFunction`),
the modelled stress follows the same cell's fluctuations, the resolved turbulence next to the
wall is damped, and the mean velocity of the cells above lies above the log law: the log-layer
mismatch (Kawai & Larsson 2012; Yang, Park & Moin 2017). The wall stress is then too low for the
flow outside. This boundary condition takes the input velocity of Spalding's law elsewhere:

- `sampleLayer n`: from the cell n layers from the wall (1 = the wall cell), found by walking
  from the wall cell through the face whose normal is closest to the wall normal (the cell
  straight above in prism layers); the walk stops early at a processor or domain boundary and
  the start of the run reports on how many faces it did;
- `Tf`: through an exponential time filter, U_f ← (1 − ε) U_f + ε U_s with ε = Δt/T_f, once per
  time step; U_f is written with the field so that a restart continues it.

```
uτ: y_s uτ/ν = y⁺(|U_in|/uτ) (Spalding),   ν_t,w = max(uτ²/|∂U/∂n|_w − ν, 0)
```

with y_s the distance of the sample cell centre from the wall and U_in its wall-parallel velocity
relative to the wall (filtered if `Tf` > 0). With `sampleLayer 1` and `Tf 0` it is
`nutUSpaldingWallFunction`. No coefficient is added; κ = 0.41, E = 9.8 as upstream.

```
wall { type nutSpaldingLESWallFunction; sampleLayer 2; Tf 4; value uniform 0; }   // 0/nut

libs ("libkabutoWallFunctions.so");                                                // controlDict
```

Check: plane channel at Re_τ ≈ 2000 (Dean) built from OpenFOAM-14's `channel395` (4 × 2 × 2) on a
coarse uniform grid of 40 × 30 × 30 cells (wall-cell centre at y⁺ ≈ 67, Δx⁺ ≈ 200, Δz⁺ ≈ 133),
16 000 steps (≈ 20 h/u_τ), Re_τ from the driving pressure gradient and the log-layer mismatch
(U₂ − U₁ against the log law, relative to U₂; Yang et al.) from the mean profile over the second
half. T_f = 4 s is about 4 times the time the flow takes to cross a wall cell.

| wall model | SGS | Re_τ vs Dean | τ_w / Dean | mismatch |
|---|---|---|---|---|
| `nutUSpaldingWallFunction` (upstream) | WALE | −3.3 % | 0.935 | +7.7 % |
| `sampleLayer 1` (same as upstream) | WALE | −3.7 % | 0.927 | +8.0 % |
| `Tf 4` | WALE | −2.4 % | 0.953 | +5.2 % |
| `Tf 15` | WALE | −2.6 % | 0.949 | +4.6 % |
| `sampleLayer 2` | WALE | +1.4 % | 1.028 | +6.3 % |
| `sampleLayer 3` | WALE | −2.1 % | 0.959 | +5.2 % |
| `sampleLayer 2; Tf 4` | WALE | +0.2 % | 1.003 | +4.8 % |
| `sampleLayer 1` | σ | −6.1 % | 0.881 | +10.2 % |
| `sampleLayer 2; Tf 4` | σ | −0.6 % | 0.989 | +7.7 % |

The sample layer moves the wall stress most, the time filter adds a little, and together they
remove the friction deficit; σ, whose eddy viscosity vanishes in near-wall shear, makes the
mismatch worse with the wall cell as input and needs the remedy more. Differences of 1–2 % are
within the scatter of a run of this length. `tutorials/wallFunctions/channelWMLES/Allrun -full`
reruns the `sampleLayer 2; Tf 4` / WALE case and fails when Re_τ is more than 3 % off.

References: Kawai & Larsson (2012), Phys. Fluids 24, 015105. Yang, Park & Moin (2017), Phys. Rev.
Fluids 2, 104601. Spalding (1961), J. Appl. Mech. 28, 455.

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
  is more than 3 % off the DNS of Moser, Kim & Mansour (1999).
- `Test-wallLaw` checks the laws of the wall used by the wall functions.
- `Test-sigma` checks the σ model's operator; `tutorials/sigma/channel395/Allrun [steps]` runs
  the wall-resolved channel with σ and fails when it was not selected or gives no subgrid
  viscosity.
- `tutorials/wallFunctions/channelWMLES/Allrun [steps | -full]` runs the wall-modelled LES channel
  with `nutSpaldingLESWallFunction` on 2 processors, each holding one of the walls, and fails
  when the wall function did not sample both walls or write its filtered velocity; `-full` (about 45 min) also fails when Re_τ is more than 3 %
  off Dean's correlation.
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
