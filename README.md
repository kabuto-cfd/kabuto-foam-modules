# kabuto-foam-modules

Add-on modules for [OpenFOAM-14](https://openfoam.org) (OpenFOAM Foundation release) by Kabuto CFD.
Every module builds with `wmake` against an unmodified OpenFOAM-14 and is loaded at run time
through `libs (...)`; nothing in OpenFOAM itself is changed.

## Modules

| Module | Library | What it does |
|---|---|---|
| [`realizableR`](src/fvConstraints/realizableR/realizableR.H) | `libkabutoFvConstraints.so` | fvConstraint that keeps a Reynolds-stress field realizable: in cells where the smallest eigenvalue of **R** drops below a floor, the eigenvalues are clipped to that floor and **R** is rebuilt from its unchanged eigenvectors |

Planned: ω-based Reynolds-stress models, non-equilibrium wall functions, and solver modules
for `foamRun`.

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

GitHub Actions builds every module against the openfoam.org OpenFOAM-14 package on
Ubuntu 24.04 and runs both for every change.

## Licence

GNU General Public License v3 or later (`COPYING`), as OpenFOAM itself.
Copyright (C) 2026 Kabuto CFD.

## Trademarks

This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the
OpenFOAM software via www.openfoam.com, and owner of the OPENFOAM® and OpenCFD® trade marks,
nor by the OpenFOAM Foundation.
