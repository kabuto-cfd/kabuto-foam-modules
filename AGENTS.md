# kabuto-foam-modules — notes for AI coding agents

Public repository of Kabuto CFD: add-on libraries for OpenFOAM-14 (OpenFOAM Foundation),
GPL-3.0-or-later. Everything here links OpenFOAM and is GPL.

## Scope and boundaries
- Modules build with `wmake` against an **unmodified** OpenFOAM-14 and are loaded with
  `libs (...)`. A change that needs OpenFOAM's own sources edited does not go here; it goes to
  a separate fork (`kabuto-foam-14`, not created yet). Solver modules for `foamRun` are fine here.
- Never copy code from Kabuto CFD's closed repositories (polymesher / Kabuto Meshing, kabuto-web)
  into this one, and never make those link anything from here: this repository is GPL.
- Target: OpenFOAM-14 only. Older versions are not tested or supported.

## Layout
- `src/<area>/` one library per area, `Make/files` with
  `LIB = $(FOAM_USER_LIBBIN)/libkabuto<Area>` (e.g. `libkabutoFvConstraints`). A new area is one
  more `wmake $targetType src/<area>` line in `Allwmake` (and `wclean` in `Allwclean`).
  - `src/fvConstraints/realizableR` Reynolds-stress realizability limiter (eigenvalue clipping)
  - `src/momentumTransportModels` RAS and LES models, `libkabutoMomentumTransportModels`: the
    model templates under `RAS/<model>/` and `LES/<model>/`, registered for incompressible and
    compressible flows by `incompressible/` and `compressible/` (one `makeRASModel(<model>)` or
    `makeLESModel(<model>)` line each). Cases load it with
    `libs ("libkabutoMomentumTransportModels.so");` in `system/controlDict`.
    - `RAS/stressOmega` Wilcox's stress-omega Reynolds-stress model (2006)
    - `LES/sigma` Nicoud et al.'s sigma SGS model (2011); its point-wise operator in
      `sigmaOperator.H` (inline, unit-tested by `test/sigma`)
  - `src/wallFunctions` wall functions, `libkabutoWallFunctions`: the point-wise laws in
    `wallLaw/wallLaw.H` (inline, unit-tested by `test/wallLaw`), two families of nut / epsilon /
    omega boundary conditions and the `twoLayerEpsilon` fvConstraint:
    - `*NonEquilibriumWallFunction` Kim & Choudhury (log layer, two-layer averages in the wall cell)
    - `*NonEquilibriumSpaldingWallFunction` y+-insensitive (Spalding's law with the same
      pressure-gradient correction); epsilon models add `twoLayerEpsilon` for cells inside y+ ~ 30
    - `nutSpaldingLESWallFunction` wall-stress model for wall-modelled LES: Spalding's law fed
      with a velocity sampled `sampleLayer` cells from the wall and / or filtered in time (`Tf`,
      default `auto`), the remedies for the log-layer mismatch; the filtered velocity is written
      with the field. Its optional `nutFloor` acts through the `wallModelledNut` fvConstraint
    The epsilon / omega classes derive from copies of OpenFOAM's epsilon / omegaWallFunction with
    `calculate` made virtual; the Spalding ones override only `calculate`. All walls of a case use
    one family (OpenFOAM handles the wall cells of all wall-function patches from one master patch).
- `test/<module>/` small `Test-<module>` applications that check the numerics without a case;
  exit status 1 on failure. Built by `Allwmake`, run by CI.
- `tutorials/<module>/<case>/Allrun` a short run that shows the module working. Prefer copying an
  upstream tutorial from `$FOAM_TUTORIALS` in `Allrun` and adding to it over committing a copy;
  `Allrun` must fail (non-zero exit) when the module did not act. Working copies (`case/`) are
  ignored by git. A case of our own (no upstream equivalent) is committed in full and its
  `Allrun` checks a number against a reference (`stressOmega/channel395`: Re_tau against the DNS
  of Moser, Kim & Mansour); its outputs are ignored by git.
- `.github/workflows/build.yml` builds against the openfoam.org `openfoam14` package (Ubuntu
  24.04) and runs the tests and tutorials. Markdown-only and `docs/` changes start no run.

## Building in the Claude workspace
The workspace cannot reach openfoam.org. Install the repackaged openfoam.org package from this
repository's release `openfoam14-runtime-<version>` into `/opt/openfoam14`, as
`docs/workspace-openfoam14.md` shows (seconds, not a build); then
```bash
export USER=${USER:-root}
set +e; source /opt/openfoam14/etc/bashrc; set -e
./Allwmake -j 2 && Test-realizableR && tutorials/realizableR/ductSecondaryFlow/Allrun 100
tutorials/stressOmega/ductSecondaryFlow/Allrun 300 && tutorials/stressOmega/channel395/Allrun
Test-wallLaw && tutorials/wallFunctions/channelYPlus/Allrun && tutorials/wallFunctions/pitzDailySteady/Allrun nes2l
Test-sigma && tutorials/sigma/channel395/Allrun && tutorials/wallFunctions/channelWMLES/Allrun
```

## Code conventions
- OpenFOAM's style (4-space indent, `// * * *` section banners, `Foam::fv::` qualified
  definitions, `TypeName`, run-time selection). Take the closest upstream class as the template
  and say so in the header.
- File header: the Kabuto CFD block of `src/fvConstraints/realizableR/realizableR.H`
  (Copyright (C) 2026 Kabuto CFD, GPL v3 or later). Keep upstream copyright lines for code
  derived from OpenFOAM.
- Every module documents its equations, its dictionary entries with defaults, and an example in
  the header's `Usage`, and gets a row in the README's table.
- Never name a product, library or type "OpenFOAM ..." (trade mark). "for OpenFOAM-14" in prose
  is fine.

## Commit identity
Every commit is authored as
```bash
git -c user.name="Kayne Sakuragi" \
    -c user.email="60303801+kayne-sakuragi@users.noreply.github.com" commit ...
```
(or `git config user.name` / `user.email` in the clone, not `--global`). Never commit under
another name or address: no real name, no personal mail address, and not the session's default
identity. This repository is public. Trailers such as `Co-Authored-By:` stay as the session
adds them.
