#!/usr/bin/env python3
"""Write a plane-channel case (half channel h = 1, Ubar = 1 held by
meanVelocityForce, nu = 2.5e-5: Re_b = 80 000, Re_tau = 1864 from Dean's
correlation Cf = 0.073 Re_b^-1/4) whose first cell centre sits at a given y+.

    makeCase.py <dir> <kEpsilon|kOmegaSST|...> <treatment> <y+>

treatment: std    nutkWallFunction + epsilon/omegaWallFunction (OpenFOAM)
           ne     nutNonEquilibriumWallFunction + ...NonEquilibriumWallFunction
           nes    nutNonEquilibriumSpaldingWallFunction + ...Spalding...
           nes2l  nes + the twoLayerEpsilon constraint (epsilon models)
The cells grow geometrically from the wall by at most 1.12.
"""
import sys, os, math

d, model, wf, yp = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4])
twoLayer = wf.endswith('2l'); wf = wf[:-2] if twoLayer else wf
nu = 2.5e-5; ut = 0.0466
h1 = 2*yp*nu/ut
# geometric grading: smallest N with ratio <= 1.12
for N in range(2, 2000):
    if h1*N >= 1: r = 1.0; break
    lo, hi = 1.0, 3.0
    for _ in range(200):
        m = (lo+hi)/2
        s = h1*(m**N-1)/(m-1)
        lo, hi = (m, hi) if s < 1 else (lo, m)
    r = (lo+hi)/2
    if r <= 1.12: break
if h1*N >= 1: N = max(1, round(1/h1))
for sub in ['0', 'constant', 'system']: os.makedirs(f'{d}/{sub}', exist_ok=True)
H = lambda c, o: f'FoamFile {{ format ascii; class {c}; object {o}; }}\n'
def w(p, s): open(f'{d}/{p}', 'w').write(s)
w('system/blockMeshDict', H('dictionary','blockMeshDict') + f'''convertToMeters 1;
vertices ( (0 0 0) (0.1 0 0) (0.1 1 0) (0 1 0) (0 0 0.1) (0.1 0 0.1) (0.1 1 0.1) (0 1 0.1) );
blocks ( hex (0 1 2 3 4 5 6 7) (1 {N} 1) simpleGrading (1 {r**(N-1) if N>1 else 1} 1) );
boundary (
 inlet {{ type cyclic; neighbourPatch outlet; faces ((0 4 7 3)); }}
 outlet {{ type cyclic; neighbourPatch inlet; faces ((1 2 6 5)); }}
 wall {{ type wall; faces ((0 1 5 4)); }}
 top {{ type symmetryPlane; faces ((3 7 6 2)); }}
 frontAndBack {{ type empty; faces ((0 3 2 1) (4 5 6 7)); }} );
''')
libs = '' if wf == 'std' else 'libs ("libkabutoWallFunctions.so");\n'
w('system/controlDict', H('dictionary','controlDict') + 'solver incompressibleFluid;\nstartFrom startTime; startTime 0; stopAt endTime; endTime 6000; deltaT 1;\nwriteControl timeStep; writeInterval 6000; writeFormat ascii; writePrecision 8;\n' + libs)
w('system/fvSchemes', H('dictionary','fvSchemes') + '''ddtSchemes { default steadyState; }
gradSchemes { default Gauss linear; }
divSchemes { default none; div(phi,U) bounded Gauss linearUpwind grad(U);
 div(phi,k) bounded Gauss upwind; div(phi,epsilon) bounded Gauss upwind; div(phi,omega) bounded Gauss upwind;
 div((nuEff*dev2(T(grad(U))))) Gauss linear; }
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes { default corrected; }
wallDist { method meshWave; }
''')
w('system/fvSolution', H('dictionary','fvSolution') + '''solvers { p { solver PCG; preconditioner DIC; tolerance 1e-10; relTol 0.01; }
 "(U|k|epsilon|omega)" { solver smoothSolver; smoother symGaussSeidel; tolerance 1e-10; relTol 0.1; } }
SIMPLE { nNonOrthogonalCorrectors 0; consistent yes; pRefCell 0; pRefValue 0;
 residualControl { p 1e-8; U 1e-9; k 1e-9; epsilon 1e-9; omega 1e-9; } }
relaxationFactors { equations { U 0.9; k 0.7; epsilon 0.7; omega 0.7; } }
''')
w('system/fvConstraints', H('dictionary','fvConstraints') + 'momentumSource { type meanVelocityForce; cellZone all; Ubar (1 0 0); }\n' + ('twoLayerEpsilon { type twoLayerEpsilon; libs ("libkabutoWallFunctions.so"); }\n' if twoLayer else ''))
w('constant/physicalProperties', H('dictionary','physicalProperties') + f'viscosityModel constant; nu {nu};\n')
w('constant/momentumTransport', H('dictionary','momentumTransport') + f'simulationType RAS;\nRAS {{ model {model}; turbulence on; }}\n')
bc = lambda wall: 'boundaryField { inlet { type cyclic; } outlet { type cyclic; } top { type symmetryPlane; } frontAndBack { type empty; } wall { ' + wall + ' } }\n'
nutwf = {'std': 'nutkWallFunction', 'ne': 'nutNonEquilibriumWallFunction', 'nes': 'nutNonEquilibriumSpaldingWallFunction'}[wf]
sfx = {'std': 'WallFunction', 'ne': 'NonEquilibriumWallFunction', 'nes': 'NonEquilibriumSpaldingWallFunction'}[wf]
w('0/U', H('volVectorField','U') + 'dimensions [0 1 -1 0 0]; internalField uniform (1 0 0);\n' + bc('type noSlip;'))
w('0/p', H('volScalarField','p') + 'dimensions [0 2 -2 0 0]; internalField uniform 0;\n' + bc('type zeroGradient;'))
w('0/k', H('volScalarField','k') + 'dimensions [0 2 -2 0 0]; internalField uniform 0.004;\n' + bc('type kqRWallFunction; value uniform 0.004;'))
w('0/epsilon', H('volScalarField','epsilon') + 'dimensions [0 2 -3 0 0]; internalField uniform 4e-4;\n' + bc(f'type epsilon{sfx}; value uniform 4e-4;'))
w('0/omega', H('volScalarField','omega') + 'dimensions [0 0 -1 0 0]; internalField uniform 1.1;\n' + bc(f'type omega{sfx}; value uniform 1.1;'))
w('0/nut', H('volScalarField','nut') + 'dimensions [0 2 -1 0 0]; internalField uniform 0;\n' + bc(f'type {nutwf}; value uniform 0;'))
print(N, r)
