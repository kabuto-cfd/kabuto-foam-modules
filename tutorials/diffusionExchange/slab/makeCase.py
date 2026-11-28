#!/usr/bin/env python3
"""Write one of the three 1-D checks of the diffusionExchange solver module.

    makeCase.py <directory> twoLayer | closed | exchange

A slab of length 1 mm in x, 100 cells, C = 1 at x = 0 unless stated.
  twoLayer  D = 1e-10 in the left half, 1e-11 in the right half, C = 0 at
            x = L, steady state: the value at the interface is
            D1/(D1 + D2) (two resistances in series)
  closed    D = 1e-10 in the left half, 0 in the right half, zero gradient at
            x = L, 20000 s: the left half fills, the right half stays empty
  exchange  uniform D, kExchange = 1e-3 1/s, Cref = 1, zero gradient at both
            ends, 1000 s: C = 1 - exp(-kExchange t) everywhere
Pure Python, no third-party modules.
"""
import os
import sys

case, kind = sys.argv[1], sys.argv[2]
n = 100
head = ("FoamFile\n{\n    format      ascii;\n    class       %s;\n"
        "    object      %s;\n}\n\n")


def write(path, cls, body):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(head % (cls, os.path.basename(path)) + body)


def field(name, dims, internal, left, right):
    write(os.path.join(case, "0", name), "volScalarField",
          "dimensions      %s;\n\ninternalField   %s;\n\nboundaryField\n{\n"
          "    left\n    {\n        %s\n    }\n    right\n    {\n        %s\n    }\n"
          "    sides\n    {\n        type            empty;\n    }\n}\n"
          % (dims, internal, left, right))


def nonuniform(values):
    return "nonuniform List<scalar> %d(%s)" % (len(values), " ".join("%g" % v for v in values))


zg = "type            zeroGradient;"
fixed = "type            fixedValue;\n        value           uniform %g;"

write(os.path.join(case, "system", "blockMeshDict"), "dictionary",
      "convertToMeters 1e-3;\n\nvertices\n(\n    (0 0 0) (1 0 0) (1 0.1 0) (0 0.1 0)\n"
      "    (0 0 0.1) (1 0 0.1) (1 0.1 0.1) (0 0.1 0.1)\n);\n\n"
      "blocks\n(\n    hex (0 1 2 3 4 5 6 7) (%d 1 1) simpleGrading (1 1 1)\n);\n\n"
      "boundary\n(\n    left\n    {\n        type patch;\n        faces ((0 4 7 3));\n    }\n"
      "    right\n    {\n        type patch;\n        faces ((1 2 6 5));\n    }\n"
      "    sides\n    {\n        type empty;\n        faces ((0 1 5 4) (3 7 6 2) (0 3 2 1) (4 5 6 7));\n    }\n);\n" % n)

ddt, end, dt, props = "Euler", 20000, 100, "Cref            0;\n"
if kind == "twoLayer":
    ddt, end, dt = "steadyState", 1, 1
    field("D", "[0 2 -1 0 0 0 0]", nonuniform([1e-10] * (n // 2) + [1e-11] * (n // 2)), zg, zg)
    field("C", "[0 0 0 0 0 0 0]", "uniform 0", fixed % 1, fixed % 0)
elif kind == "closed":
    field("D", "[0 2 -1 0 0 0 0]", nonuniform([1e-10] * (n // 2) + [0] * (n // 2)), zg, zg)
    field("C", "[0 0 0 0 0 0 0]", "uniform 0", fixed % 1, zg)
elif kind == "exchange":
    end, dt = 1000, 1
    props = "Cref            1;\nD               1e-10;\nkExchange       1e-3;\n"
    field("C", "[0 0 0 0 0 0 0]", "uniform 0", zg, zg)
else:
    sys.exit("makeCase.py: unknown case " + kind)

write(os.path.join(case, "constant", "physicalProperties"), "dictionary", props)
write(os.path.join(case, "system", "controlDict"), "dictionary",
      "solver          diffusionExchange;\n\nstartFrom       startTime;\nstartTime       0;\n"
      "stopAt          endTime;\nendTime         %d;\ndeltaT          %d;\n\n"
      "writeControl    runTime;\nwriteInterval   %d;\nwriteFormat     ascii;\n"
      "writePrecision  10;\ntimeFormat      general;\nrunTimeModifiable no;\n" % (end, dt, end))
write(os.path.join(case, "system", "fvSchemes"), "dictionary",
      "ddtSchemes\n{\n    default         %s;\n}\n\ngradSchemes\n{\n    default         Gauss linear;\n}\n\n"
      "divSchemes\n{\n    default         none;\n}\n\n"
      "laplacianSchemes\n{\n    default         none;\n    laplacian(D,C)  Gauss harmonic corrected;\n}\n\n"
      "interpolationSchemes\n{\n    default         linear;\n}\n\nsnGradSchemes\n{\n    default         corrected;\n}\n" % ddt)
write(os.path.join(case, "system", "fvSolution"), "dictionary",
      "solvers\n{\n    \"C.*\"\n    {\n        solver          PCG;\n        preconditioner  DIC;\n"
      "        tolerance       1e-14;\n        relTol          0;\n    }\n}\n\n"
      "PIMPLE\n{\n    nNonOrthogonalCorrectors 0;\n}\n")
