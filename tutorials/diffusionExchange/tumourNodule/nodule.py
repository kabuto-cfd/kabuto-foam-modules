#!/usr/bin/env python3
"""Delivery of nanoparticles from the blood to a small tumour nodule in liver tissue:
an illustration of the diffusionExchange solver module fed by a pore-network lookup.

    nodule.py make <directory> [particle radius in nm, default 50]
    nodule.py check <directory>
    nodule.py reference [particle radius in nm]     (explicit scheme in pure Python)

Domain: a cube of 8 mm, 32^3 cells, a spherical nodule of radius 3 mm at its centre.

    ddt(C) - div(D grad(C)) = kExchange (Cref(t) - C)

  C         interstitial concentration, relative to the plasma concentration at t = 0
  Cref      plasma concentration, halved every 6 h after one dose
  D         D_free(R) * 0.915 * sigma(R/rbar): Stokes-Einstein in water at 37 C times the
            lookup of poreNetworkLookup.csv with the local mean throat radius rbar of the
            interstitium; zero where the particle is above the percolation threshold
  kExchange supply from the capillaries the mesh does not resolve, two resistances in series:
            the vessel wall, P S/V, and the diffusion over the half distance delta between
            capillaries, pi^2 D/(4 delta^2)

EVERYTHING ABOUT THE TISSUE IS ASSUMED, for illustration only; none of it is measured:
  rbar      20 nm at the centre of the nodule rising to 60 nm at its rim, 80 nm in the liver
            tissue around it, with a smooth +-25 % variation
  delta     15 um in the liver tissue, 40 um at the rim of the nodule, growing inwards; no
            vessels inside 40 % of the nodule radius
  P         power law through 1.2 nm/s at R = 3.5 nm and 0.2 nm/s at R = 50 nm, ten times
            that in the liver tissue; zero for particles above the wall pore radius (100 nm
            in the nodule, 60 nm in the liver tissue)
The lookup itself was compared with diffusion measurements in 1.5 % agarose gel only
(Fatin-Rouge, Starchev & Buffle, Biophys. J. 86, 2004), not with tissue. Interstitial flow,
uptake by cells and clearance are left out. Not for any clinical use.

Pure Python, no third-party modules.
"""
import math
import os
import re
import sys

N = 32
LENGTH = 8.0e-3
RT = 3.0e-3
DX = LENGTH/N
PLATEAU = 0.915
T_HALF = 6*3600.0
END = 86400
HERE = os.path.dirname(os.path.abspath(__file__))



def lookup():
    xs, ss, xc = [], [], None
    for line in open(os.path.join(HERE, "poreNetworkLookup.csv")):
        m = re.search(r"above x = ([0-9.]+)", line)
        if m:
            xc = float(m.group(1).rstrip("."))
        if line[0] in "#x":
            continue
        x, s = line.split(",")
        xs.append(float(x)); ss.append(float(s))
    return xs, ss, xc


def sigma(x, table):
    xs, ss, xc = table
    if x >= xc:
        return 0.0
    if x <= xs[0]:
        return ss[0]
    for i in range(1, len(xs)):
        if x <= xs[i]:
            return ss[i - 1] + (ss[i] - ss[i - 1])*(x - xs[i - 1])/(xs[i] - xs[i - 1])
    return 0.0


def fields(R):
    """Cell lists (x fastest, as blockMesh numbers a single block) of s, rbar, D, kExchange."""
    table = lookup()
    kT = 1.380649e-23*310.15
    D0 = kT/(6*math.pi*6.92e-4*R*1e-9)
    n = math.log(6.0)/math.log(50.0/3.5)
    perm = 1.2e-9*(R/3.5)**(-n)
    two_pi = 2*math.pi
    S, RB, D, K = [], [], [], []
    for k in range(N):
        z = (k + 0.5)*DX - LENGTH/2
        for j in range(N):
            y = (j + 0.5)*DX - LENGTH/2
            for i in range(N):
                x = (i + 0.5)*DX - LENGTH/2
                s = math.sqrt(x*x + y*y + z*z)/RT
                xr, yr, zr = x/LENGTH, y/LENGTH, z/LENGTH
                noise = 0.25/math.sqrt(1.5)*(
                    math.cos(two_pi*(1.3*xr + 0.7*yr))
                  + math.cos(two_pi*(0.9*yr - 1.1*zr) + 1.0)
                  + math.cos(two_pi*(1.7*zr + 0.4*xr) + 2.0))
                inside = s < 1
                rbar = (20.0 + 40.0*s*s if inside else 80.0)*math.exp(noise)
                d = D0*PLATEAU*sigma(R/rbar, table)
                delta = 40e-6*(1.0/min(max(s, 0.4), 1.0))**3 if inside else 15e-6
                sv = 0.0 if (inside and s < 0.4) else 2*5e-6/delta**2
                wall = 100.0 if inside else 60.0
                kw = (perm*(1.0 if inside else 10.0) if R < wall else 0.0)*sv
                kd = math.pi**2*d/(4*delta**2)
                S.append(s); RB.append(rbar); D.append(d)
                K.append(1.0/(1.0/kw + 1.0/kd) if kw > 0 and kd > 0 else 0.0)
    return S, RB, D, K


HEAD = "FoamFile\n{\n    format      ascii;\n    class       %s;\n    object      %s;\n}\n\n"


def write(path, cls, body):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(HEAD % (cls, os.path.basename(path)) + body)


def write_field(case, name, dims, values):
    internal = ("uniform %g" % values if isinstance(values, float) else
                "nonuniform List<scalar>\n%d\n(\n%s\n)\n" % (len(values), "\n".join("%.6e" % v for v in values)))
    write(os.path.join(case, "0", name), "volScalarField",
          "dimensions      %s;\n\ninternalField   %s;\n\nboundaryField\n{\n    outer\n    {\n"
          "        type            zeroGradient;\n    }\n}\n" % (dims, internal))


def make(case, R):
    S, RB, D, K = fields(R)
    write_field(case, "C", "[0 0 0 0 0 0 0]", 0.0)
    write_field(case, "D", "[0 2 -1 0 0 0 0]", D)
    write_field(case, "kExchange", "[0 0 -1 0 0 0 0]", K)
    write_field(case, "rbar", "[0 0 0 0 0 0 0]", RB)
    h = 1e3*LENGTH/2
    write(os.path.join(case, "system", "blockMeshDict"), "dictionary",
          "convertToMeters 1e-3;\n\nvertices\n(\n    (-%g -%g -%g) (%g -%g -%g) (%g %g -%g) (-%g %g -%g)\n"
          "    (-%g -%g %g) (%g -%g %g) (%g %g %g) (-%g %g %g)\n);\n\n"
          "blocks\n(\n    hex (0 1 2 3 4 5 6 7) (%d %d %d) simpleGrading (1 1 1)\n);\n\n"
          "boundary\n(\n    outer\n    {\n        type patch;\n        faces\n        (\n"
          "            (0 4 7 3) (1 2 6 5) (0 1 5 4) (3 7 6 2) (0 3 2 1) (4 5 6 7)\n        );\n    }\n);\n"
          % ((h,)*24 + (N, N, N)))
    table = " ".join("(%d %.10f)" % (t, 0.5**(t/T_HALF)) for t in range(0, END + 1, 300))
    write(os.path.join(case, "constant", "physicalProperties"), "dictionary",
          "// Plasma concentration: one dose at t = 0, half-life 6 h\n"
          "Cref\n{\n    type    table;\n    values  (%s);\n}\n" % table)
    write(os.path.join(case, "system", "controlDict"), "dictionary",
          "solver          diffusionExchange;\n\nstartFrom       startTime;\nstartTime       0;\n"
          "stopAt          endTime;\nendTime         %d;\ndeltaT          300;\n\n"
          "writeControl    runTime;\nwriteInterval   %d;\nwriteFormat     ascii;\n"
          "writePrecision  8;\ntimeFormat      general;\nrunTimeModifiable no;\n" % (END, END))
    write(os.path.join(case, "system", "fvSchemes"), "dictionary",
          "ddtSchemes\n{\n    default         Euler;\n}\n\ngradSchemes\n{\n    default         Gauss linear;\n}\n\n"
          "divSchemes\n{\n    default         none;\n}\n\n"
          "laplacianSchemes\n{\n    default         none;\n    laplacian(D,C)  Gauss harmonic corrected;\n}\n\n"
          "interpolationSchemes\n{\n    default         linear;\n}\n\nsnGradSchemes\n{\n    default         corrected;\n}\n")
    write(os.path.join(case, "system", "fvSolution"), "dictionary",
          "solvers\n{\n    \"C.*\"\n    {\n        solver          PCG;\n        preconditioner  DIC;\n"
          "        tolerance       1e-12;\n        relTol          0;\n    }\n}\n\n"
          "PIMPLE\n{\n    nNonOrthogonalCorrectors 0;\n}\n")
    print("R = %g nm: %d of %d cells closed (D = 0), %d without supply (kExchange = 0)"
          % (R, sum(d == 0 for d in D), N**3, sum(k == 0 for k in K)))


def averages(C, S, D):
    def mean(sel):
        v = [c for c, ok in zip(C, sel) if ok]
        return sum(v)/len(v)
    closed = [abs(c) for c, d in zip(C, D) if d == 0]
    return {"nodule": mean([s < 1 for s in S]), "rim": mean([0.7 < s < 1 for s in S]),
            "tissue": mean([s >= 1 for s in S]), "core": mean([s < 0.4 for s in S]), "closed": max(closed) if closed else 0.0}


def reference(R):
    """Explicit finite volume scheme, harmonic mean of D at the faces, zero flux outside."""
    S, RB, D, K = fields(R)
    C = [0.0]*N**3
    dt = min(0.15*DX*DX/max(D), 0.2/max(K), 600.0)
    steps = int(math.ceil(END/dt)); dt = END/steps
    nb = []
    for off, lim in ((1, lambda i, j, k: i < N - 1), (N, lambda i, j, k: j < N - 1), (N*N, lambda i, j, k: k < N - 1)):
        for k in range(N):
            for j in range(N):
                for i in range(N):
                    if lim(i, j, k):
                        a = i + N*(j + N*k); b = a + off
                        if D[a] > 0 and D[b] > 0:
                            nb.append((a, b, 2*D[a]*D[b]/(D[a] + D[b])/DX/DX))
    for n in range(steps):
        cv = 0.5**((n + 0.5)*dt/T_HALF)
        dC = [k*(cv - c) for k, c in zip(K, C)]
        for a, b, g in nb:
            f = g*(C[b] - C[a])
            dC[a] += f; dC[b] -= f
        C = [c + dt*d for c, d in zip(C, dC)]
    return averages(C, S, D)


def check(case):
    text = open(os.path.join(case, str(END), "C")).read()
    C = [float(v) for v in text.split("\n(\n", 1)[1].split("\n)\n", 1)[0].split()]
    S, RB, D, K = fields(50.0)
    a = averages(C, S, D)
    # Independent solution: explicit scheme, written without the module
    ref_all = reference(50.0)
    ok = a["closed"] < 1e-12 and sum(d == 0 for d in D) > 0
    print("closed cells: max |C| = %.1e" % a["closed"])
    for name in ("nodule", "rim", "tissue"):
        ref = ref_all[name]
        dev = 100*(a[name]/ref - 1)
        print("%-7s C = %.5f (reference %.5f, %+.2f %%)" % (name, a[name], ref, dev))
        ok = ok and abs(dev) < 3
    return ok


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "make":
        make(sys.argv[2], float(sys.argv[3]) if len(sys.argv) > 3 else 50.0)
    elif cmd == "check":
        sys.exit(0 if check(sys.argv[2]) else 1)
    elif cmd == "reference":
        print(reference(float(sys.argv[2]) if len(sys.argv) > 2 else 50.0))
    else:
        sys.exit(__doc__)
