#!/usr/bin/env python3
"""Compare the result of a slab case with its exact solution; exit 1 if off.

    check.py <directory> twoLayer | closed | exchange | convection | source
"""
import math
import os
import re
import sys

case, kind = sys.argv[1], sys.argv[2]
times = [d for d in os.listdir(case) if re.fullmatch(r"[0-9.e+-]+", d) and float(d) > 0]
last = max(times, key=float)
text = open(os.path.join(case, last, "C")).read()
m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\(([^)]*)\)", text)
if m:
    C = [float(v) for v in m.group(2).split()]
else:
    C = [float(re.search(r"internalField\s+uniform\s+(\S+);", text).group(1))] * 100
n = len(C)
x = [(i + 0.5) / n for i in range(n)]

if kind == "twoLayer":
    Ci = 1e-10 / (1e-10 + 1e-11)
    exact = [1 - (1 - Ci) * xi / 0.5 if xi < 0.5 else Ci * (1 - xi) / 0.5 for xi in x]
    err = max(abs(a - b) for a, b in zip(C, exact))
    print("twoLayer: interface cells %.5f %.5f (exact %.5f %.5f), max error %.2e"
          % (C[n // 2 - 1], C[n // 2], exact[n // 2 - 1], exact[n // 2], err))
    ok = err < 1e-6
elif kind == "closed":
    left, right = min(C[: n // 2]), max(abs(c) for c in C[n // 2:])
    print("closed: min of the open half %.5f (-> 1), max of the closed half %.2e (-> 0)" % (left, right))
    ok = left > 0.999 and right < 1e-9
elif kind == "convection":
    exact = [1 - (math.exp(5 * xi) - 1) / (math.exp(5) - 1) for xi in x]
    err = max(abs(a - b) for a, b in zip(C, exact))
    print("convection: C at mid length %.5f (exact %.5f), max error %.2e" % (C[n // 2], exact[n // 2], err))
    ok = err < 1e-3
elif kind == "source":
    err = max(abs(ci - 1) for ci in C)
    print("source: C = %.6f (exact 1)" % C[0])
    ok = err < 1e-6
else:
    exact = 1 - math.exp(-1.0)
    err = max(abs(c / exact - 1) for c in C)
    print("exchange: C = %.5f (exact %.5f, %+.3f %%)" % (C[0], exact, 100 * (C[0] / exact - 1)))
    ok = err < 2e-3
sys.exit(0 if ok else 1)
