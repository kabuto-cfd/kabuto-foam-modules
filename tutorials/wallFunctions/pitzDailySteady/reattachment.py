#!/usr/bin/env python3
"""Reattachment length behind the pitzDaily step (h = 25.4 mm): where the
x-velocity of the first cell row above the lower wall turns from negative to
positive, last such crossing before x = 0.2 m. Reads the latest time of
<case> written in ascii with cell centres (writeCellCentres)."""
import os, re, sys

case = sys.argv[1]
h = 0.0254


def read(path):
    s = open(path).read()
    m = re.search(r'internalField\s+nonuniform\s+List<\w+>\s*\n?(\d+)\s*\n?'
                  r'\((.*?)\n\)\s*;', s, re.S)
    return [tuple(map(float, t.split()))
            for t in re.findall(r'\(([^()]*)\)', m.group(2))]


t = max((d for d in os.listdir(case) if re.fullmatch(r'\d+', d)), key=int)
C, U = read(f'{case}/{t}/C'), read(f'{case}/{t}/U')
row = {}
for c, u in zip(C, U):
    if 0 < c[0] < 0.2 and c[1] < -h + 0.004:
        x = round(c[0], 6)
        if x not in row or c[1] < row[x][0]:
            row[x] = (c[1], u[0])
xs = sorted(row)
xr = None
for a, b in zip(xs, xs[1:]):
    if row[a][1] < 0 <= row[b][1]:
        xr = a + (b - a)*(-row[a][1])/(row[b][1] - row[a][1])
if xr is None:
    print('no reattachment')
    sys.exit(1)
print(f'reattachment x_r = {1000*xr:.1f} mm = {xr/h:.2f} h (iteration {t})')
