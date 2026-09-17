#!/usr/bin/env python3
#
# Writes system/lineLocations for a periodic wave case: vertical sampling lines
# at a few cells across the wave length, each at the cell centre, a quarter of
# a cell off the centre and 2 % of a cell from a cell face, so that the effect
# of the sampling position relative to the cells can be seen.
#
# Usage: makeLines.py <caseDir>   (reads L and nx from system/blockMeshDict)

import re
import sys

CELLS = (3, 12, 22, 31)          # cell columns, as fractions of 38 cells
OFFSETS = (0.5, 0.25, 0.02)      # position within the cell, in cell widths
Y, Z = 0.0, 0.0                  # line location (z only sets the height origin


def entry(text, key):
    return float(re.search(r"^\s*%s\s+([-\d.eE+]+)\s*;" % key, text, re.M).group(1))


case = sys.argv[1]
text = open(f"{case}/system/blockMeshDict").read()
L, nx = entry(text, "L"), int(entry(text, "nx"))
dx = L/nx

lines = []
for c in CELLS:
    i = round(c*nx/38)
    for f in OFFSETS:
        lines.append(((i + f)*dx, i, f))

with open(f"{case}/system/lineLocations", "w") as out:
    out.write("// Written by makeLines.py: nx = %d, dx = %.12g\n" % (nx, dx))
    out.write("// line: cell column, offset in cell widths\n")
    for n, (x, i, f) in enumerate(lines):
        out.write("//   %2d: %3d %.2f\n" % (n, i, f))
    out.write("lineLocations\n(\n")
    for x, i, f in lines:
        out.write("    (%.12g %g %g)\n" % (x, Y, Z))
    out.write(");\n")

print(f"{case}: {len(lines)} lines, nx = {nx}, dx = {dx:.6g}")
