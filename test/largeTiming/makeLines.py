#!/usr/bin/env python3
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Writes system/lineLocations with a grid of vertical sampling lines in the
# cube of the large timing case, placed off the cell centres.
#
# Usage: makeLines.py <caseDir> [nx ny]   (default 5 x 10 = 50 lines)

import re
import sys

case = sys.argv[1]
nLx, nLy = (int(sys.argv[2]), int(sys.argv[3])) if len(sys.argv) > 3 else (5, 10)

text = open(f"{case}/system/blockMeshDict").read()
L = float(re.search(r"^\s*L\s+([-\d.eE+]+)\s*;", text, re.M).group(1))
n = int(re.search(r"^\s*n\s+(\d+)\s*;", text, re.M).group(1))
dx = L/n

points = []
for i in range(nLx):
    for j in range(nLy):
        # Cell columns spread over the cube, positions 0.3 and 0.7 cells off
        ci = round((i + 0.5)*n/nLx)
        cj = round((j + 0.5)*n/nLy)
        points.append(((ci + 0.3)*dx, (cj + 0.7)*dx, 0.5*L))

with open(f"{case}/system/lineLocations", "w") as out:
    out.write("// Written by makeLines.py: %d x %d lines, n = %d\n" % (nLx, nLy, n))
    out.write("lineLocations\n(\n")
    for p in points:
        out.write("    (%.12g %.12g %.12g)\n" % p)
    out.write(");\n")

print(f"{case}: {len(points)} lines, n = {n}")
