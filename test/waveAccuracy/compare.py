#!/usr/bin/env python3
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Compares the interface elevation of interfaceHeight and surfacePosition with
# the exact stream function wave, per line offset within the cells, and the
# serial with the parallel run.
#
# Usage: compare.py <case> [<case> ...]
# A case <name> is compared with <name>_parallel if that exists.

import json
import os
import re
import sys

import numpy as np

METHODS = [
    ("interfaceHeight cellPoint", "interfaceHeightCellPoint", "height.dat"),
    ("interfaceHeight cell", "interfaceHeightCell", "height.dat"),
    ("surfacePosition", "surfacePosition", "height.dat"),
]


def table(case, fo, fname):
    return np.genfromtxt(os.path.join(case, "postProcessing", fo, "0", fname), comments="#")


def elevation(case, fo, fname, nLines, zBottom):
    d = table(case, fo, fname)
    t = d[:, 0]
    # Columns per line: hB (height above the lower boundary), hL
    z = d[:, 1:1 + 2*nLines:2] + zBottom
    return t, z


def lines(case):
    text = open(os.path.join(case, "system", "lineLocations")).read()
    info = re.findall(r"^//\s+(\d+):\s+(\d+)\s+([\d.]+)", text, re.M)
    return [(int(i), float(f)) for _, i, f in info]


def zBottom(case):
    text = open(os.path.join(case, "system", "blockMeshDict")).read()
    return float(re.search(r"^\s*z1\s+([-\d.eE+]+)\s*;", text, re.M).group(1))


def compare(case):
    info = lines(case)
    n = len(info)
    zb = zBottom(case)
    ex = np.genfromtxt(os.path.join(case, "postProcessing", "exactSurface", "exactSurface.dat"),
                       comments="#")
    exact = {round(t, 9): row for t, row in zip(ex[:, 0], ex[:, 1:])}
    offsets = sorted({f for _, f in info}, reverse=True)

    print(f"\n=== {case}: {n} lines; error of the interface elevation against the exact wave"
          " [mm]: rms, mean (max abs) over one period")
    print(f"{'method':28s}" + "".join(f"{'offset %.2f' % f:>24s}" for f in offsets) + f"{'all lines':>24s}")

    results = {}
    summary = {"lines": n, "offsets": offsets, "samples": 0, "methods": {}, "parallel": {}}
    for label, fo, fname in METHODS:
        t, z = elevation(case, fo, fname, n, zb)
        keep = [i for i, tt in enumerate(t) if round(tt, 9) in exact]
        ze = np.array([exact[round(t[i], 9)] for i in keep])
        err = (z[keep] - ze)*1e3
        results[label] = z
        summary["samples"] = len(keep)
        stats = summary["methods"][label] = {}
        cells = []
        for f in offsets + [None]:
            cols = [k for k, (_, ff) in enumerate(info) if f is None or ff == f]
            e = err[:, cols]
            stats["all" if f is None else f"{f:.2f}"] = {
                "rms": float(np.sqrt(np.nanmean(e**2))),
                "mean": float(np.nanmean(e)),
                "max": float(np.nanmax(np.abs(e))),
                "nan": int(np.isnan(e).sum()),
            }
            cells.append(f"{np.sqrt(np.nanmean(e**2)):6.3f} {np.nanmean(e):+6.3f} ({np.nanmax(np.abs(e)):6.3f})")
        print(f"{label:28s}" + "".join(f"{c:>24s}" for c in cells))

    par = case.rstrip("/") + "_parallel"
    if os.path.isdir(par):
        print(f"--- serial against parallel ({par}): max |dz| [m]")
        for label, fo, fname in METHODS:
            _, zp = elevation(par, fo, fname, n, zb)
            m = min(len(zp), len(results[label]))
            dz = float(np.nanmax(np.abs(zp[:m] - results[label][:m])))
            summary["parallel"][label] = dz
            print(f"    {label:28s} {dz:.3e}")

    return summary


out = {c: compare(c) for c in sys.argv[1:]}
with open("results.json", "w") as f:
    json.dump(out, f, indent=1)
