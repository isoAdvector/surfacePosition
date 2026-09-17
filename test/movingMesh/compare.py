#!/usr/bin/env python3
#
# Checks of surfacePosition on changing meshes:
#  - the line cells after each rebuild against a full search (log lines of
#    checkLineCells),
#  - the interface elevation against interfaceHeight, and that surfacePosition
#    finds every line in the mesh at every time (hL = GREAT marks a line
#    missing the mesh),
#  - serial against parallel.
#
# Usage: compare.py <case> [<case> ...]   (<case>_parallel is used if present)

import json
import os
import re
import sys

import numpy as np

GREAT = 1e15        # OpenFOAM's GREAT (double precision)


def rebuilds(log):
    text = open(log).read()
    rows = re.findall(r"^surfacePosition (\S+): rebuild (\d+) at time \S+: line cells (\d+), missing (\d+), extra (\d+)",
                      text, re.M)
    out = {}
    for name, _, cells, missing, extra in rows:
        n, c, m, e = out.get(name, (0, 0, 0, 0))
        out[name] = (n + 1, c + int(cells), m + int(missing), e + int(extra))
    return out, len(re.findall(r"^Time = ", text, re.M))


def heights(case, fo, nLines):
    """Times, hB (above the bottom z = 0) and hL per line"""
    p = os.path.join(case, "postProcessing", fo, "0", "height.dat")
    d = np.genfromtxt(p, comments="#")
    return d[:, 0], d[:, 1:1 + 2*nLines:2], d[:, 2:2 + 2*nLines:2]


def nLines(case):
    text = open(os.path.join(case, "system", "lineLocations")).read()
    return len(re.findall(r"^\s+\(", text, re.M))


failures = []
results = {}

for case in sys.argv[1:]:
    n = nLines(case)
    for run in (case, case + "_parallel"):
        if not os.path.isdir(run):
            continue
        app = re.search(r"^application\s+(\S+);", open(os.path.join(run, "system", "controlDict")).read(), re.M).group(1)
        stats, steps = rebuilds(os.path.join(run, f"log.{app}"))
        summary = results[run] = {"steps": steps, "rebuilds": {}, "elevation": {}}
        print(f"\n=== {run}: {steps} time steps")
        for name, (k, c, m, e) in stats.items():
            status = "OK" if m == 0 and e == 0 else "MISMATCH"
            summary["rebuilds"][name] = {"rebuilds": k, "cells": c, "missing": m, "extra": e}
            if status != "OK":
                failures.append(f"{run} {name}: line cells")
            print(f"  {name}: {k} rebuilds after mesh changes, {c} line cells in total,"
                  f" missing {m}, extra {e}: {status}")
        t, ih, _ = heights(run, "interfaceHeight1", n)
        for fo in ("surfacePosition1",):
            ts, sp, hL = heights(run, fo, n)
            m = min(len(t), len(ts))
            d = (sp[:m] - ih[:m])*1e3
            nan = np.argwhere(hL >= 0.5*GREAT)
            summary["elevation"][fo] = {
                "mean": float(np.nanmean(d)), "rms": float(np.sqrt(np.nanmean(d**2))),
                "max": float(np.nanmax(np.abs(d))), "steps": m, "missing": len(nan),
            }
            print(f"  {fo} - interfaceHeight1: mean {np.nanmean(d):+.3f} mm, rms {np.sqrt(np.nanmean(d**2)):.3f} mm,"
                  f" max {np.nanmax(np.abs(d)):.3f} mm over {m} steps")
            if len(nan):
                lines = sorted(set(nan[:, 1]))
                print(f"  {fo}: LINE MISSING THE MESH in {len(nan)} outputs, lines {lines}, from time {ts[nan[0, 0]]}")
                failures.append(f"{run} {fo}: line missing the mesh")

    par = case + "_parallel"
    if os.path.isdir(par):
        print(f"--- {case} serial against parallel at the serial times, max |dz|:")
        for fo in ("surfacePosition1", "interfaceHeight1"):
            ts, zs, _ = heights(case, fo, n)
            tp, zp, _ = heights(par, fo, n)
            keep = (ts >= tp[0]) & (ts <= tp[-1])
            d = [np.abs(zs[keep, k] - np.interp(ts[keep], tp, zp[:, k])) for k in range(n)]
            results.setdefault(case + " parallel", {})[fo] = float(np.nanmax(d))
            print(f"    {fo:22s} {np.nanmax(d):.3e} m  (time steps: serial {len(ts)}, parallel {len(tp)})")

with open("results.json", "w") as f:
    json.dump({"passed": not failures, "failures": failures, "runs": results}, f, indent=1)

if failures:
    print("\nSURFACEPOSITION CHECKS FAILED: " + "; ".join(failures))
    sys.exit(1)
print("\nALL SURFACEPOSITION CHECKS PASSED")
