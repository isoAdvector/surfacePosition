#!/usr/bin/env python3
#
# Checks of surfacePosition on changing meshes:
#  - the line cells after each rebuild against a full search (log lines of
#    checkLineCells),
#  - the interface elevation against interfaceHeight, and that surfacePosition
#    finds an interface on every line at every time (no nan),
#  - serial against parallel.
#
# Usage: compare.py <case> [<case> ...]   (<case>_parallel is used if present)

import json
import os
import re
import sys

import numpy as np


def rebuilds(log):
    text = open(log).read()
    rows = re.findall(r"^surfacePosition (\S+): rebuild (\d+) at time \S+: line cells (\d+), missing (\d+), extra (\d+)",
                      text, re.M)
    out = {}
    for name, _, cells, missing, extra in rows:
        n, c, m, e = out.get(name, (0, 0, 0, 0))
        out[name] = (n + 1, c + int(cells), m + int(missing), e + int(extra))
    return out, len(re.findall(r"^Time = ", text, re.M))


def elevation(case, fo, fname, nLines):
    p = os.path.join(case, "postProcessing", fo, "0", fname)
    d = np.genfromtxt(p, comments="#")
    if fname == "height.dat":
        return d[:, 0], d[:, 1:1 + 2*nLines:2]          # hB above z = 0
    return d[:, 0], d[:, 5:1 + 5*nLines:5]


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
        t, ih = elevation(run, "interfaceHeight1", "height.dat", n)
        for fo in ("surfacePosition1", "surfacePositionPlic1"):
            ts, sp = elevation(run, fo, "surfacePosition.dat", n)
            m = min(len(t), len(ts))
            d = (sp[:m] - ih[:m])*1e3
            nan = np.argwhere(np.isnan(sp))
            summary["elevation"][fo] = {
                "mean": float(np.nanmean(d)), "rms": float(np.sqrt(np.nanmean(d**2))),
                "max": float(np.nanmax(np.abs(d))), "steps": m, "nan": len(nan),
            }
            print(f"  {fo} - interfaceHeight1: mean {np.nanmean(d):+.3f} mm, rms {np.sqrt(np.nanmean(d**2)):.3f} mm,"
                  f" max {np.nanmax(np.abs(d)):.3f} mm over {m} steps")
            if len(nan):
                lines = sorted(set(nan[:, 1]))
                print(f"  {fo}: NO INTERFACE in {len(nan)} outputs, lines {lines}, from time {ts[nan[0, 0]]}")
                failures.append(f"{run} {fo}: no interface")

    par = case + "_parallel"
    if os.path.isdir(par):
        print(f"--- {case} serial against parallel at the serial times, max |dz|:")
        for fo, fname in (("surfacePosition1", "surfacePosition.dat"),
                          ("surfacePositionPlic1", "surfacePosition.dat"),
                          ("interfaceHeight1", "height.dat")):
            ts, zs = elevation(case, fo, fname, n)
            tp, zp = elevation(par, fo, fname, n)
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
