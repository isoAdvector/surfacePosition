#!/usr/bin/env python3
#
# Time per step of the flow solution and of each function object in the large
# timing case, from OpenFOAM's profiling written with the fields, and
# agreement of the interface elevations between the function objects.
#
# The profiling times are cumulative: the difference between the first and
# the last written time leaves out the steps before the first write (start-up
# and line search). Each processor's time includes waiting for the others;
# the mean over the processors is given.
#
# Usage: timing.py <caseDir>

import json
import os
import re
import sys

import numpy as np

FUNCTION_OBJECTS = ["surfacePosition1", "interfaceHeight1"]

case = sys.argv[1]


def profiling(path):
    """Total time and number of calls of each profiling entry, by description"""
    out = {}
    for block in re.findall(r"trigger\d+\s*\{(.*?)\}", open(path).read(), re.S):
        name = re.search(r'description\s+"(.*?)"', block).group(1)
        total = float(re.search(r"totalTime\s+(\S+);", block).group(1))
        calls = int(re.search(r"calls\s+(\d+);", block).group(1))
        out[name] = (total, calls)
    return out


def entry(prof, name):
    for key, value in prof.items():
        if key.startswith(name):
            return value
    return (0.0, 0)


procs = sorted(d for d in os.listdir(case) if re.fullmatch(r"processor\d+", d))
times = sorted(
    (d for d in os.listdir(os.path.join(case, procs[0]))
     if os.path.isfile(os.path.join(case, procs[0], d, "uniform", "profiling"))),
    key=float
)
if len(times) < 2:
    sys.exit(f"Need profiling at two written times in {case}/{procs[0]}, found {times}")

text = open(os.path.join(case, "system", "lineLocations")).read()
nLines = len(re.findall(r"^\s+\(", text, re.M))
m = re.search(r"nCells:\s*(\d+)", open(os.path.join(case, "log.blockMesh")).read())
nCells = int(m.group(1)) if m else None

rows = {"flow solution": "time.run()"}
rows.update({fo: f"functionObject::{fo}:write" for fo in FUNCTION_OBJECTS})

cost = {row: [] for row in rows}
start = {row: [] for row in rows}
nSteps = None
for proc in procs:
    first = profiling(os.path.join(case, proc, times[0], "uniform", "profiling"))
    last = profiling(os.path.join(case, proc, times[-1], "uniform", "profiling"))
    for row, name in rows.items():
        t0, c0 = entry(first, name)
        t1, c1 = entry(last, name)
        if row == "flow solution":
            nSteps = c1 - c0
        cost[row].append((t1 - t0)/max(c1 - c0, 1))
        start[row].append(t0 - c0*cost[row][-1])

print(f"Case {case}: {nCells} cells, {len(procs)} processors, {nLines} lines,"
      f" {nSteps} steps timed (times {times[0]} to {times[-1]})")
print("ms/step: mean time per step over the processors, spread: its range over"
      " the processors,\nstart-up: time of the steps before the first write"
      " beyond ms/step (line search)")
summary = {
    "cells": nCells, "processors": len(procs), "lines": nLines, "steps": nSteps,
    "times": [float(times[0]), float(times[-1])], "cost": {}, "elevation": {},
}
print(f"{'':22s} {'ms/step':>9s} {'ms/line':>9s} {'spread':>9s} {'start-up [ms]':>14s}")
for row, values in cost.items():
    v = np.array(values)*1e3
    perLine = "" if row == "flow solution" else f"{v.mean()/nLines:9.3f}"
    summary["cost"][row] = {
        "msPerStep": float(v.mean()), "spread": float(v.max() - v.min()),
        "startup": float(np.mean(start[row])*1e3),
    }
    print(f"{row:22s} {v.mean():9.2f} {perLine:>9s} {v.max() - v.min():9.2f}"
          f" {np.mean(start[row])*1e3:14.1f}")


def elevations(fo, fname):
    """Times, hB (above the bottom z = 0) and hL per line"""
    p = os.path.join(case, "postProcessing", fo, "0", fname)
    if not os.path.isfile(p):
        return None
    d = np.genfromtxt(p, comments="#")
    return d[:, 0], d[:, 1:1 + 2*nLines:2], d[:, 2:2 + 2*nLines:2]


ih = elevations("interfaceHeight1", "height.dat")
for fo in ("surfacePosition1",):
    sp = elevations(fo, "height.dat")
    if ih is None or sp is None:
        continue
    m = min(len(ih[0]), len(sp[0]))
    d = (sp[1][:m] - ih[1][:m])*1e3
    summary["elevation"][fo] = {
        "mean": float(np.nanmean(d)), "rms": float(np.sqrt(np.nanmean(d**2))),
        "max": float(np.nanmax(np.abs(d))), "steps": m, "missing": int((sp[2] >= 0.5e15).sum()),
    }
    print(f"{fo} - interfaceHeight1: mean {np.nanmean(d):+.4f} mm, rms {np.sqrt(np.nanmean(d**2)):.4f} mm,"
          f" max {np.nanmax(np.abs(d)):.4f} mm over {m} steps and {nLines} lines")

with open("results.json", "w") as f:
    json.dump(summary, f, indent=1)
