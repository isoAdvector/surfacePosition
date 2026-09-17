#!/usr/bin/env python3
#
# Compares the heights and positions written by surfacePosition (and
# interfaceHeight) in the special cases with the exact values, serial and
# parallel. Both write, per location, the height hB of the liquid above the
# lower boundary (z = 0 here), the height hL of the interface above the
# location, and the position of the interface; a line missing the mesh gives
# hB = 0 and hL = GREAT. surfacePosition is checked in its three modes:
# accumulated (liquid length), highest (top of the highest liquid) and lowest
# (top of the lowest liquid part).
#
# Usage: check.py <scenario> ...   (cases <scenario> and <scenario>_parallel)

import json
import os
import re
import sys

import numpy as np

TOL = 1e-8          # tolerance for the heights and positions [m]
GREAT = 1e15        # OpenFOAM's GREAT (double precision)
EPS = 1e-9


def locations(case, fo):
    """The locations from the header of the height file"""
    text = open(os.path.join(case, "postProcessing", fo, "0", "height.dat")).read()
    return [tuple(float(v) for v in m.split())
            for m in re.findall(r"^# Location \d+\s*:\s*\(([^)]*)\)", text, re.M)]


MODES = {"accumulated": "vertical", "highest": "vertical_highest", "lowest": "vertical_lowest"}


def expected_hB(scenario, loc, mode="accumulated"):
    """Exact height above the bottom in the mode, None for a line missing the
    mesh (the unit cube)"""
    x, y, _ = loc
    if not (-EPS <= x <= 1 + EPS and -EPS <= y <= 1 + EPS):
        return None
    if scenario.startswith("flat"):
        return int(scenario[4:])/1000
    if scenario == "layers":
        # Water below z = 0.35 and at 0.6 < z < 0.75
        return {"accumulated": 0.5, "highest": 0.75, "lowest": 0.35}[mode]
    # Slab of water for x < 0.5, cyclic in x: a line on the plane x = 0.5
    # lies between a full and an empty cell; their liquid lengths are
    # averaged, and for the other modes the line counts as in the liquid
    if abs(x - 0.5) < EPS:
        return 0.5 if mode == "accumulated" else 1.0
    return 1.0 if x < 0.5 else 0.0


def last_rows(case, fo):
    """Last row of the height and position files: (hB, hL) and (x, y, z) per
    location"""
    d = os.path.join(case, "postProcessing", fo, "0")
    if not os.path.isfile(os.path.join(d, "height.dat")):
        return None, None
    h = np.atleast_2d(np.genfromtxt(os.path.join(d, "height.dat"), comments="#"))[-1]
    lines = [l for l in open(os.path.join(d, "position.dat")) if not l.startswith("#")]
    p = [tuple(float(v) for v in m.split()) for m in re.findall(r"\(([^)]*)\)", lines[-1])]
    heights = [(h[1 + 2*i], h[2 + 2*i]) for i in range((len(h) - 1)//2)]
    return heights, p


def errors(scenario, locs, heights, points, mode="accumulated"):
    """Largest error of hB, hL and the position per line; None if the line
    should miss the mesh and correctly does, NaN if wrong"""
    out = []
    for loc, (hB, hL), p in zip(locs, heights, points):
        e = expected_hB(scenario, loc, mode)
        if e is None:
            out.append(None if (hB == 0 and hL >= 0.5*GREAT) else float("nan"))
            continue
        exact_p = (loc[0], loc[1], e)
        out.append(max(abs(hB - e), abs(hL - (e - loc[2])),
                       max(abs(a - b) for a, b in zip(p, exact_p))))
    return out


allOk = True
results = {}
for scenario in sys.argv[1:]:
    print(f"\n=== {scenario}")
    for run in (scenario, scenario + "_parallel"):
        if not os.path.isdir(run):
            continue
        print(f"--- {run}")
        summary = results[run] = {"checks": {}, "interfaceHeight": {}}

        locs = locations(run, "vertical")
        for mode, fo in MODES.items():
            heights, points = last_rows(run, fo)
            label = f"surfacePosition {mode}"
            if heights is None:
                print(f"  {label:28s} no output")
                allOk = False
                continue
            errs = errors(scenario, locs, heights, points, mode)
            bad = [i for i, e in enumerate(errs) if e is not None and not e <= TOL]
            maxErr = max((e for e in errs if e is not None and e <= TOL), default=0.0)
            if bad:
                allOk = False
            summary["checks"][mode] = {
                "maxError": maxErr, "nLines": len(errs), "failed": bad,
                "missing": sum(e is None for e in errs),
            }
            status = "ok" if not bad else "FAILED lines " + ", ".join(map(str, bad))
            print(f"  {label:28s} max error {maxErr:.2e}  {status}")

        # interfaceHeight on the same locations, for comparison
        ihHeights, ihPoints = last_rows(run, "vertical_interfaceHeight")
        if ihHeights is not None:
            ihErrs = errors(scenario, locs, ihHeights, ihPoints)
            differing = []
            for i, (loc, (hB, hL)) in enumerate(zip(locs, ihHeights)):
                e = expected_hB(scenario, loc)
                summary["interfaceHeight"][i] = {"hB": float(hB), "hL": float(hL), "exact": e}
                if ihErrs[i] is not None and not ihErrs[i] <= 1e-6:
                    differing.append(f"{i}:hB {hB:.3g}" + ("" if e is None else f" (exact {e:.3g})"))
            print(f"  {'interfaceHeight':28s} lines differing: {', '.join(differing) if differing else 'none'}")

with open("results.json", "w") as f:
    json.dump({"passed": allOk, "tol": TOL, "runs": results}, f, indent=1)

print("\nALL SURFACEPOSITION CHECKS PASSED" if allOk else "\nSOME SURFACEPOSITION CHECKS FAILED")
sys.exit(0 if allOk else 1)
