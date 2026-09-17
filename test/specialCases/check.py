#!/usr/bin/env python3
#
# Compares the points reported by surfacePosition (and interfaceHeight) in the
# special cases with the exact interface positions, serial and parallel.
#
# Usage: check.py <scenario> ...   (cases <scenario> and <scenario>_parallel)

import json
import math
import os
import sys

import numpy as np

NAN = float("nan")
TOL = 1e-8      # tolerance for the geometric methods [m]
# The alpha method interpolates linearly between cell centres, so it is only
# expected within half a cell (0.05 m) for an interface inside a cell


def expected(scenario):
    """Exact (x, z) per line, or None where no interface should be found."""
    if scenario.startswith("flat"):
        z = int(scenario[4:])/1000
        vertical = [(None, z)]*12
        vertical[10] = None
        lines = [(0.33, z), (0.33, z), None, None, None, (0.25, z)]
        return vertical, lines
    # Slab of water for x < 0.5, cyclic in x
    vertical = [None]*12
    lines = [None, None, (0.0, None), (0.5, None), (1.0, None), None]
    return vertical, lines


def last_row(case, fo, fname):
    p = os.path.join(case, "postProcessing", fo, "0", fname)
    if not os.path.isfile(p):
        return None
    d = np.atleast_2d(np.genfromtxt(p, comments="#"))
    return d[-1]


def points(row, nLines, cols):
    """(x, z) per line from a surfacePosition row with cols per line."""
    out = []
    for i in range(nLines):
        x, z = row[1 + i*cols + cols - 3], row[1 + i*cols + cols - 1]
        out.append((x, z))
    return out


def verdict(got, exp, tol):
    gx, gz = got
    if exp is None:
        return ("ok", 0.0) if (math.isnan(gx) and math.isnan(gz)) else ("FOUND", NAN)
    if math.isnan(gz):
        return ("NAN", NAN)
    err = max(abs(gx - exp[0]) if exp[0] is not None else 0.0,
              abs(gz - exp[1]) if exp[1] is not None else 0.0)
    return ("ok" if err <= tol else "WRONG", err)


allOk = True
results = {}
for scenario in sys.argv[1:]:
    vExp, lExp = expected(scenario)
    print(f"\n=== {scenario}")
    for run in (scenario, scenario + "_parallel"):
        if not os.path.isdir(run):
            continue
        print(f"--- {run}")
        summary = results[run] = {"checks": {}, "interfaceHeight": {}}
        rows = []
        for m in ("reconstructed", "plic", "alpha"):
            rows.append((f"vertical {m}", last_row(run, f"vertical_{m}", "surfacePosition.dat"), 5, vExp,
                         TOL if m != "alpha" else 0.05, True))
        for m in ("reconstructed", "plic"):
            rows.append((f"lines {m}", last_row(run, f"lines_{m}", "surfacePosition.dat"), 4, lExp, TOL, True))
        for label, row, cols, exp, tol, strict in rows:
            if row is None:
                print(f"  {label:24s} no output")
                allOk = False
                continue
            res = [verdict(p, e, tol) for p, e in zip(points(row, len(exp), cols), exp)]
            bad = [i for i, (v, _) in enumerate(res) if v != "ok"]
            maxErr = max((e for v, e in res if v == "ok"), default=0.0)
            status = "ok" if not bad else "FAILED lines " + ", ".join(f"{i} ({res[i][0]})" for i in bad)
            if bad and strict:
                allOk = False
            summary["checks"][label] = {
                "maxError": maxErr, "nLines": len(exp), "failed": bad,
                "noInterface": sum(e is None for e in exp),
            }
            print(f"  {label:24s} max error {maxErr:.2e}  {status}")
        # interfaceHeight: hB is the height above the lower boundary z = 0
        row = last_row(run, "vertical_interfaceHeight", "height.dat")
        if row is not None:
            errs = []
            for i, e in enumerate(vExp):
                hB = row[1 + 2*i]
                summary["interfaceHeight"][i] = {
                    "hB": float(hB), "exact": None if e is None else e[1]
                }
                if e is None:
                    errs.append(f"{i}:{hB:.3g}")
                else:
                    err = abs(hB - e[1])
                    if err > 1e-6:
                        errs.append(f"{i}:{hB - e[1]:+.3g}")
            print(f"  {'vertical interfaceHeight':24s} lines differing: {', '.join(errs) if errs else 'none'}")

with open("results.json", "w") as f:
    json.dump({"passed": allOk, "tol": TOL, "runs": results}, f, indent=1)

print("\nALL SURFACEPOSITION CHECKS PASSED" if allOk else "\nSOME SURFACEPOSITION CHECKS FAILED")
sys.exit(0 if allOk else 1)
