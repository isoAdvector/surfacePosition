#!/usr/bin/env python3
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Summarises the compatibility checks of Allrun: per OpenFOAM version, whether
# it is built and runs, whether surfacePosition compiles (with the distinct
# compiler errors), and whether the specialCases and movingMesh tests pass
# (with the first fatal error of a failing run). Compile errors are also
# listed for a second build with fileFormats on the include path, where
# Allrun made one. Writes results.json.
#
# Usage: summarise.py <workDir>

import glob
import json
import os
import re
import sys

work = sys.argv[1]


def read(path, default=None):
    try:
        return open(path).read().strip()
    except OSError:
        return default


def compile_errors(log):
    """Distinct compiler errors, without file positions and templates"""
    out = []
    for line in open(log, errors="replace"):
        m = re.search(r"(?:fatal )?error: (.*)", line)
        if not m:
            continue
        msg = m.group(1)
        msg = re.sub(r"‘class Foam::DynamicList<[^’]*’", "DynamicList", msg)
        msg = re.sub(r"<[^<>]*>", "", msg)
        msg = re.sub(r"<[^<>]*>", "", msg)
        if msg not in out:
            out.append(msg)
    return out


def fatal_error(testDir):
    """The first fatal error in the logs of a test: the first line of a FOAM
    FATAL ERROR, or a crash with the innermost OpenFOAM function of its stack
    trace and whether surfacePosition is in the trace"""
    for log in sorted(glob.glob(os.path.join(testDir, "*", "log.*"))):
        case = os.path.basename(os.path.dirname(log))
        lines = open(log, errors="replace").read().splitlines()
        for i, line in enumerate(lines):
            if "FOAM FATAL" in line:
                for text in lines[i + 1:]:
                    if text.strip():
                        return f"{case}: {text.strip()}"
            m = re.search(r"Foam::(sigFpe|sigSegv)::sigHandler", line)
            if m:
                kind = "floating point exception" if m.group(1) == "sigFpe" else "segmentation fault"
                where = ""
                for text in lines[i + 1:i + 6]:
                    f = re.search(r"#\d+\s+(Foam::[A-Za-z_:]+)", text)
                    if f:
                        where = f" in {f.group(1)}"
                        break
                trace = "\n".join(lines[i + 1:i + 30])
                origin = "surfacePosition" if "surfacePosition" in trace else "solver"
                return f"{case}: {kind}{where} ({origin})"
    return None


results = {}
for w in sorted(glob.glob(os.path.join(work, "v*"))):
    v = os.path.basename(w)
    r = results[v] = {
        "install": read(os.path.join(w, "install")),
        "state": read(os.path.join(w, "state"), "unknown"),
        "compiled": None,
        "errors": [],
        "tests": {},
    }
    status = read(os.path.join(w, "compiled"))
    if status is not None:
        r["compiled"] = (status == "0")
        if not r["compiled"]:
            r["errors"] = compile_errors(os.path.join(w, "log.wmake"))
            log2 = os.path.join(w, "log.wmake2")
            if os.path.isfile(log2):
                r["errorsWithFileFormats"] = compile_errors(log2)
    for t in ("specialCases", "movingMesh"):
        s = read(os.path.join(w, t, "status"))
        if s is not None:
            error = None
            if s == "124" or s == "137":
                error = "time limit reached"
            elif s != "0":
                error = fatal_error(os.path.join(w, t))
            without = read(os.path.join(w, t, "withoutFO", "result"))
            r["tests"][t] = {
                "passed": s == "0", "error": error,
                "withoutFunctionObjects": without.splitlines() if without else [],
            }

with open("results.json", "w") as f:
    json.dump(results, f, indent=1)

print(f"{'version':8s} {'installation':14s} {'compiles':9s} {'specialCases':13s} {'movingMesh':11s}")
for v, r in results.items():
    comp = "-" if r["compiled"] is None else ("yes" if r["compiled"] else "no")
    tests = [("passed" if r["tests"][t]["passed"] else "FAILED") if t in r["tests"] else "-"
             for t in ("specialCases", "movingMesh")]
    print(f"{v:8s} {r['state']:14s} {comp:9s} {tests[0]:13s} {tests[1]:11s}")
    for e in r["errors"]:
        print(f"         compile error: {e}")
    for e in r.get("errorsWithFileFormats", []):
        print(f"         with fileFormats included: {e}")
    for t, s in r["tests"].items():
        if s["error"]:
            print(f"         {t}: {s['error']}")
        for line in s["withoutFunctionObjects"]:
            print(f"         {t} without function objects: {line}")
