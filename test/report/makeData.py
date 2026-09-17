#!/usr/bin/env python3
#
# Writes the numbers, tables and plot data of the report from the results of
# the executed test suite (results.json of each test and the postProcessing
# data of its cases) into generated/:
#   results.tex        macros with single numbers
#   tab*.tex           tables
#   wave*.csv          time series for the plots
#   results.json       the macros and table rows in Markdown, for the README
#
# Usage: makeData.py   (from the report directory, after ../Allrun)

import json
import os
import re
import subprocess
import time

import numpy as np

TESTS = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "generated")
os.makedirs(OUT, exist_ok=True)

macros = {}
tables = {}


def load(test):
    with open(os.path.join(TESTS, test, "results.json")) as f:
        return json.load(f)


def macro(name, value):
    assert re.fullmatch(r"[A-Za-z]+", name), name
    macros[name] = value


def num(x, digits=3):
    """Number for siunitx \\num"""
    return f"\\num{{{x:.{digits}g}}}" if x != 0 else "\\num{0}"


def fixed(x, decimals):
    return f"{x:.{decimals}f}"


def write(name, text):
    with open(os.path.join(OUT, name), "w") as f:
        f.write(text)


def table(name, rows, nColumns):
    """Table rows: lists of cells, or a string for a group heading"""
    lines = []
    for row in rows:
        if isinstance(row, str):
            if lines:
                lines.append(r"\addlinespace")
            lines.append(rf"\multicolumn{{{nColumns}}}{{l}}{{{row}}} \\")
        else:
            lines.append(" & ".join(row) + r" \\")
    write(f"{name}.tex", "\n".join(lines) + "\n")
    tables[name] = [
        "| " + (tomd(row) + " |"*nColumns if isinstance(row, str)
                else " | ".join(tomd(c) for c in row) + " |")
        for row in rows
    ]


SUPERSCRIPT = str.maketrans("-0123456789", "\u207b\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079")
UNITS = {"m": " m", "mm": " mm", "s": " s", "ms": " ms", "\\%": " %"}


def tomd(text):
    """Markdown version of the LaTeX of a macro or table cell"""
    def number(x):
        m = re.fullmatch(r"([-\d.]+)e([-+]?)0*(\d+)", x)
        if not m:
            return x
        return f"{m.group(1)}\u00b710" + (m.group(2).replace("+", "") + m.group(3)).translate(SUPERSCRIPT)

    text = re.sub(r"\\SI\{([^{}]*)\}\{([^{}]*)\}", lambda m: number(m.group(1)) + UNITS[m.group(2)], text)
    text = re.sub(r"\\num\{([^{}]*)\}", lambda m: number(m.group(1)), text)
    text = re.sub(r"\\texttt\{([^{}]*)\}", r"`\1`", text)
    text = re.sub(r"\\textbf\{([^{}]*)\}", r"**\1**", text)
    text = re.sub(r"\\emph\{([^{}]*)\}", r"*\1*", text)
    text = re.sub(r"\$([^$]*)\$", lambda m: m.group(1).replace("\\times", "\u00d7"), text)
    return text.replace("\\,", "\u202f").replace("--", "\u2013")


# --------------------------------------------------------------------------- #
# waveAccuracy

wave = load("waveAccuracy")
METHODS = [
    ("interfaceHeight cellPoint", r"\texttt{interfaceHeight}, cellPoint"),
    ("interfaceHeight cell", r"\texttt{interfaceHeight}, cell"),
    ("surfacePosition plic", r"\texttt{surfacePosition}, plic"),
    ("surfacePosition alpha", r"\texttt{surfacePosition}, alpha"),
]
MESH = {"coarse": (38, 8), "medium": (95, 20)}

rows = []
for case in ("coarse", "medium"):
    r = wave[case]
    nx, nz = MESH[case]
    rows.append(rf"\emph{{{case}}}, ${nx}\times{nz}$ cells")
    for key, label in METHODS:
        st = r["methods"][key]
        rows.append([label] + [fixed(st[o]["rms"], 3) for o in ("0.50", "0.25", "0.02")]
                    + [fixed(st["all"]["mean"], 3), fixed(st["all"]["max"], 3)])
table("tabWave", rows, 6)

for case, cap in (("coarse", "Coarse"), ("medium", "Medium")):
    r = wave[case]
    macro(f"wave{cap}Samples", str(r["samples"]))
    macro(f"wave{cap}Lines", str(r["lines"]))
    macro(f"wave{cap}PlicRms", fixed(r["methods"]["surfacePosition plic"]["all"]["rms"], 3))
    macro(f"wave{cap}CellPointRms", fixed(r["methods"]["interfaceHeight cellPoint"]["all"]["rms"], 3))
    macro(f"wave{cap}CellMean", fixed(r["methods"]["interfaceHeight cell"]["all"]["mean"], 1))
    macro(f"wave{cap}AlphaRms", fixed(r["methods"]["surfacePosition alpha"]["all"]["rms"], 2))
    macro(f"wave{cap}ParallelMax", f"{max(r['parallel'].values()):.3g}")
    macro(f"wave{cap}Nan", str(sum(st["all"]["nan"] for st in r["methods"].values())))

pmax = max(max(wave[c]["parallel"].values()) for c in ("coarse", "medium"))
macro("waveParallel", "are identical" if pmax == 0 else f"differ by up to \\SI{{{pmax:.3g}}}{{m}}")

# Time series of one line per offset in the medium case, for the plots
case = os.path.join(TESTS, "waveAccuracy", "medium")
ex = np.genfromtxt(os.path.join(case, "postProcessing", "exactSurface", "exactSurface.dat"), comments="#")
n = wave["medium"]["lines"]
zb = float(re.search(r"^\s*z1\s+([-\d.eE+]+)\s*;",
                     open(os.path.join(case, "system", "blockMeshDict")).read(), re.M).group(1))


def series(fo, fname):
    d = np.genfromtxt(os.path.join(case, "postProcessing", fo, "0", fname), comments="#")
    if fname == "height.dat":
        return d[:, 0], d[:, 1:1 + 2*n:2] + zb
    return d[:, 0], d[:, 5:1 + 5*n:5]


t, ihp = series("interfaceHeightCellPoint", "height.dat")
_, ihc = series("interfaceHeightCell", "height.dat")
_, spp = series("surfacePositionPlic", "surfacePosition.dat")
_, spa = series("surfacePositionAlpha", "surfacePosition.dat")
exact = {round(tt, 9): row for tt, row in zip(ex[:, 0], ex[:, 1:])}
keep = [i for i, tt in enumerate(t) if round(tt, 9) in exact]
for line in (6, 8):
    with open(os.path.join(OUT, f"waveLine{line}.csv"), "w") as f:
        f.write("t exact ihCellPoint ihCell spPlic spAlpha\n")
        for i in keep:
            e = exact[round(t[i], 9)][line]
            f.write(f"{t[i]:.4f} {e*1e3:.4f} {(ihp[i, line] - e)*1e3:.4f} {(ihc[i, line] - e)*1e3:.4f}"
                    f" {(spp[i, line] - e)*1e3:.4f} {(spa[i, line] - e)*1e3:.4f}\n")
info = re.findall(r"^//\s+(\d+):\s+(\d+)\s+([\d.]+)",
                  open(os.path.join(case, "system", "lineLocations")).read(), re.M)
for line, name in ((6, "A"), (8, "B")):
    macro(f"waveLine{name}Offset", info[line][2])
    macro(f"waveLine{name}", str(line))

# --------------------------------------------------------------------------- #
# specialCases

special = load("specialCases")
CHECKS = ["vertical reconstructed", "vertical plic", "vertical alpha",
          "lines reconstructed", "lines plic"]
rows = []
for scenario, label in (("flat537", r"flat, $z=0.537$"), ("flat500", r"flat, $z=0.5$"),
                        ("flat300", r"flat, $z=0.3$"), ("slab", "vertical slab")):
    for run, procs in ((scenario, "1"), (scenario + "_parallel", "8")):
        r = special["runs"][run]["checks"]
        cells = []
        for c in CHECKS:
            st = r[c]
            err = num(st["maxError"], 2) if st["maxError"] > 0 else "0"
            cells.append(err if not st["failed"] else r"\textbf{failed}")
        rows.append([label if procs == "1" else "", procs] + cells)
table("tabSpecial", rows, 7)
macro("specialAlphaError", fixed(special["runs"]["flat537"]["checks"]["vertical alpha"]["maxError"]*1e3, 1))
macro("specialPassed", "passed" if special["passed"] else "FAILED")
macro("specialTol", num(special["tol"]))

# interfaceHeight on the special vertical lines of the flat interface at z = 0.5
rows = []
names = {3: "processor face", 4: "processor edge", 7: "cyclic/processor edge",
         8: "far above", 9: "far below", 10: "outside the domain"}
for i, name in names.items():
    s = special["runs"]["flat500"]["interfaceHeight"][str(i)]["hB"]
    p = special["runs"]["flat500_parallel"]["interfaceHeight"][str(i)]["hB"]
    exact = special["runs"]["flat500"]["interfaceHeight"][str(i)]["exact"]
    ex = "--" if exact is None else fixed(exact, 3)
    rows.append([str(i), name, ex, fixed(s, 3), fixed(p, 3)])
table("tabSpecialIH", rows, 5)

# --------------------------------------------------------------------------- #
# movingMesh

moving = load("movingMesh")
rows = []
for run, label in (("amr", "amr"), ("amr_parallel", "amr"), ("motion", "motion"),
                   ("motion_parallel", "motion")):
    r = moving["runs"][run]
    procs = "4" if run.endswith("parallel") else "1"
    rb = r["rebuilds"]["surfacePosition1"]
    rbp = r["rebuilds"]["surfacePositionPlic1"]
    assert rb == rbp, (rb, rbp)
    el = r["elevation"]
    rows.append([
        label, procs, str(r["steps"]), str(rb["rebuilds"]), str(rb["cells"]),
        str(rb["missing"]), str(rb["extra"]),
        str(el["surfacePosition1"]["nan"] + el["surfacePositionPlic1"]["nan"]),
        fixed(el["surfacePosition1"]["rms"], 2), fixed(el["surfacePositionPlic1"]["rms"], 2),
    ])
table("tabMoving", rows, 10)
macro("movingPassed", "passed" if moving["passed"] else "FAILED")
macro("amrRebuilds", str(moving["runs"]["amr"]["rebuilds"]["surfacePosition1"]["rebuilds"]))
macro("amrParallelRebuilds", str(moving["runs"]["amr_parallel"]["rebuilds"]["surfacePosition1"]["rebuilds"]))
macro("amrSteps", str(moving["runs"]["amr"]["steps"]))
macro("amrMaxDiff", fixed(moving["runs"]["amr"]["elevation"]["surfacePosition1"]["max"], 0))
macro("amrParallelDiffIH", fixed(moving["runs"]["amr parallel"]["interfaceHeight1"]*1e3, 1))
macro("motionRebuilds", str(moving["runs"]["motion"]["rebuilds"]["surfacePosition1"]["rebuilds"]))
macro("motionCells", str(moving["runs"]["motion_parallel"]["rebuilds"]["surfacePosition1"]["cells"]))
macro("motionParallelDiffSP", fixed(moving["runs"]["motion parallel"]["surfacePosition1"]*1e3, 2))
macro("motionParallelDiffIH", fixed(moving["runs"]["motion parallel"]["interfaceHeight1"]*1e3, 2))

# Time step of the amr run, from its log
log = open(os.path.join(TESTS, "movingMesh", "amr", "log.interIsoFoam")).read()
dts = [float(x) for x in re.findall(r"^deltaT = (\S+)", log, re.M)]
macro("amrDeltaT", f"{float(np.median(dts)):.2g}")
cells = [int(x) for x in re.findall(r"(?:Refined|Unrefined) from \d+ to (\d+) cells", log)]
macro("amrCellsMax", str(max(cells)))
macro("amrEndTime", fixed(float(re.findall(r"^Time = (\S+)", log, re.M)[-1]), 4))

# --------------------------------------------------------------------------- #
# largeTiming

timing = load("largeTiming")
cost = timing["cost"]
rows = []
solver = re.search(r"^application\s+(\w+);", open(os.path.join(
    TESTS, "largeTiming", "case", "system", "controlDict")).read(), re.M).group(1)
macro("timingSolver", solver)
for key, label in (("flow solution", f"flow solution (\\texttt{{{solver}}})"),
                   ("surfacePosition1", r"\texttt{surfacePosition}, reconstructed"),
                   ("surfacePositionPlic1", r"\texttt{surfacePosition}, plic"),
                   ("interfaceHeight1", r"\texttt{interfaceHeight}")):
    c = cost[key]
    perLine = "" if key == "flow solution" else fixed(c["msPerStep"]/timing["lines"], 3)
    startup = fixed(c["startup"], 0) if key.startswith("surfacePosition") else ""
    rows.append([label, fixed(c["msPerStep"], 1), perLine, fixed(c["spread"], 2), startup])
table("tabTiming", rows, 5)
macro("timingCells", f"{timing['cells']:,}".replace(",", "\\,"))
macro("timingProcs", str(timing["processors"]))
macro("timingLines", str(timing["lines"]))
macro("timingSteps", str(timing["steps"]))
macro("timingSP", fixed(cost["surfacePosition1"]["msPerStep"], 1))
macro("timingSPPlic", fixed(cost["surfacePositionPlic1"]["msPerStep"], 1))
macro("timingIH", fixed(cost["interfaceHeight1"]["msPerStep"], 0))
macro("timingFlow", fixed(cost["flow solution"]["msPerStep"], 0))
macro("timingRatio", fixed(cost["interfaceHeight1"]["msPerStep"]/cost["surfacePosition1"]["msPerStep"], 0))
macro("timingIHFraction", fixed(100*cost["interfaceHeight1"]["msPerStep"]/cost["flow solution"]["msPerStep"], 0))
macro("timingSPFraction", fixed(100*cost["surfacePosition1"]["msPerStep"]/cost["flow solution"]["msPerStep"], 1))
macro("timingStartup", fixed(cost["surfacePosition1"]["startup"], 0))
macro("timingDiffRms", fixed(timing["elevation"]["surfacePosition1"]["rms"], 3))
macro("timingDiffMax", fixed(timing["elevation"]["surfacePosition1"]["max"], 3))
macro("timingDiffPlicMax", fixed(timing["elevation"]["surfacePositionPlic1"]["max"], 3))

# Wall time of the whole suite, from ../log.Allrun
m = re.search(r"^real\s+(\d+)m([\d.]+)s", open(os.path.join(TESTS, "log.Allrun")).read(), re.M)
macro("suiteTime", f"{int(m.group(1))}\\,min\\,{float(m.group(2)):.0f}\\,s" if m else "--")

# The run: time of the suite log, and the last revision of the source and the
# tests (the report and the README are left out, as they are generated after
# the run and committed with its results)
stamp = time.strftime("%Y-%m-%d %H:%M", time.localtime(os.path.getmtime(os.path.join(TESTS, "log.Allrun"))))
tested = [".", ":(exclude)test/report", ":(exclude)README.md"]


def git(*args):
    return subprocess.run(["git", "-C", os.path.join(TESTS, ".."), *args],
                          capture_output=True, text=True).stdout.strip()


rev = git("log", "-1", "--format=%h", "--", *tested)
if git("status", "--porcelain", "--untracked-files=no", "--", *tested):
    rev += " with uncommitted changes"
macro("testRun", f"{stamp} (source and tests at revision \\texttt{{{rev}}})")

write("results.tex", "".join(f"\\newcommand{{\\{k}}}{{{v}}}\n" for k, v in sorted(macros.items())))
write("results.json", json.dumps(
    {"macros": {k: tomd(v) for k, v in macros.items()}, "tables": tables},
    indent=1, ensure_ascii=False))
print(f"Wrote {len(macros)} macros and the tables to {OUT}")
