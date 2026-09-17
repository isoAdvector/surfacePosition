#!/usr/bin/env pvbatch
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
#
# ParaView figures of the test cases of the executed test suite, written to
# figures/. Run with pvbatch from the report directory, after ../Allrun.
# The processor decompositions are the cellDist fields written by
# decomposePar -cellDist in the undecomposed directories of the parallel
# cases.

import os
import re
import shutil

import numpy as np
from paraview.simple import (
    AppendDatasets, ColorBy, Contour, CreateRenderView, Delete, ExtractSurface,
    GetColorTransferFunction, GetOpacityTransferFunction, Glyph, Hide, Line, ResetCamera,
    OpenFOAMReader, PointSource, Render, SaveScreenshot, Show, Slice,
    Threshold, Tube,
)

HERE = os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.abspath(os.path.join(HERE, ".."))
FIGURES = os.path.join(HERE, "figures")
GENERATED = os.path.join(HERE, "generated")
os.makedirs(FIGURES, exist_ok=True)

WATER = [0.25, 0.55, 0.85]
AIR = [0.96, 0.96, 0.96]
LINE = [0.85, 0.1, 0.1]
POINT = [0.0, 0.0, 0.0]
# Categorical colours for the processors
PROCS = [
    [0.12, 0.47, 0.71], [0.68, 0.78, 0.91], [1.00, 0.50, 0.05], [1.00, 0.73, 0.47],
    [0.17, 0.63, 0.17], [0.60, 0.87, 0.54], [0.84, 0.15, 0.16], [1.00, 0.60, 0.59],
    [0.58, 0.40, 0.74], [0.77, 0.69, 0.84], [0.55, 0.34, 0.29], [0.77, 0.61, 0.58],
    [0.89, 0.47, 0.76], [0.97, 0.71, 0.82], [0.50, 0.50, 0.50], [0.78, 0.78, 0.78],
    [0.74, 0.74, 0.13], [0.86, 0.86, 0.55], [0.09, 0.75, 0.81], [0.62, 0.85, 0.90],
]


def case(*path):
    return os.path.join(TESTS, *path)


def reader(caseDir, fields, time=None):
    r = OpenFOAMReader(FileName=os.path.join(caseDir, "system", "controlDict"))
    # The time directories are listed again with the zero time included
    r.SkipZeroTime = 0
    r.Refresh()
    r.UpdatePipelineInformation()
    r.MeshRegions = ["internalMesh"]
    r.CellArrays = fields
    r.Createcelltopointfiltereddata = 1
    times = list(r.TimestepValues) if r.TimestepValues else [0.0]
    t = times[-1] if time is None else min(times, key=lambda x: abs(x - time))
    r.UpdatePipeline(t)
    return r, t


def view(size, time=0.0):
    v = CreateRenderView()
    v.ViewSize = size
    v.OrientationAxesVisibility = 0
    v.Background = [1, 1, 1]
    v.UseColorPaletteForBackground = 0
    v.ViewTime = time
    v.AxesGrid.Visibility = 0
    return v


def show(src, v):
    d = Show(src, v)
    d.DataAxesGrid.GridAxesVisibility = 0
    return d


def alpha_lut(name, lo, hi):
    """From air (alpha = 0) to water (alpha = 1)"""
    lut = GetColorTransferFunction(name)
    lut.RGBPoints = [0.0] + AIR + [1.0] + WATER
    lut.RescaleTransferFunction(lo, hi)
    lut.AutomaticRescaleRangeMode = "Never"
    return lut


def proc_lut(name, n):
    lut = GetColorTransferFunction(name)
    lut.InterpretValuesAsCategories = 1
    lut.Annotations = sum(([str(i), str(i)] for i in range(n)), [])
    lut.IndexedColors = sum((PROCS[i % len(PROCS)] for i in range(n)), [])
    return lut


def show_lines(v, segments, colour=LINE, radius=0.004):
    for p0, p1 in segments:
        src = Line(Point1=list(p0), Point2=list(p1))
        tube = Tube(Input=src, Radius=radius, NumberofSides=12)
        d = show(tube, v)
        d.ColorArrayName = ["POINTS", ""]
        d.AmbientColor = colour
        d.DiffuseColor = colour


def show_points(v, points, radius):
    for p in points:
        if np.any(np.isnan(p)):
            continue
        src = PointSource(Center=list(p), NumberOfPoints=1, Radius=0)
        g = Glyph(Input=src, GlyphType="Sphere")
        g.ScaleArray = ["POINTS", "No scale array"]
        g.ScaleFactor = 2*radius
        g.GlyphMode = "All Points"
        d = show(g, v)
        d.ColorArrayName = ["POINTS", ""]
        d.AmbientColor = POINT
        d.DiffuseColor = POINT


def camera(v, position, focal, up, scale=None):
    Render(v)
    v.CameraPosition = position
    v.CameraFocalPoint = focal
    v.CameraViewUp = up
    if scale is not None:
        v.CameraParallelProjection = 1
        v.CameraParallelScale = scale


def last_positions(case, fo):
    """Interface positions of the last write (position.dat, as written by
    interfaceHeight); None for a line missing the mesh"""
    path = os.path.join(case, "postProcessing", fo, "0", "position.dat")
    last = [l for l in open(path) if not l.startswith("#")][-1]
    out = []
    for m in re.findall(r"\(([^)]*)\)", last):
        p = tuple(float(x) for x in m.split())
        out.append(None if max(abs(x) for x in p) > 1e10 else p)
    return out


def locations(path):
    text = open(path).read()
    block = re.search(r"lineLocations\s*\((.*?)\n\);", text, re.S).group(1)
    return [tuple(float(x) for x in m.split())
            for m in re.findall(r"\(([^()]*)\)", block)]


def save(v, name, size):
    Render(v)
    SaveScreenshot(os.path.join(FIGURES, name), v, ImageResolution=size,
                   TransparentBackground=0)
    print("wrote", name)


# --------------------------------------------------------------------------- #
# waveAccuracy: medium case at the end of the period, with the lines and the
# positions found by surfacePosition (plic)

def figure_wave():
    c = case("waveAccuracy", "medium")
    r, t = reader(c, ["alpha.water"])
    size = [2400, 620]
    v = view(size, t)
    d = show(r, v)
    ColorBy(d, ("CELLS", "alpha.water"))
    d.LookupTable = alpha_lut("alphaWave", 0, 1)
    d.SetRepresentationType("Surface With Edges")
    d.EdgeColor = [0.6, 0.6, 0.6]

    locs = locations(os.path.join(c, "system", "lineLocations"))
    show_lines(v, [((x, -0.11, -0.4), (x, -0.11, 0.4)) for x, _, _ in locs], radius=0.003)
    pts = [(p[0], -0.12, p[2]) for p in last_positions(c, "surfacePosition") if p]
    show_points(v, pts, 0.008)

    camera(v, [1.894, -10, 0], [1.894, 0, 0], [0, 0, 1], scale=0.5)
    save(v, "wave.png", size)
    Delete(v)


# --------------------------------------------------------------------------- #
# specialCases: 3D views of the flat interface at z = 0.5 and the slab, on 8
# processors, with the lines and the points found

def special_locations(c):
    text = open(os.path.join(c, "system", "controlDict")).read()
    block = text[text.index("\n    vertical\n"):]
    locs = re.search(r"locations\s*\((.*?)\n\s*\);", block, re.S).group(1)
    return [tuple(float(x) for x in m.split()) for m in re.findall(r"\(([^()]*)\)", locs)]


def figure_special(name):
    c = case("specialCases", name + "_parallel")
    r, t = reader(c, ["alpha.water", "cellDist"], 0.0)
    size = [1400, 1300]
    v = view(size, t)

    # Water, coloured by processor
    water = Threshold(Input=r, Scalars=["CELLS", "alpha.water"])
    water.LowerThreshold = 0.5
    water.UpperThreshold = 2
    d = show(water, v)
    ColorBy(d, ("CELLS", "cellDist"))
    d.LookupTable = proc_lut("procSpecial", 8)
    d.SetRepresentationType("Surface With Edges")
    d.EdgeColor = [0.3, 0.3, 0.3]
    d.Opacity = 0.55

    # Outline of the domain
    o = show(r, v)
    o.SetRepresentationType("Outline")
    o.ColorArrayName = ["POINTS", ""]
    o.AmbientColor = [0, 0, 0]
    o.DiffuseColor = [0, 0, 0]

    vertical = special_locations(c)
    show_lines(v, [((x, y, -0.35), (x, y, 1.35)) for x, y, _ in vertical], radius=0.008)
    show_points(v, [p for p in last_positions(c, "vertical") if p], 0.02)

    camera(v, [3.6, -3.4, 2.6], [0.95, 0.5, 0.5], [0, 0, 1])
    v.CameraViewAngle = 30
    save(v, f"special_{name}.png", size)
    Delete(v)


# --------------------------------------------------------------------------- #
# movingMesh: refined mesh at the end of the amr run, and the moving mesh of
# the motion case coloured by processor at two times

def figure_amr():
    c = case("movingMesh", "amr")
    r, t = reader(c, ["alpha.water"])
    size = [1300, 1300]
    v = view(size, t)
    d = show(r, v)
    ColorBy(d, ("CELLS", "alpha.water"))
    d.LookupTable = alpha_lut("alphaAmr", 0, 1)
    d.SetRepresentationType("Surface With Edges")
    d.EdgeColor = [0.35, 0.35, 0.35]

    locs = locations(os.path.join(c, "system", "lineLocations"))
    show_lines(v, [((x, -0.011, 0), (x, -0.011, 1)) for x, _, _ in locs], radius=0.003)
    show_points(v, [(p[0], -0.012, p[2]) for p in last_positions(c, "surfacePosition1") if p], 0.01)

    camera(v, [0.5, -10, 0.5], [0.5, 0, 0.5], [0, 0, 1], scale=0.52)
    save(v, "amr.png", size)
    Delete(v)


def motion_view_case():
    """Serial motion case with the cellDist of the parallel run"""
    src = case("movingMesh", "motion")
    dist = case("movingMesh", "motion_parallel", "0", "cellDist")
    dst = os.path.join(GENERATED, "motionView")
    shutil.rmtree(dst, ignore_errors=True)
    os.makedirs(dst)
    for d in ("constant", "system"):
        os.symlink(os.path.join(src, d), os.path.join(dst, d))
    times = sorted((d for d in os.listdir(src) if re.fullmatch(r"[\d.]+", d)), key=float)
    for tm in times:
        os.makedirs(os.path.join(dst, tm))
        for f in os.listdir(os.path.join(src, tm)):
            p = os.path.join(src, tm, f)
            if f in ("alpha.water", "polyMesh"):
                os.symlink(p, os.path.join(dst, tm, f))
        shutil.copy(dist, os.path.join(dst, tm, "cellDist"))
    return dst


def figure_motion(time, tag, highlight):
    c = motion_view_case() if not os.path.isdir(os.path.join(GENERATED, "motionView")) \
        else os.path.join(GENERATED, "motionView")
    r, t = reader(c, ["alpha.water", "cellDist"], time)
    size = [1100, 1100]
    v = view(size, t)
    d = show(r, v)
    ColorBy(d, ("CELLS", "cellDist"))
    d.LookupTable = proc_lut("procMotion", 4)
    d.SetRepresentationType("Surface With Edges")
    d.EdgeColor = [0.3, 0.3, 0.3]

    # Interface as a contour of alpha
    s = Slice(Input=r)
    s.SliceType.Normal = [0, 1, 0]
    s.SliceType.Origin = [0, 0, 0]
    cont = Contour(Input=s, ContourBy=["POINTS", "alpha.water"], Isosurfaces=[0.5])
    dc = show(Tube(Input=cont, Radius=0.004), v)
    dc.ColorArrayName = ["POINTS", ""]
    dc.AmbientColor = WATER
    dc.DiffuseColor = WATER

    locs = locations(os.path.join(c, "system", "lineLocations"))
    others = [((x, -0.011, -0.05), (x, -0.011, 1.05)) for i, (x, _, _) in enumerate(locs) if i != highlight]
    show_lines(v, others, colour=[0.25, 0.25, 0.25], radius=0.003)
    x = locs[highlight][0]
    show_lines(v, [((x, -0.012, -0.05), (x, -0.012, 1.05))], radius=0.004)

    camera(v, [0.5, -10, 0.5], [0.5, 0, 0.5], [0, 0, 1], scale=0.56)
    save(v, f"motion_{tag}.png", size)
    Delete(v)


# --------------------------------------------------------------------------- #
# largeTiming: the cube on 18 processors, water coloured by processor, with
# the 50 lines

def figure_timing():
    c = case("largeTiming", "case")
    r, t = reader(c, ["alpha.water", "cellDist"], 0.0)
    nProcs = len([d for d in os.listdir(c) if d.startswith("processor")])
    size = [1400, 1300]
    v = view(size, t)

    water = Threshold(Input=r, Scalars=["CELLS", "alpha.water"])
    water.LowerThreshold = 0.5
    water.UpperThreshold = 2
    surf = ExtractSurface(Input=water)
    d = show(surf, v)
    ColorBy(d, ("CELLS", "cellDist"))
    d.LookupTable = proc_lut("procTiming", nProcs)

    o = show(r, v)
    o.SetRepresentationType("Outline")
    o.ColorArrayName = ["POINTS", ""]
    o.AmbientColor = [0, 0, 0]
    o.DiffuseColor = [0, 0, 0]

    locs = locations(os.path.join(c, "system", "lineLocations"))
    show_lines(v, [((x, y, 0), (x, y, 1.25)) for x, y, _ in locs], radius=0.004)

    camera(v, [2.9, -2.4, 2.4], [0.5, 0.5, 0.45], [0, 0, 1])
    v.CameraViewAngle = 30
    save(v, "timing.png", size)
    Delete(v)


shutil.rmtree(os.path.join(GENERATED, "motionView"), ignore_errors=True)
figure_wave()
figure_special("flat500")
figure_special("slab")
figure_amr()
figure_motion(0.0, "t0", 3)
figure_motion(0.2, "t2", 3)
figure_timing()
