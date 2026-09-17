# surfacePosition

An OpenFOAM (v2512) function object that reports where a volume-of-fluid
interface crosses user-specified lines, as an alternative to
`interfaceHeight`. The cells crossed by each line are found once and followed
through mesh changes (adaptive refinement, mesh motion, also across processor
boundaries) by a flood fill, so the cost per write does not depend on the mesh
size. The crossing is computed geometrically in the interface cell, from the
solver's reconstructed interface plane or from a PLIC plane.

It links only OpenFOAM's `finiteVolume`, `meshTools`, `sampling` and
`lagrangian` libraries. The reconstructed planes are read from the registered
fields `interfaceNormal.<phase>` and `interfaceCentre.<phase>`, so it works
with `interIsoFoam` and with solvers using a modified copy of `geometricVoF`,
without linking to either.

## Build

```sh
./Allwmake        # libsurfacePosition in $FOAM_USER_LIBBIN
```

## Usage

```
surfacePosition1
{
    type            surfacePosition;
    libs            (surfacePosition);

    // Vertical lines through locations, across the whole mesh
    locations       ((1.0001 0.0001 0.0001));
    direction       (0 0 1);            // optional, default -g
    // or explicit lines; the first interface from the start is reported
    // lines        (((1 0 0) (1 0 1)));

    alpha           alpha.water;        // optional
    method          auto;               // reconstructed, plic, alpha, auto
    alphaTol        1e-8;               // optional

    checkLineCells      false;          // tests: compare line cells with all cells
    fullSearchInterval  0;              // full search every N mesh changes

    writeControl    timeStep;
    writeInterval   1;
}
```

The output `postProcessing/<name>/<time>/surfacePosition.dat` holds, per
line, the distance along the line, the height above the location (locations
only) and the point; a line without an interface is written as `nan`. See the
header of `surfacePosition.H` for details.

## Tests and report

`test/` holds a test suite (accuracy against an exact wave, special cases on
faces, edges, processor and cyclic patches, moving meshes, cost in a large
parallel case) that runs in about two minutes, and `test/report/` a LaTeX
report on the implementation and the tests, generated from the results of the
suite. See `test/README.md`.
