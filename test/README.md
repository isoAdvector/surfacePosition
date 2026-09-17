# Tests of the surfacePosition function object

`surfacePosition` (in the directory above) reports the position of the
interface along user-specified lines. These tests compare it
with OpenFOAM's `interfaceHeight` for accuracy, robustness and cost.

```sh
./Allrun      # all tests, about 2 minutes; each writes a summary log (see below)
./Allclean
```

Each test can also be run on its own with its `Allrun`. Its cases are
copies of base cases (`baseCase`, or `amrBase` and `motionBase` in
`movingMesh`), each a complete case with its own `Allrun` (serial, or
`Allrun -parallel`) and `Allclean`; a base case runs on its own with the
test's settings. A test's `Allrun` copies the base cases, sets the parameters
that vary between the copies (mesh, water, number of processors) and runs the
copies with their `Allrun`; a test's `Allclean` removes the copies. The
OpenFOAM v2512 environment must be sourced and `libsurfacePosition` built
(`../Allwmake`).
All flow cases use OpenFOAM's `interIsoFoam`.

| Test | What it checks | Summary |
|---|---|---|
| `waveAccuracy` | Exact stream function wave (`streamFunFoam`) on a coarse (38x8) and a medium (95x20) periodic flume; 12 vertical lines at cell centres, a quarter cell off and 2 % of a cell from a face; 40 samples over one wave period; serial and 4 processors | `log.compare` |
| `specialCases` | Flat interfaces inside cells, on internal faces and on processor faces, a vertical water slab with interfaces on cyclic and processor faces, and two layers of water; the three modes accumulated, highest and lowest; lines on faces, edges, processor and cyclic patches, locations far outside the mesh, explicit lines entering through cyclic faces; serial and 8 processors; exact answers | `log.check` |
| `movingMesh` | Standing wave on a 20x20 mesh with adaptive refinement (`amrBase`: 3D, 2 cells across, `plicRDF`, 400 time steps with the mesh changing at most of them; serial and 2 processors) and with a solid body oscillation of the mesh (`motionBase`: 2D, 50 steps; serial and 4 processors), serial and parallel at the same time; the lines lie in the plane y = 0, on internal faces; the line cells after every mesh change are compared with all cells (`checkLineCells`). On 4 processors the amr case loses liquid through an OpenFOAM fault unrelated to `surfacePosition`, hence 2 | `log.compare` |
| `largeTiming` | 3D standing wave in a cube with 100^3 cells on 18 processors (the last 18 cores; the other tests run on the remaining cores at the same time), 50 lines; one run with `surfacePosition` and `interfaceHeight`, timed per function object with OpenFOAM's profiling over 10 steps after 5 start-up steps | `log.timing` |

`compatibility/Allrun` (not run by `./Allrun`) builds `surfacePosition` with
every OpenFOAM.com installation with `geometricVoF` in `~/OpenFOAM` (or
`FOAM_INSTALLS`) and runs `specialCases` and `movingMesh` with each; it writes
`compatibility/log.compatibility` and takes about eight minutes.

`waveAccuracy` needs `streamFunFoam` and `libstreamFunctionWaveFvPatchField`
from the `streamFunctionWave` repository, found next to this repository or
next to a directory above it (or `STREAMFUNCTIONWAVE_SRC` pointing to its
`src/streamFunctionWave`).
`largeTiming` binds one process per core (through a rank file) with local
memory; set `NPROCS`, `NCELLS`, `NSKIP`, `NSTEPS` and `CPUSET` to change the defaults.

## Report

`report/` holds a LaTeX report on the implementation and the tests, and the
template of the repository's `README.md`, a Markdown version of the report.
After `./Allrun`, `report/Allrun` writes the numbers and tables from the
`results.json` files and case data (`makeData.py`), renders the figures from
the cases with ParaView (`makeFigures.py`, `pvbatch`; set `PVBATCH` if it is
not on the path), builds `report/report.pdf` and writes `../README.md`
(`makeReadme.py` from `report/README.md.in`). The figures, the PDF and the
README are committed, so commit them after a new run.
