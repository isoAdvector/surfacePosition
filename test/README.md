# Tests of the surfacePosition function object

`surfacePosition` (in the directory above) reports the position of the
interface along user-specified lines. These tests compare it
with OpenFOAM's `interfaceHeight` for accuracy, robustness and cost.

```sh
./Allrun      # all tests, about 2 minutes; each writes a summary log (see below)
./Allclean
```

Each test can also be run on its own with its `Allrun`. The OpenFOAM v2512
environment must be sourced and `libsurfacePosition` built (`../Allwmake`).
All flow cases use OpenFOAM's `interIsoFoam`.

| Test | What it checks | Summary |
|---|---|---|
| `waveAccuracy` | Exact stream function wave (`streamFunFoam`) on a coarse (38x8) and a medium (95x20) periodic flume; 12 vertical lines at cell centres, a quarter cell off and 2 % of a cell from a face; 40 samples over one wave period; serial and 4 processors | `log.compare` |
| `specialCases` | Flat interfaces inside cells, on internal faces and on processor faces, and a vertical water slab with interfaces on cyclic and processor faces; lines on faces, edges, processor and cyclic patches, locations far outside the mesh, explicit lines entering through cyclic faces; serial and 8 processors; exact answers | `log.check` |
| `movingMesh` | 2D standing wave on a 20x20 mesh with adaptive refinement (`amr`, 400 time steps, some 180 mesh changes) and with a solid body oscillation of the mesh (`motion`, 50 steps), serial and 4 processors at the same time; the lines lie in the plane y = 0, on the internal faces that refinement creates; the line cells after every mesh change are compared with all cells (`checkLineCells`) | `log.compare` |
| `largeTiming` | 3D standing wave in a cube with 100^3 cells on 18 processors (the last 18 cores; the other tests run on the remaining cores at the same time), 50 lines; one run with `surfacePosition` (reconstructed and plic) and `interfaceHeight`, timed per function object with OpenFOAM's profiling over 10 steps after 5 start-up steps | `log.timing` |

`waveAccuracy` needs `streamFunFoam` and `libstreamFunctionWaveFvPatchField`
from the `streamFunctionWave` repository, found next to this repository or
next to a directory above it (or `STREAMFUNCTIONWAVE_SRC` pointing to its
`src/streamFunctionWave`).
`largeTiming` binds one process per core (through a rank file) with local
memory; set `NPROCS`, `NCELLS`, `NSKIP`, `NSTEPS` and `CPUSET` to change the defaults.

## Report

`report/` holds a LaTeX report on the implementation and the tests. After
`./Allrun`, `report/Allrun` writes its numbers and tables from the
`results.json` files and case data (`makeData.py`), renders the figures from
the cases with ParaView (`makeFigures.py`, `pvbatch`; set `PVBATCH` if it is
not on the path) and builds `report/report.pdf`.
