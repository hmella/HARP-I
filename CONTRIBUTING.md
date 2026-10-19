# Contributing

Issues and pull requests are welcome. There are two implementations in this
repository and they are checked differently.

## The Python package

```sh
pip install -e 'python[gui,dev]'
cd python && pytest -q
```

621 tests, about two minutes. Tests that need the 3D dataset, the tagged
DICOM exports or a MATLAB reference skip with a message naming what is
missing, so a fresh clone sees no failures.

The interface is tested without a display. On a headless machine set
`QT_QPA_PLATFORM=offscreen`, which is also what continuous integration
should use rather than relying on Qt's fallback.

`HARPI_NO_NUMBA=1` exercises the plain NumPy unwrapper instead of the
compiled one.

Style is checked with ruff, configured in `python/pyproject.toml`:

```sh
cd python && ruff check src tests
```

## The MATLAB toolbox

```matlab
HARPIsetup
runtests('tests')
```

`HARPIsetup` puts DENSEanalysis before HARP-I on the path and reports any
HARP-I function still shadowed. See the search path note in the README.

## Keeping the two in step

The Python port reproduces the MATLAB toolbox to about 1e-9 at every stage,
and the regression harness in `tools/regression/` compares them with a
tolerance of exactly zero. A change to a shared code path has to keep that
green:

```matlab
maxNumCompThreads(1); addpath('tools'); harpi_paths;
addpath('tools/regression'); capture_reference('Output', '/tmp/r.mat')
```

then `compare_reference` against the stored reference with `Tol 0`.

Two departures from the toolbox are deliberate and documented in
`python/README.md`: the temporal fit, and the default rule for phase jumps.
Anything else that differs is a bug in one of them.

## Licence

Contributions are accepted under the Mozilla Public License 2.0, the licence
of this repository. New source files should carry the short notice the other
files have.
