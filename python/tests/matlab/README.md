# MATLAB references

`gen_refs.m` produces the stage-by-stage snapshot the Python tests compare
against. It is not committed (~60 MB); the tests skip cleanly without it.

```sh
matlab -batch "run('$(git rev-parse --show-toplevel)/python/tests/matlab/gen_refs.m')"
```

Run single-threaded. MATLAB is only bit-reproducible at a fixed thread count,
and `gen_refs.m` sets `maxNumCompThreads(1)` for that reason.

The tests also accept `tools/regression/refs/py_ref_2D.mat` if you generated
one there with the MATLAB regression harness.
