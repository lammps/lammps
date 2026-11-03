---
applyTo: "tools/regression-tests/**,examples/**"
---

# Regression Tests and Example Inputs

The inputs in `examples/` double as regression tests: `tools/regression-tests/run_tests.py`
runs them and compares their thermo output (PotEng, TotEng, Press, Temp, E_vdwl) with
the reference log files next to them.  CI runs a quick subset after code review; local
runs are rarely needed.  Full documentation: `tools/regression-tests/README` and
`tools/regression-tests/REPORTING.md`.

```bash
python3 -m venv testenv && source testenv/bin/activate
pip install numpy pyyaml
python3 tools/regression-tests/run_tests.py --lmp-bin=build/lmp \
    --config-file=tools/regression-tests/config_quick.yaml --examples-top-level=examples
```

## Conventions for example inputs

- Input scripts must be named `in.<name>`; other names are never collected.
- Reference logs are named `log.<date>.<name>.<compiler>.<N>` (e.g.
  `log.8Apr21.melt.g++.4`).  `<name>` must match the input script, or the input is
  silently left unchecked; `<N>` counts MPI ranks times OpenMP threads.
- Make results independent of the number of MPI ranks: use `velocity ... loop geom` when
  atoms come from `create_atoms` (the default `loop all` is fine with `read_data`).
- Keep runs short; production-sized runs hit the test timeout.
- Examples for styles in a package go under `examples/PACKAGES/<name>/`.

## Gotchas

- Local runs execute inside `examples/` and overwrite tracked output files written by the
  inputs (logs, images, data files): run `git restore examples/` and remove the untracked
  `log.*` debris before committing.
- `progress.yaml`/`reference.yaml` are keyed by the input's basename, so duplicate names in
  different folders collide there (the JUnit and JSON results use full paths).
