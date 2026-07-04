# Coverity Scan modeling for LAMMPS

This directory holds the [Coverity](https://scan.coverity.com/) *modeling file*
used to reduce false positives in the automated static-analysis runs of LAMMPS
(project `LAMMPS` on scan.coverity.com, driven by
`.github/workflows/coverity.yml`).

## What a modeling file is

`coverity_model.cpp` is **not compiled into LAMMPS**. It contains stub bodies
that describe, to the Coverity analysis engine, the *relevant behavior* of a few
functions the engine would otherwise misinterpret -- using Coverity's
`__coverity_*__` modeling intrinsics (e.g. `__coverity_alloc__`,
`__coverity_free__`, `__coverity_panic__`). The analyzer matches a model to the
real function by signature (namespace + class + name + parameter types), so each
stub must declare exactly the same signature as the corresponding function in the
LAMMPS sources.

A useful mental model: Coverity analyzes call graphs, and a model edits the
analyzer's belief about what a *called function* does. Consequently a model can
only suppress a false positive that flows through a **function call** (allocation,
freeing, non-return, tainting, buffer initialization). It cannot affect a false
positive that arises from inline arithmetic the analyzer sees directly (for
example a "division by zero" on a local variable) -- triage those in the Scan UI
or add a real guard in the code.

## What is (and is not) modeled here

* **`Memory::smalloc` / `srealloc` / `sfree`** (`src/memory.h`) are modeled as
  allocate / reallocate / free. The templated `Memory::create`, `grow`, and
  `destroy` helpers funnel through these three primitives, so modeling the
  primitives covers the whole custom allocator and suppresses spurious
  "resource leak" / "use after free" reports that stem from the analyzer not
  recognizing it as an allocator.

* **`Error::all` / `Error::one` / `Error::done` are deliberately NOT modeled.**
  They are already declared `[[noreturn]]` in `src/error.h`, which `cov-analyze`
  honors directly. An in-source attribute is preferable to a model: it also
  informs the compiler, clang-tidy, and every other tool, and it is version
  controlled and reviewed. **Rule of thumb: when we own the code, annotate it in
  the source (`[[noreturn]]`, `[[nodiscard]]`, `assert`, explicit checks) rather
  than adding a model.** Reserve this file for code we cannot annotate or where no
  suitable attribute exists.

## Building the model database

Use the same Coverity analysis toolkit that the CI downloads:

```bash
cov-make-library -of user_models.xmldb coverity_model.cpp
```

This produces `user_models.xmldb`.

## Deploying to Coverity Scan

LAMMPS uses the hosted Scan service, which runs `cov-analyze` *server-side* after
each uploaded build (see `.github/workflows/coverity.yml`: the workflow does
`cov-build` -> tar -> upload, not a local analysis). The model is therefore
registered with the **project**, not referenced from the workflow:

1. Sign in to <https://scan.coverity.com/> as a project administrator for `LAMMPS`.
2. Open the project's analysis settings and upload the modeling file
   (`coverity_model.cpp`; Scan compiles it for you). If a compiled database is
   requested instead, upload the `user_models.xmldb` produced above.
3. The next analysis after a build upload applies the model.

For a local/on-prem analysis (not the current LAMMPS setup) the equivalent is:

```bash
cov-analyze --dir cov-int --user-model-file user_models.xmldb ...
```

## Maintenance

* **Keep signatures in sync.** If a modeled function's signature changes in the
  LAMMPS sources, update the stub here and re-upload. A drifted signature stops
  matching silently -- with no error -- and the model goes inert.
* **Keep the file minimal and accurate.** A model is production logic *for the
  analyzer*: an inaccurate model (marking a function non-returning when it can
  return, or freeing without reallocating) silently hides real defects. Review
  every change here as carefully as shipped code.
* **Prefer in-source annotations.** Before adding a model for LAMMPS-owned code,
  check whether a `[[noreturn]]` / `[[nodiscard]]` attribute, an `assert`, or an
  explicit check expresses the same fact in the source instead.
