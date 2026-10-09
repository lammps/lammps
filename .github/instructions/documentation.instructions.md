---
applyTo: "doc/**"
---

# LAMMPS Documentation Guidelines

## Language and format

- reStructuredText (`.rst`) in `doc/src/`; American English spelling; 7-bit ASCII only.
- Wrap code examples in `.. code-block::` with the appropriate language
  (`LAMMPS`, `bash`, `c++`, `python`); use `.. note::` for important remarks and
  `.. warning::` for critical warnings.
- Avoid computer-science jargon and ambiguous abbreviations in user-facing text; the
  audience is researchers.  Never abbreviate to "FP" (ambiguous between false positive,
  file pointer, and floating-point) -- spell the term out.
- Always write "KOKKOS package" (all uppercase), never "Kokkos library" or "KOKKOS
  library" -- even where the library is technically meant.  Single exception: the build
  documentation (`Build_extras` KOKKOS section) that explicitly discusses the bundled,
  downloaded, or external Kokkos *library*.
- Never add legacy-status disclaimers ("no longer actively developed", "superseded
  by ...") to manual pages: they undermine user trust and confuse inexperienced users.
  Steer users toward a modern alternative by improving that style and its docs, not by
  talking down the old one.  Legacy status belongs only in developer-facing notes and
  test nomenclature (e.g. a `legacy` filename token or tag).

- Do not describe internal build machinery (patches applied to downloaded library
  sources, timestamp or cache workarounds) in the user manual: the build integration
  should be transparent, and such notes only create concerns.  Keep these details in
  CMake comments, commit messages, and pull request notes; document only what users
  must act on (settings, variables, warnings they may see).

## Command-line conventions

Use one consistent form for command lines in the manual and in README files:

- CMake configuration is a single command from the top-level LAMMPS folder:
  `cmake -S cmake -B build [-C cmake/presets/NAME.cmake ...] [-D NAME=value ...]`.
  Do not write `mkdir build`, `cd build`, `cmake [...] ../cmake`; refer to presets
  as `cmake/presets/...`.  The default build folder is `build`; variants are named
  `build-<purpose>` (e.g. `build-kokkos`).
- To change settings of an existing build, repeat the full
  `cmake -S cmake -B build -D NAME=value` command.  `cmake -B build ...` without
  `-S` fails for LAMMPS, because the source folder then defaults to the current
  folder, which has no `CMakeLists.txt`.
- Write CMake variables as `-D NAME=value` (with a space) and boolean values as
  `on`/`off`.
- Compile with a plain `cmake --build build` (`--target NAME` for other targets),
  install with `cmake --install build`, and run tests with
  `ctest --test-dir build`.  No bare `make`, `make -j N`, or `make -C build` in
  CMake instructions (they fail with the Ninja generator).  Parallel compilation
  (`-j N` / `--parallel N`) is explained once in the CMake build docs instead of
  being repeated in every example.
- `-G Ninja` appears only in the section on faster compilation with Ninja and
  ccache; all other examples stay generator-neutral.  Prefer explicit command-line
  flags over environment variables or shell-profile settings.
- The LAMMPS executable is `lmp`.  Use `lmp_<machine>` names only where a legacy
  make build is the topic (`make mpi` produces `lmp_mpi`), for the
  `LAMMPS_MACHINE` setting, or for actual binary names of distribution packages.
- Always pass input files with `-in in.file` (long form).  Use I/O redirection
  (`lmp < in.file`) only where reading from standard input is the topic.

## versionadded / versionchanged policy

- New publicly visible commands, styles, and added keywords require
  `.. versionadded:: TBD`; changed behavior of existing commands/keywords requires
  `.. versionchanged:: TBD`.  `TBD` is replaced with the release version string during
  release preparation.
- Place the directive in front of the paragraphs documenting the new or changed
  functionality, or following the "Description" headline for completely new commands.
- The directives are standalone markers with an EMPTY body: explanatory prose goes in
  the *following* normal paragraph at column 0, NOT indented under the directive.
  Do not "fix" existing directives by indenting the following paragraph.
- This does NOT apply to internal styles (upper-case style names, which need no
  documentation) or to merely adding an accelerated variant of an existing style.
  For accelerator variants, add the code letter (`g`/`i`/`k`/`o`/`t`) to the
  respective `Commands_<type>.rst` page instead -- in alphabetical order with no
  commas (`ko`, not `o,k`) -- and add the `/suffix` index entry plus an
  `Accelerator Variants:` line to the per-style `.rst` page.
- Restrictions and accelerator-variant notes belong in the per-style `.rst` file.
- When commands or keywords are added or changed, check whether example inputs use
  them and need updating.

## Building and validating

```bash
cd doc
make html          # must complete without sphinx/docutils warnings or bad references
make pdf           # requires pdflatex; must also complete cleanly
make spelling      # must report no issues
make anchor_check  # duplicate anchors
make style_check   # style lists complete
make check         # anchor/style/package/char/role checks in one target
```

- Spelling exceptions (author surnames, acronyms, code fragments) go into
  `doc/utils/sphinx-config/false_positives.txt` (note: underscores in the file name).
- `make spelling` resolves the word list relative to the SOURCE dir by first copying
  `false_positives.txt` into `doc/src/`.  A bare `sphinx-build -b spelling` without
  that copy uses NO word list and flags everything -- always go through the Makefile.
- The doc build runs sphinx with parallel reading (`-j`); the spelling target keeps a
  dedicated `doctrees-spelling` cache directory.  If sphinx claims "configuration has
  changed" on every run, suspect a config value mutated at `builder-inited` time
  (handled in `conf.py.in`; keep new config values immutable).

## Figures and diagrams

- Editable SVG sources for generated figures live in `doc/src/SVG/`; they rasterize to
  `doc/src/JPG/*.png`.  The re-render command is recorded in each SVG's header comment
  (Inkscape with `--export-area-drawing` and transparent background).
- Inkscape does not support the `feDropShadow` filter shorthand -- build drop shadows
  from primitives (`feGaussianBlur` + `feOffset` + `feFlood` + `feComposite` +
  `feMerge`), and avoid `--` inside XML comments.
- Verify figures are ASCII-clean (`grep -P '[^\x00-\x7F]'`) and well-formed
  (`xmllint --noout`) before committing.
