---
applyTo: "src/**,lib/**,unittest/**"
---

# LAMMPS C++ Source Conventions

Rules for C++ code in `src/`, `lib/`, and `unittest/`, both when writing code and when
reviewing it.  General background: https://docs.lammps.org/Modify_requirements.html and
https://docs.lammps.org/Modify_style.html.  Implementation lessons for styles (MPI,
restart, buffers, parsing): `.github/dev-docs/style-implementation-notes.md`.

## Hard rules (enforced by CI or rejected in review)

- **C++17**; format with `.clang-format` in `src/`; 7-bit ASCII only.
- **No variable-length arrays**; use `memory->create()` or `std::vector`.
- **No alternative logical-operator tokens:** use `&&`, `||`, `!`, `^` -- never `and`,
  `or`, `not`, `xor` (breaks MSVC).
- **No `printf()`, `fprintf(screen, ...)`/`fprintf(logfile, ...)`, or C++ iostreams
  (`std::cout`, `std::cerr`) in new code:** use `utils::logmesg()` or `utils::print()`,
  which both support `std::format`-style formatting.  Remove output that was added for
  debugging entirely.
- **No new `error->all(FLERR, "Illegal XXX command")` messages:** use
  `utils::missing_cmd_args()` or a specific message with the error-pointer argument that
  highlights the offending argument, plus `utils::errorurl()` to point to the
  explanations on https://docs.lammps.org/Errors_details.html (see
  https://docs.lammps.org/Developer_notes.html#errors-warnings-and-informational-messages
  and compare with code that is already converted).
- **No `strtok()`, `sscanf()`, `atoi()`, `atof()` and similar for parsing:** use the
  `Tokenizer`/`ValueTokenizer` classes (and a file reader class where possible);
  convert arguments with `utils::numeric()`, `inumeric()`, `bnumeric()`, `tnumeric()`.
- **String formatting with fmtlib** (`fmt::format()`), not `sprintf`.
- **No commented-out code** (debugging leftovers, disabled features, unused
  alternatives) unless the pull request explains it as a placeholder for a planned
  feature.
- **Package file bookkeeping:** new files in package directories go into
  `src/.gitignore` (the legacy make build copies them into `src/`, and those copies
  must not be committed); renamed or removed package files go into `src/Purge.list` so
  `make purge` removes stale copies.

## Conventions

- Parenthesize each operand of chained `&&`/`||` conditionals.
- No two-trip loops for trivial initialization: assign pairs directly
  (`xstyle[0] = xstyle[1] = NONE;`); keep the loop when the body is substantial.
- **Portable numerics:** never compute with, hash, or compare plain `char` values (signed
  on x86-64, unsigned on ARM64 Linux) -- use `signed char`/`unsigned char`/`uint8_t`;
  never rely on exact cancellation like `a*a - b*b == 0.0`, since compilers may contract
  it into a fused multiply-add (the default on ARM64).
- **Unused parameters:** silence `-Wunused-parameter` by commenting out the name
  (`int /*nlocal*/`), never by dropping it; the name documents the argument.
- **clang-format:** keep a style's own header in its own include block (blank line after
  it) or clang-format sorts it into the other includes.  When substantially editing a
  small, simply structured legacy file marked `// clang-format off`, remove the marker
  and reformat; keep it for complex, deeply nested code.
- **Error handling:** `error->all()` when all MPI ranks hit the error, `error->one()` for
  a single rank; `error->warning()` prints on every rank, so guard with `comm->me == 0`
  where a single message is wanted.
- **RAII for C resources:** prefer `SafeFilePtr` (`src/safe_pointers.h`) over raw
  `FILE *`/`fopen` when touching such code.
- **`delete[]` before `utils::strdup()`:** when storing a copied name (variable, region,
  group ID, ...) in a class member, always `delete[]` the member immediately before
  re-assigning it -- even when it is provably still `nullptr`.  Static analysis
  (Coverity) flags the bare assignment as a leak, and the idiom is defensive against
  keywords being parsed twice.
- **MPI stubs:** if a serial build misses an MPI symbol, add it to `src/STUBS/mpi.h`
  instead of special-casing the caller.
- **Block comments:** an embedded `*/` (e.g. in a glob like `gb_*/ga_*`) silently ends a
  `/* ... */` comment; reword or use `//` comments.

## Adding new styles

1. Place `<prefix>_<name>.cpp`/`.h` in `src/` or the package directory, starting from a
   similar existing style (https://docs.lammps.org/Modify_style.html).
2. Update `src/.gitignore` and `src/Purge.list` as described above.
3. Create or update the matching `doc/src/*.rst` page; new public commands and keywords
   need `.. versionadded:: TBD` (see `.github/instructions/documentation.instructions.md`).
   Internal styles (upper-case style names) need no documentation.
4. Examples for styles in a package go under `examples/PACKAGES/<name>/`; top-level
   `examples/` folders are for core styles only.
