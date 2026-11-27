# Backporting Fixes to the `maintenance` Branch

The `maintenance` branch collects bug fixes for the current stable release; the
maintainers turn it into the stable update releases (see "LAMMPS Stable Update
Release" in `.github/release_steps.md`).  The maintainers decide what is backported.
This file describes how to prepare and check such backports.

## What qualifies

- Fixes for defects that the stable release has: wrong results, crashes, memory
  errors, build failures, wrong documentation.  No new features, refactoring,
  performance work, or new accelerator variants.
- Do not discount a fix because the branch will be retired with the next stable
  release: many users keep running a stable release long after its last update.
- When a fix builds on code that the stable release does not have (a new class, a
  changed interface), adapt it to the old code, or leave it out and report it.  Do not
  backport the dependency only to make the patch apply.
- Binary restart files written by the stable release must stay readable.  Do not
  change what a style writes to or reads from a restart file, even when the fix on
  `develop` does: adapt the fix to the old layout.

## How

- One commit per commit of `develop`, made with `git cherry-pick -x <hash>` wherever
  possible.  This keeps the original author and message and adds the line "(cherry
  picked from commit ...)", by which a later review recognizes what is already
  backported.  If the original message has AI attribution trailers, drop them (see
  the policy in `.github/copilot-instructions.md`).
- Unit tests are not backported: the test infrastructure changes during a development
  cycle, so test files of `develop` may not work with the older one.  Leave the
  `unittest/` part of a commit out and say so in the commit message.  The exception:
  when a fix changes the result of a test that `maintenance` already has, the updated
  reference data goes into the same commit.  Check it with the test driver of
  `maintenance`; the tolerances in the two branches can differ.
- `.. versionadded::` and `.. versionchanged::` directives of a backported change use
  `TBD` or `<stable version>-backport` (e.g. `30Sep2026-backport`), never the version
  of the feature release that first contained the change.
- The pull request that merges `maintenance` into `stable` is the draft of the release
  notes: add a line for each backported fix to its description, with the number of
  the pull request it comes from.

## Finding candidates

Compare the commits that `develop` gained since the last review with `maintenance`:

1. `git cherry maintenance develop`, identical subject lines, and the "(cherry picked
   from commit ...)" lines of `maintenance` identify what is already there.
2. For every other commit, test file by file against an export of `maintenance`
   (`git archive maintenance | tar -x -C <scratch>`) whether its change is already
   present (`git apply --check -R`) or applies (`git apply --check`).  A commit that
   only touches files which `maintenance` does not have belongs to a new feature.
3. Read the subjects of the rest and check each fix against the sources of
   `maintenance`: does the defect exist there?  A fix for a regression that entered
   `develop` after the stable release is not a candidate.
4. Report the candidates with the defect each one fixes, and separately those that
   need a decision (behavior changes, fixes that need adapting).

## Verification

- Apply the selected fixes to a scratch copy first, in the order of `develop`.  Build
  with as many packages as possible, without new warnings, and run the unit tests of
  `maintenance` itself.  Build the manual when documentation changed (`make html` and
  `make spelling` in `doc/`).
- Do not run the regression tests of `develop` (`tools/regression-tests` with the
  current `examples/`) against a `maintenance` executable: the inputs have changed
  along with `develop`, and the differences are not defects.
- For a fix that touches restart files, write them with an unpatched build of
  `maintenance` and read them with the patched one.
