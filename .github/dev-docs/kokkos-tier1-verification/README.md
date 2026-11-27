# KOKKOS Tier-1 verification decks

Hand-written input decks that exercise every style ported to KOKKOS in the
Tier-1 pass, for styles whose behavior the `unittest/force-styles` YAML
fixtures do not cover (fixes, computes, regions, minimizers, granular pair
styles) and, as a second opinion, for those they do.

The point of running them is that the `package kokkos` defaults are not the
settings a GPU run gets.  Without a GPU, LAMMPS defaults to `comm no`,
`sort no`, `atom/map no`, `neigh half`, `newton on` and `gpu/aware off`,
which routes around the device communication, sorting and atom-map paths.
`run_checks.sh` therefore drives every KOKKOS run with

    -pk kokkos neigh full newton off comm device sort device atom/map device gpu/aware on

which is what a GPU build would pick, and falls back to `neigh half newton on`
(third field of a `CASES` line) for the styles that require a half neighbor
list.

## Building

`tier1-packages.cmake` turns on the packages these decks need and nothing
else.  The two executables the passes below want, from the top of the
repository:

    cmake -S cmake -B build-sync -G Ninja \
      -C .github/dev-docs/kokkos-tier1-verification/tier1-packages.cmake \
      -D PKG_KOKKOS=on -D Kokkos_ENABLE_SERIAL=on -D BUILD_MPI=on \
      -D CMAKE_BUILD_TYPE=Release -D KOKKOS_DEBUG_SYNC=on \
      -D BUILD_SHARED_LIBS=on -D DOWNLOAD_POTENTIALS=off -D ENABLE_TESTING=on

    cmake -S cmake -B build-poison -G Ninja \
      -C .github/dev-docs/kokkos-tier1-verification/tier1-packages.cmake \
      -D PKG_KOKKOS=on -D Kokkos_ENABLE_SERIAL=on -D BUILD_MPI=on \
      -D CMAKE_BUILD_TYPE=RelWithDebInfo \
      -D KOKKOS_DEBUG_SYNC=on -D KOKKOS_DEBUG_SYNC_ASAN=on \
      -D BUILD_SHARED_LIBS=on -D DOWNLOAD_POTENTIALS=off

`step0` and `mpi` only compare styles against each other and run in any
build; the detector passes need `build-sync`, and `poison` needs
`build-poison`.

## Usage

    ./run_checks.sh step0       # thermo: plain CPU styles vs KOKKOS styles, same binary
    ./run_checks.sh detect      # LMP_KOKKOS_WATCH / LMP_KOKKOS_STALE detectors
    ./run_checks.sh mpi         # step0 again under two MPI ranks
    ./run_checks.sh detect_mpi  # detect again under two MPI ranks
    ./run_checks.sh poison      # LMP_KOKKOS_POISON, needs the ASan build

`step0` runs each deck twice from the same executable, once with the plain
styles and once with `-sf kk`, and compares every thermo column with
`cmp.py`.  A column is equal when `|a-b| <= 1e-8 + 1e-6*max(|a|,|b|)`; the
reported number is the ratio to that tolerance, so anything above 1.0 fails.

The `mpi` and `detect_mpi` passes repeat the same cases under
`mpirun -np 2`, which is what exercises the exchange and border
communication; `MPIRUN` can be overridden in the environment.

`detect` needs an executable built with `-D KOKKOS_DEBUG_SYNC=on`; see
`../kokkos-sync-debugging.md` for what the reports mean and for the
poison-mode build that pinpoints a stale access.  `poison` runs that build
(`build-poison/lmp`, or `POISON_LMP`) with `LMP_KOKKOS_POISON=1`, which
makes the non-current side of every view unreadable so AddressSanitizer
stops at the offending instruction; it is the only detector that sees a
read through a plain pointer.  Judge its output by
diffing the reported labels (array plus the routine that touched it)
against a run of a style of the same shape that is known to be clean, not
by their absolute number: the KOKKOS infrastructure reports a fixed set of
labels for every input.  `LMP` and `OUT` can be set in the environment;
`LMP` defaults to `build-sync/lmp` at the top of the repository.

## Layout

`common_lj.mod` (LJ melt), `common_sphere.mod` (sphere-style melt) and
`mol_head.mod` (`data.fourmol`, read from `unittest/force-styles/tests`) are
the shared setups; each `in.*` deck includes one of them and cycles through a
family of styles with an index variable, so one deck covers a whole work
package.  `common_lj.mod` leaves the boundary in an index variable `bnd` so a
deck that needs a non-periodic or shrink-wrapped z can set it before the
include.
