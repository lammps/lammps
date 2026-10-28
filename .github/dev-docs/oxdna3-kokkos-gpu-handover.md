# oxDNA3 KOKKOS GPU handover: build, run, and A/B the performance branch

Hands-on guide for benchmarking branch `oxdna3KK-kk-perf` (stanmoore1/lammps)
on a GPU.  What each commit does, its expected effect and exactness are in
`oxdna3-kokkos-perf-commits.md` (same directory); this file only says how to
build, run and compare.  Goal: close the gap to the standalone oxDNA CUDA code
on the dilute oligomer benchmark (about 2.7x on an RTX 4090 before this work)
without losing the dense polybrick case.

## 0. Checklist

1. Build the standalone oxDNA (section 2) and LAMMPS at the baseline and at
   the branch HEAD (section 1), in the same precision.
2. Run oligomer and polybrick for all three; save the logs (section 3).
3. Check that LAMMPS HEAD matches the baseline energies (section 5).
4. If HEAD is faster: A/B the compile-time switches and runtime options
   (section 4), then bisect by commit groups (section 6).
5. Profile the best configuration against standalone (section 7).
6. Report per configuration: time per step (or ns/day), and the Pair, Bond,
   Neigh, Comm, Modify and Other lines of the LAMMPS timing breakdown.

## 1. LAMMPS with KOKKOS for GPUs

### Get the code

```bash
git clone https://github.com/stanmoore1/lammps.git lammps-oxdna
cd lammps-oxdna
git fetch origin oxdna3KK-kk-perf oxdna3KK-kk-fixes
# one worktree per version, so builds can coexist
git worktree add ../lmp-base origin/oxdna3KK-kk-fixes    # baseline
git worktree add ../lmp-head origin/oxdna3KK-kk-perf     # this branch
```

The branch sits on top of kk-fixes 2afff57fb5 (neighbor list trimming on by
default, fp32/fp64 conversions removed).  Always compare against that base,
not an older kk-fixes, or its speedups are credited to this branch.

### Configure and build (CUDA, RTX 4090 = Ada, sm_89)

```bash
cd ../lmp-head
cmake -S cmake -B build-cuda -G Ninja \
  -C cmake/presets/kokkos-cuda.cmake \
  -D Kokkos_ARCH_ADA89=on \
  -D KOKKOS_PREC=mixed \
  -D PKG_CG-DNA=on -D PKG_MOLECULE=on -D PKG_ASPHERE=on \
  -D BUILD_MPI=off -D BUILD_OMP=off \
  -D CMAKE_BUILD_TYPE=Release \
  -D Kokkos_ENABLE_LIBDL=on \
  -D CMAKE_CXX_COMPILER_LAUNCHER=ccache
cmake --build build-cuda -j 16
# executable: build-cuda/lmp
```

- Other GPUs: RTX 3050 (Ampere) `-D Kokkos_ARCH_AMPERE86=on`; RX 7900 XT
  use `-C cmake/presets/kokkos-hip.cmake -D Kokkos_ARCH_AMD_GFX1100=on`
  instead of the CUDA preset and arch.
- `KOKKOS_PREC`: `double` (default), `mixed` or `single`.  Standalone oxDNA
  runs CUDA in mixed precision by default, so `mixed` is the fair comparison;
  use the same setting for base and head, and state it in every report.
  Consumer GPUs run FP64 at 1/64 rate, so double vs mixed matters a lot.
- `Kokkos_ENABLE_LIBDL=on` is needed for the Kokkos profiling tools (section 7).
- `BUILD_MPI=off` is fine for one GPU; with MPI use `-D BUILD_MPI=on` and run
  one rank per GPU.
- Build the baseline worktree `../lmp-base` with the identical command.

### Run

```bash
build-cuda/lmp -k on g 1 -sf kk -pk kokkos neigh full newton off -in in.oligo
```

`neigh full newton off` is the KOKKOS default on GPUs; also try
`neigh half newton on` (section 4).  The CG-DNA KOKKOS styles need the
internal fixes they create themselves; nothing else changes in the input.

## 2. Standalone oxDNA (CUDA)

```bash
git clone https://github.com/lorenzo-rovigatti/oxDNA.git
cd oxDNA && mkdir build && cd build
cmake .. -DCUDA=ON -DCUDA_COMMON_ARCH=OFF    # OFF: compile for the local GPU
make -j 16
# executable: build/bin/oxDNA ; run: build/bin/oxDNA input_file
```

- It is compiled with `-use_fast_math`; CUDA runs in mixed precision
  (`backend_precision` defaults to mixed with `backend = CUDA`).
- CUDA input keys that change performance: `use_edge` (one thread per pair;
  default false), `CUDA_list = verlet`, `verlet_skin`,
  `cells_auto_optimisation`, `CUDA_sort_every` (Hilbert sorting, default 0 =
  off), `threads_per_block`.  Use the settings of the published comparison
  and keep them fixed; note them in the report.
- The oligomer and polybrick inputs for both codes are in the benchmarking
  repository lrussell676/oxDNA-KOKKOS-LAMMPS-Benchmarking, directory `oxDNA3`.

## 3. Benchmark runs

- Run each case at least 3 times and take the median; discard the first run
  after a build (JIT/driver warm-up).
- Use enough steps that setup is negligible (e.g. 10^4 to 10^5 steps) and
  the same thermo/dump frequency in both codes (output steps sync to host).
- LAMMPS reports `Performance: ... timesteps/s` and the `MPI task timing
  breakdown` (Pair/Bond/Neigh/Comm/Modify/Other).  Kernels run
  asynchronously, so the breakdown is only accurate with the environment
  variable `CUDA_LAUNCH_BLOCKING=1` (see the KOKKOS section of the manual);
  that slows the run, so use it only for a separate breakdown run and make
  sure it is unset for the headline number.

Check the LAMMPS benchmark inputs for these settings, which cost a lot on
GPUs independent of the pair styles (details in the commits guide):
`fix nve/dotc/langevin` (no KOKKOS version: host round trips every step;
use `fix nve/asphere` + `fix langevin ... angmom`), atom sorting (add
`atom_modify sort 0 0.0`: sorting falls back to the host for this atom
style), `fix balance` and `fix print` on one GPU, `neigh_modify every 1 check
yes` (host wait every step), and `bond_style hybrid` with a single sub-style.

## 4. Options to A/B on the branch

### Compile-time switches (need a rebuild)

All in `src/KOKKOS/mf_oxdna_kokkos.h`, set through `CMAKE_CXX_FLAGS`.  Use
one build directory per variant so that switching does not rebuild
everything:

```bash
V=nofuse; F="-DOXDNA_KK_FUSE_HBXSTK=0"
cmake -S cmake -B build-cuda-$V -G Ninja -C cmake/presets/kokkos-cuda.cmake \
  -D Kokkos_ARCH_ADA89=on -D KOKKOS_PREC=mixed \
  -D PKG_CG-DNA=on -D PKG_MOLECULE=on -D PKG_ASPHERE=on \
  -D BUILD_MPI=off -D BUILD_OMP=off -D CMAKE_BUILD_TYPE=Release \
  -D CMAKE_CXX_COMPILER_LAUNCHER=ccache -D CMAKE_CXX_FLAGS="$F"
cmake --build build-cuda-$V -j 16
```

| Variant | Flags | What it tests |
|---------|-------|---------------|
| default | (none) | fused hbond+oxdna3/xstk kernel, one thread per screened pair |
| nofuse | `-DOXDNA_KK_FUSE_HBXSTK=0` | separate hbond and xstk kernels |
| twophase | `-DOXDNA_KK_TWO_PHASE=1` | radial prefilter + compacted hbond/xstk kernels (no fusion) |
| peratom | `-DOXDNA_KK_SCREENED_PER_ATOM=1` | one thread per atom over its screened pairs (no fusion) |
| lb-pair | `-DOXDNA_KK_PAIR_MAXT=128 -DOXDNA_KK_PAIR_MINB=2` | launch bounds of hbond/xstk/coaxstk/fused (fused uses 123 registers) |
| lb-atom | `-DOXDNA_KK_ATOM_MAXT=128 -DOXDNA_KK_ATOM_MINB=2` | launch bounds of excv/dh (default 64/1 on CUDA) |
| lb-bond | `-DOXDNA_KK_BOND_MAXT=128 -DOXDNA_KK_BOND_MINB=4` | launch bounds of stk/fene (stk uses 117 registers) |

Sweep MAXT in {64, 128, 256} and MINB in {1, 2, 4} for the kernel class
that dominates the profile.  `TWO_PHASE` and `SCREENED_PER_ATOM` cannot be
combined (compile error by design).

### Runtime options (no rebuild; add to the input or command line)

| Option | Values to try |
|--------|---------------|
| neighbor list | `-pk kokkos neigh full newton off` (default) vs `neigh half newton on` |
| trimming | default on; `pair_modify neigh/trim no` to see its effect |
| skin | `neighbor 0.3 bin` ... `neighbor 1.0 bin` (smaller skin = shorter lists, more rebuilds) |
| rebuild check | `neigh_modify every 1 delay 0 check yes` vs `every 5` / `every 10` |
| sorting | `atom_modify sort 0 0.0` |
| bonds | `bond_style oxdna3/fene` instead of `bond_style hybrid oxdna3/fene` |

## 5. Correctness check on the GPU

Before timing, confirm the branch gives the same physics as the baseline:

```bash
cd examples/PACKAGES/cgdna/examples/lj_units/oxDNA3/duplex2
../../../../../../build-cuda/lmp -k on g 1 -sf kk -in in.duplex2 -log head.log
<base>/build-cuda/lmp            -k on g 1 -sf kk -in in.duplex2 -log base.log
<any CPU build>/lmp                               -in in.duplex2 -log cpu.log
```

- Step 0 energies (all thermo columns) must agree to about 1e-6 relative in
  mixed precision (1e-12 in double).
- Later steps diverge slowly (Langevin noise, atomics order); compare
  averages, not single steps.  For a trajectory check, replace the Langevin
  thermostat by plain `fix nve/asphere` for 100 steps.
- The branch was verified on CPU with Kokkos Serial on 4 MPI ranks against the
  CPU styles (all oxDNA1/2/3 and oxRNA2 examples, oligomer and dense proxies,
  rerun, minimize, trimmed lists, both sub-style orders), so a mismatch on
  the GPU points to a device-only path: the screened-pair kernels, the fused
  kernel, or the coax list.  Report it with the deck and the thermo lines.

## 6. Bisecting the branch

The commits are ordered so that groups can be benchmarked separately (table
in `oxdna3-kokkos-perf-commits.md`).  Build these points (one worktree each),
all on the same base:

| Point | Commit (hash; subject for `git log --oneline --grep`) | Contains |
|-------|--------------------------------------------------------|----------|
| base | kk-fixes 2afff57fb5 | baseline |
| G1 | a61f378bc4 "one atomic update per atom and component in pair oxdna*/stk" | exact cleanups, early exits, coax list, stack removal (#1-10) |
| G2 | 87f3af21a7 "pack the coefficients of pair oxdna2/coaxstk" | launch bounds, atomics, layout, fused force clear, packed coefficients (#11-18) |
| G3 | 46acb30b6d "fused kernel of pair oxdna3/hbond and oxdna3/xstk" | two-phase option, npair list cutoff, dropped lists, fused kernel (#19-22) |
| head | branch HEAD | second search: excv/dh/record/npair/stk/fene (#23-29) |

```bash
H=$(git log --format=%h --grep="fused kernel of pair oxdna3/hbond" -1 origin/oxdna3KK-kk-perf)
git worktree add ../lmp-g3 $H
```

The hashes are those of the branch as pushed on top of kk-fixes 2afff57fb5;
if the branch is rebased again, use the subjects.  If one group changes the
timing a lot, bisect inside it commit by commit.

## 7. Profiling

- Kokkos simple kernel timer (per-kernel time and call count):
  ```bash
  git clone https://github.com/kokkos/kokkos-tools.git
  cmake -S kokkos-tools -B kt-build && cmake --build kt-build -j
  export KOKKOS_TOOLS_LIBS=$PWD/kt-build/profiling/simple-kernel-timer/libkp_kernel_timer.so
  build-cuda/lmp -k on g 1 -sf kk -in in.oligo
  ```
  Depending on the kokkos-tools version the summary is printed or written to a
  `*.dat` file per process (read with `kp_reader`).  LAMMPS must be built with
  `Kokkos_ENABLE_LIBDL=on`.
- Nsight Systems for launch gaps and host synchronization:
  `nsys profile --stats=true -o oligo build-cuda/lmp ...`.  Look for idle GPU
  time between kernels (launch latency, the per-step `check_distance`
  readback, fences) versus kernel time.
- Nsight Compute for single kernels, e.g. the fused kernel:
  `ncu -k regex:HbXstkFused --set full build-cuda/lmp ...`; check achieved
  occupancy, register limits, and `sm__inst_executed_pipe_fp64` (FP64 work
  left in mixed precision).
- Registers, stack and kernel parameter size without running:
  `cuobjdump -res-usage build-cuda/CMakeFiles/lammps.dir/<path>/pair_oxdna_excv_kokkos.cpp.o`.
- Standalone: `nsys profile` the same way; its step is about 4 kernels
  (first_step, forces, second_step, thermostat) versus about 13-15 for
  LAMMPS, which is what matters most in the dilute case.

## 8. What to report back

- GPU, driver, CUDA version, `KOKKOS_PREC`, standalone input keys.
- Table: case x version/variant -> timesteps/s (median of 3) and speedup vs
  base and vs standalone.
- For the best LAMMPS configuration: top 10 kernels (time, calls) from the
  kernel timer or nsys, and the GPU idle fraction per step.
- Any correctness mismatch from section 5.
