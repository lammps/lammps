---
applyTo: "src/GPU/**,lib/gpu/**"
---

# GPU Package Rules (always-on)

These rules apply to every change under `src/GPU/` and `lib/gpu/`, when writing and
when reviewing code.  The how-to for a new `/gpu` pair style (files, choice of the
template, kernel conventions, test recipe) is `.github/dev-docs/gpu-porting-guide.md`;
what remains to port is in `.github/dev-docs/accelerator-porting-backlog.md`.

## Key rules

1. **One kernel source for three back ends.**  A `lib/gpu/lal_*.cu` file is compiled
   for CUDA, HIP, and OpenCL.  OpenCL C has no references and no templates: use the
   types `numtyp`/`acctyp`, write literals as `(numtyp)1.0`, use the `ucl_*` math
   functions and `ucl_inline` helpers.  OpenCL kernels are compiled only at run time,
   so a kernel counts as portable only after an OpenCL run.
2. **Kernel argument lists are not checked at compile time.**  The `run(&a, &b, ...)`
   call in `loop()` must match the kernel signature in order and count, for the
   generic kernel and for the `_fast` kernel.  A mismatch shows up as wrong numbers or
   a failure to launch.
3. **A style needs all of its parts:** `lib/gpu/lal_<name>.{h,cpp,cu}`,
   `lib/gpu/lal_<name>_ext.cpp`, `src/GPU/pair_<style>_gpu.{h,cpp}`, the declarations
   in `lib/gpu/include/lammps_gpu.h`, two `action` lines in `src/GPU/Install.sh`, and
   the documentation entries (the `g` letter in `Commands_pair.rst`, the `.. index::`
   line and the "Accelerator Variants" line of the style's page).  New files need no
   `src/Purge.list` entry.
4. **The wrapper's `init_style()` replaces the one of the CPU style.**  It must repeat
   all error checks and set every member that the inherited `single()` uses, with the
   meaning the CPU style gives it.  Pass GPU-specific transforms (inverse
   denominators, squared radii) as arguments instead of overwriting members.
5. **`compute()` must rebuild the topology lists:** `if (atom->molecular !=
   Atom::ATOMIC && neighbor->ago == 0) neighbor->build_topology();`.  With a GPU
   pair style `Neighbor::build()` no longer does this, and without the line the bonded
   styles use stale indices after the first reneighboring.
6. **Parameters that `fix adapt` can change need a reinit path:** override
   `init_one()` to set a flag and call `<prefix>_gpu_reinit()` at the start of
   `compute()` when it is set.  A style without that path sets `reinitflag = 0`, so
   that `fix adapt` stops with an error instead of changing only the host copy.
7. **Special bond factors apply to the force and to the energy.**  Check which
   convention the cloned kernel uses (`sp_lj[...]` for cut Coulomb, `1.0 - sp_lj[...]`
   for the erfc based styles) before changing a formula.
8. `init_message()` gets the exact style name, without `/gpu`.
9. Do not copy a defect of the CPU style into the kernel: fix the CPU style and its
   `/omp` and `/kk` variants in a separate commit.

## Testing

- The YAML test of the CPU style runs the `gpu` case by itself when the variant
  exists.  A missing variant or a `skip_tests` entry passes silently: confirm from
  the test output that the case ran.
- Validate in double precision first (`-D GPU_PREC=double`), then in mixed precision
  (the default).  Long-range Coulomb on the device is always the analytic erfc
  polynomial, so the `*_table*` YAML files of such a style list `gpu` in `skip_tests`.
- A HIP build compiles all kernels without running them (`cmake --build <build>
  --target gpu`).  Run new kernels on a dedicated test machine when possible: a
  faulting kernel can hang a GPU that also drives the desktop.
