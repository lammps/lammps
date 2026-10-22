# Tests for dump image and GRAPHICS package rendering

The `test_dump_image` driver renders a scene described in a test YAML file
with the `dump image` command and compares sampled pixel data of the rendered
image against reference data stored directly in the same YAML file.  No
reference image files are needed.  The sampling scheme was validated against
deliberately introduced rendering regressions and toolchain variations before
this harness was written (see `plans/graphics-sampling-experiment/` notes in
the corresponding pull request).

## How images are compared

Two complementary sets of reference data are stored per test:

- `row_means` / `col_means`: mean RGB values of every pixel row and column
  (full-image coverage; any localized or diffuse change shows up in some
  row or column mean).  Compared against `epsilon_projection` (default 0.5
  on the 0-255 scale).
- `sample_blocks`: mean RGB values of small pixel blocks (default 3x3) on a
  uniform grid (default stride 16).  Compared against `epsilon_blocks`
  (default 2.0).  Mainly used to localize a failure to an image position.

Serial rendering is bit-deterministic and was found to be bit-identical
across compilers and optimization levels; the tolerances leave headroom for
different libm implementations.  Rendered images legitimately depend on the
MPI rank count (pixel-level depth ties composite differently and SSAO
partitions its random stream by rank), therefore each test carries an
`nprocs` key and is skipped for any other rank count.  Tests with
`nprocs: 1` run serially; tests with larger values are registered with the
MPI launcher by CMake and exercise the parallel image compositing.

## Test YAML format

```yaml
---
lammps_version: <version of LAMMPS the reference was created with>
tags: graphics
date_generated: <date>
epsilon_projection: 0.5     # tolerance for row/column mean RGB values
epsilon_blocks: 2.0         # tolerance for sampled block mean RGB values
nprocs: 1                   # required MPI rank count
prerequisites: ! |          # styles that must be present in the build
  atom full
  fix graphics/arrows
setup_commands: ! |         # scene setup incl. the commands under test
  ...
  dump viz all image 1 ${imagefile} type type size 200 200 ...
  dump_modify viz ...
run_commands: ! |           # commands that trigger the render
  run 0
image_size: 200 200         # expected image dimensions
sampling: 3 16              # block size and grid stride in pixels
row_means: ...              # reference data, created by the generator
col_means: ...
sample_blocks: ...
```

Notes for writing scenes:

- The image file name must be given as `${imagefile}`; the driver defines
  this variable (it contains the `*` placeholder required by dump image).
  The variable `${input_dir}` points to this `tests/` folder for data files.
- When multiple frames are written, the frame with the highest timestep
  number is compared.
- Keep images small (about 200x200): the tests run in well under a second
  and the reference data stays compact.
- Verify that the option under test actually changes the image: options can
  be silently ignored in some combinations (e.g. `dump_modify acolor` has no
  effect when the dump colors atoms by element).  Render the scene with and
  without the option and compare before committing a test.
- `autobond` searches the pair neighbor list, so it draws nothing on
  molecular systems with default `special_bonds 0 0 0` settings (bonded
  pairs are excluded from the neighbor list).  Test it on systems without
  a bond topology, or set `special_bonds lj/coul 1.0 1.0 1.0`.
- `dump_modify btrans` only applies when bonds are *not* colored by atom:
  with `bond atom ...` the two bond halves follow the transparency of
  their respective atoms (`atrans`) instead.
- Continuous color maps (`amap`/`bmap`/`gmap` styles `ca`/`cf`) require the
  first entry value to be the literal `min` and the last to be `max`;
  discrete maps (`da`/`df`) require a final catch-all `min max <color>`
  entry.
- Stochastic features are fine as long as they are seeded: SSAO and the
  region `points` draw style are reproducible.

## Creating or updating reference data

Build with `-D ENABLE_TESTING=on`, then run from the build folder:

```
test_dump_image <path/to/new-test.yaml> -u                 # update in place
test_dump_image <path/to/test.yaml> -g <path/to/new.yaml>  # write new file
mpirun -np 4 test_dump_image <path/to/test.yaml> -u        # 4-rank reference
```

Start from a YAML file without the reference keys (`row_means`, `col_means`,
`sample_blocks`) or from a copy of an existing test.  The generator renders
the scene, samples the image, and writes the completed file, recording the
rank count it ran with.  On a comparison failure the test keeps the rendered
image as `<testname>.failed.ppm` for visual inspection.
