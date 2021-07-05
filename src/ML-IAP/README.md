This package provides a general interface to families of
machine-learning interatomic potentials (MLIAPs).  This interface
consists of a `mliap pair style` and a `mliap compute`.  The `mliap pair
style` is used when running simulations with energies and forces
calculated by an MLIAP. The interface allows separate definitions of the
interatomic potential functional form (`model`) and the geometric
quantities that characterize the atomic positions (`descriptor`). By
defining `model` and `descriptor` separately, it is possible to use many
different models with a given descriptor, or many different descriptors
with a given model. The mliap pair_style supports the following models:

- `linear`,
- `quadratic`,
- `nn` (neural networks)
- `mliappy` (general Python interface to things like PyTorch).

It currently supports two classes of descriptors, `sna` (the SNAP
descriptors) and `so3` (the SO3 descriptor).  It is straightforward to
add new descriptor and model styles.

The `compute mliap` style provides gradients of the energy, force, and
stress tensor w.r.t. model parameters.  These are useful when training
MLIAPs to match target data.  Any `model` or `descriptor` that has been
implemented for the `mliap` pair style can also be accessed by the
`mliap` compute.  In addition to the energy, force, and stress
gradients, w.r.t.  each `model` parameter, the compute also calculates
the energy, force, and stress contributions from a user-specified
reference potential.

## Generating the model files from the third-party packages
- To train the `linear` and `quadratic` models with the SNAP descriptors, see the examples in [FitSNAP](https://github.com/FitSNAP/FitSNAP).
- To train the `nn` model with the SNAP descriptors, check the examples in [PyXtal\_FF](https://github.com/qzhu2017/PyXtal_FF).

## Building the ML-IAP package with Python support

The `mliappy` energy model requires that the MLIAP package is compiled
with Python support enabled. This extension, written by Nick Lubbers
(LANL), provides a coupling to PyTorch and other Python modules. This
should be automatically enabled by default if the prerequisite software
is installed when compiling with CMake. It can be enforced during CMake
configuration by setting the variable `MLIAP_ENABLE_PYTHON=yes` or for
conventional build by adding `-DMLIAP_PYTHON` to the `LMP_INC` variable
in your makefile and running the cythonize script on the .pyx file(s)
copied to the src folder.

Using the `mliappy` energy model also requires to install the PYTHON
package and have the [cython](https://cython.org) software
installed. During configuration/compilation the cythonize script will be
used to convert the provided .pyx file(s) to C++ code.  Please do *NOT*
run the cythonize script manually in the src/ML-IAP folder. If you have
done so by accident, you need to delete the generated .cpp and .h
file(s) in the src/ML-IAP folder or there may be problems during
compilation.

More information on building LAMMPS with this package is
[here](https://docs.lammps.org/Build_extras.html#mliap).
