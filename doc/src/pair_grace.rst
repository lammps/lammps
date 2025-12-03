.. index:: pair_style grace
.. index:: pair_style grace/fs

pair_style grace command
========================

pair_style grace/fs command
===========================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style grace keyword values ...

* one or more keyword/value pairs may be appended

  .. parsed-literal::

     keyword = *padding* or *pad_verbose* or *pair_forces* or *max_number_of_reduction* or *reduce_padding*
       *padding* value = fraction of neighbors to pad (default 0.01)
       *pad_verbose* = print messages when padding triggers recompilation
       *pair_forces* = compute pairwise forces (required for virials and MPI > 1)
       *max_number_of_reduction* value = maximum number of recompilations during padding reduction
       *reduce_padding* value = fraction to reduce padding by

.. code-block:: LAMMPS

   pair_style grace/fs keyword values ...

* one or more keyword/value pairs may be appended

  .. parsed-literal::

     keyword = *extrapolation*
       *extrapolation* = compute extrapolation grade (requires .asi file in pair_coeff)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style grace
   pair_coeff * * /path/to/saved_model Al Li

   pair_style grace padding 0.05 pad_verbose pair_forces
   pair_coeff * * /path/to/grace_cache/AlLi_model Al Li

   pair_style grace/fs
   pair_coeff * * FS_model.yaml Mo Nb Ta W

   pair_style grace/fs extrapolation
   pair_coeff * * FS_model.yaml FS_model.asi Mo Nb Ta W

Description
"""""""""""

The *grace* and *grace/fs* pair styles compute interactions using the
Graph Atomic Cluster Expansions (GRACE) framework :ref:`(Bochkarev24) <Bochkarev20241>`,
:ref:`(Lysogorskiy25) <Lysogorskiy20251>`.

**pair_style grace**

The *grace* style utilizes the *libtensorflow*  to load and execute
GRACE models, saved in the  TensorFlow ``saved_model`` format.

Only a single pair_coeff command is used with the *grace* style which
specifies the directory containing the saved model followed by N additional
arguments specifying the mapping of GRACE model elements to LAMMPS atom types,
where N is the number of LAMMPS atom types:

* path to ``saved_model`` directory
* N element names = mapping of model elements to atom types

GRACE models in TensorFlow are Just-In-Time (JIT) compiled. This means the
first evaluation will be slower than subsequent steps. To maintain performance
if the number of neighbors changes, the style uses a padding strategy.

* **padding**: Sets the fraction of neighbors to pad (default is 0.01 or 1%). Increasing this can reduce the frequency of recompilations but increases the time and memory overhead.
* **pad_verbose**: If specified, LAMMPS will output messages whenever new padding levels trigger a recompilation. By default this is false.
* **pair_forces**: By default, the GRACE model provides total atomic forces. If *pair_forces* is enabled, the model calculates pairwise forces. This is **required** for calculating atomic virials (stress) and is automatically enforced if running on more than one MPI processor.
* **max_number_of_reduction** and **reduce_padding**: Control the heuristics for reducing the padding buffer size dynamically during the simulation. This is a type of nice-to-have optimization.

**pair_style grace/fs**

The *grace/fs* style provides a native C++ implementation (product evaluator)
for the GRACE/FS family of models. It is lightweight, does not require
TensorFlow or GPUs, and supports standard MPI parallelization efficiently.

Only a single pair_coeff command is used with the *grace/fs* style:

* GRACE/FS coefficient file (.yaml format)
* (Optional) Active Set Inverted file (.asi format) if *extrapolation* keyword is used
* N element names = mapping of elements to atom types

Extrapolation grade
"""""""""""""""""""

Calculation of extrapolation grade is implemented in `pair_style grace/fs`
via the *extrapolation* keyword. It is based on the MaxVol algorithm.
In order to compute the extrapolation grade one needs to provide:

#. GRACE/FS potential in `.yaml` format
#. Active Set Inverted (ASI) file for the corresponding potential (`.asi` format)

Calculation of extrapolation grades requires matrix-vector multiplication
for each atom and is slower than the standard evaluation. The extrapolation
grade is accessed via `fix pair`, which requests to compute `gamma`.

Example of monitoring extrapolation warnings:

.. code-block:: LAMMPS

    pair_style  grace/fs extrapolation
    pair_coeff  * * FS_model.yaml FS_model.asi Mo Nb Ta W

    # Compute gamma every 100 steps, store in f_grace_gamma
    fix grace_gamma all pair 100 grace/fs gamma 1

    compute max_grace_gamma all reduce max f_grace_gamma
    variable dump_skip equal "c_max_grace_gamma < 5"

    dump grace_dump all custom 20 extrapolative_structures.dump id type x y z f_grace_gamma
    dump_modify grace_dump skip v_dump_skip

    variable max_grace_gamma equal c_max_grace_gamma
    fix extreme_extrapolation all halt 10 v_max_grace_gamma > 25

Here extrapolation grade gamma is computed every 100 steps and is stored
in the `f_grace_gamma` per-atom variable. The largest value of extrapolation
grade among all atoms in a structure is reduced to the `c_max_grace_gamma`
variable. Only if this value exceeds extrapolation threshold 5 will the
structure be dumped.

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This pair style does not support the :doc:`pair_modify <pair_modify>`
shift, table, and tail options.

This pair style does not write its information to :doc:`binary restart
files <restart>`, since it is stored in potential files/directories. Thus,
you need to re-specify the pair_style and pair_coeff commands in an input
script that reads a restart file.

----------

Restrictions
""""""""""""

These pair styles are part of the ML-PACE package. They are only enabled if
LAMMPS was built with that package.

Both styles require `metal` units.

*pair_style grace* relies on the TensorFlow library. While GPU usage is optional,
TensorFlow is significantly less efficient when running solely on the CPU.
Currently, only 1-LAYER models can be parallelized using MPI with domain
decomposition.

*pair_style grace/fs* does not require TensorFlow.

Further read
""""""""""""""""
See `gracemaker.readthedocs.io <https://gracemaker.readthedocs.io>`_ for more details.

Related commands
""""""""""""""""

:doc:`pair_style pace  <pair_pace>`,
:doc:`fix pair  <fix_pair>`

Default
"""""""

For *grace*: padding = 0.01, pair_forces is OFF, pad_verbose is OFF, max_number_of_reduction=10, reduce_padding=0.2.

For *grace/fs*: extrapolation is OFF.

.. _Bochkarev20241:

**(Bochkarev24)** Bochkarev, Lysogorskiy, Drautz, Phys Rev X, 14, 021036 (2024).

.. _Lysogorskiy20251:

**(Lysogorskiy25)** Lysogorskiy, Bochkarev, Drautz, arXiv:2508.17936 (2025).