.. index:: pair_style pace
.. index:: pair_style pace/kk
.. index:: pair_style pace/extrapolation
.. index:: pair_style pace/extrapolation/kk

pair_style pace command
=======================

Accelerator Variants: *pace/kk*

pair_style pace/extrapolation command
=====================================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style pace ... keyword values ...

* one or more keyword/value pairs may be appended

  .. parsed-literal::

     keyword = *product* or *recursive* or *chunksize*
       *product* = use product algorithm for basis functions
       *recursive* = use recursive algorithm for basis functions
       *chunksize* value = number of atoms in each pass

.. code-block:: LAMMPS

   pair_style pace/extrapolation

Examples
""""""""

.. code-block:: LAMMPS

   pair_style pace
   pair_style pace product chunksize 2048
   pair_coeff * * Cu-PBE-core-rep.ace Cu

   pair_style pace/extrapolation
   pair_coeff * * Cu.yaml Cu.asi Cu

Description
"""""""""""

Pair style *pace* computes interactions using the Atomic Cluster
Expansion (ACE), which is a general expansion of the atomic energy in
multi-body basis functions. :ref:`(Drautz) <Drautz20191>`.  The *pace*
pair style provides an efficient implementation that is described in
this paper :ref:`(Lysogorskiy) <Lysogorskiy20211>`.

In ACE, the total energy is decomposed into a sum over atomic
energies. The energy of atom *i* is expressed as a linear or non-linear
function of one or more density functions.  By projecting the density
onto a local atomic base, the lowest order contributions to the energy
can be expressed as a set of scalar polynomials in basis function
contributions summed over neighbor atoms.

Only a single pair_coeff command is used with the *pace* style which
specifies an ACE coefficient file followed by N additional arguments
specifying the mapping of ACE elements to LAMMPS atom types, where N is
the number of LAMMPS atom types:

* ACE coefficient file
* N element names = mapping of ACE elements to atom types

Only a single pair_coeff command is used with the *pace* style which
specifies an ACE file that fully defines the potential.  Note that
unlike for other potentials, cutoffs are not set in the pair_style or
pair_coeff command; they are specified in the ACE file.

The pair_style *pace* command may be followed by the optional keyword
*product* or *recursive*, which determines which of two algorithms is
used for the calculation of basis functions and derivatives.  The
default is *recursive*.

The keyword *chunksize* is only applicable when using the pair style
*pace* with the KOKKOS package on GPUs and is ignored otherwise.  This
keyword controls the number of atoms in each pass used to compute the
atomic cluster expansion and is used to avoid running out of memory.
For example if there are 8192 atoms in the simulation and the
*chunksize* is set to 4096, the ACE calculation will be broken up into
two passes (running on a single GPU).

Extrapolation grade
"""""""""""""""""""

Calculation of extrapolation grade in PACE is implemented in `pair_style
pace/extrapolation`.  It is based on the MaxVol algorithm similar to
Moment Tensor Potential (MTP) by Shapeev et al.  and is described in
:ref:`(Lysogorskiy2) <Lysogorskiy2022>`.  In order to compute
extrapolation grade one needs to provide:

#. ACE potential in B-basis form (`.yaml` format) and
#. Active Set Inverted (ASI) file for corresponding potential (`.asi` format)

Calculation of extrapolation grades requires matrix-vector
multiplication for each atom and is slower than the usual `pair_style
pace recursive`, therefore it is *not* computed by default.
Extrapolation grade calculation is involved by `fix pair`, which
requests to compute `gamma`, as shown in example below:

.. code-block:: LAMMPS

    pair_style  pace/extrapolation
    pair_coeff  * * Cu.yaml Cu.asi Cu

    fix pace_gamma all pair 10 pace/extrapolation gamma 1

    compute max_pace_gamma all reduce max f_pace_gamma
    variable dump_skip equal "c_max_pace_gamma < 5"

    dump pace_dump all custom 20 extrapolative_structures.dump id type x y z f_pace_gamma
    dump_modify pace_dump skip v_dump_skip

    variable max_pace_gamma equal c_max_pace_gamma
    fix extreme_extrapolation all halt 10 v_max_pace_gamma > 25

Here extrapolation grade gamma is computed every 10 steps and is stored
in `f_pace_gamma` per-atom variable.  The largest value of extrapolation
grade among all atoms in a structure is reduced to `c_max_pace_gamma`
variable.  Only if this value exceeds extrapolation threshold 5, then
the structure will be dumped into `extrapolative_structures.dump` file,
but not more often than every 20 steps.

On all other steps `pair_style pace recursive` will be used.

When using the pair style *pace/extrapolation* with the KOKKOS package on GPUs
product B-basis evaluator is always used and only *linear* ASI is supported.

----------

See the :doc:`pair_coeff <pair_coeff>` page for alternate ways
to specify the path for the ACE coefficient file.

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For atom type pairs I,J and I != J, where types I and J correspond to
two different element types, mixing is performed by LAMMPS with
user-specifiable parameters as described above.  You never need to
specify a pair_coeff command with I != J arguments for this style.

This pair style does not support the :doc:`pair_modify <pair_modify>`
shift, table, and tail options.

This pair style does not write its information to :doc:`binary restart
files <restart>`, since it is stored in potential files.  Thus, you need
to re-specify the pair_style and pair_coeff commands in an input script
that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*, *middle*, *outer* keywords.

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

This pair style is part of the ML-PACE package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`pair_style snap  <pair_snap>`,
:doc:`fix pair  <fix_pair>`

Default
"""""""

recursive, chunksize = 4096,

.. _Drautz20191:

**(Drautz)** Drautz, Phys Rev B, 99, 014104 (2019).

.. _Lysogorskiy20211:

**(Lysogorskiy)** Lysogorskiy, van der Oord, Bochkarev, Menon, Rinaldi, Hammerschmidt, Mrovec, Thompson, Csanyi, Ortner, Drautz, npj Comp Mat, 7, 97 (2021).

.. _Lysogorskiy2022:

**(Lysogorskiy2022)** Lysogorskiy, Bochkarev, Mrovec, Drautz, arXiv:2212.08716 (2022).
