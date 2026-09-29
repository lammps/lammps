.. index:: pair_style pace
.. index:: pair_style pace/kk
.. index:: pair_style pace/extrapolation
.. index:: pair_style pace/extrapolation/kk

pair_style pace command
=======================

Accelerator Variants: *pace/kk*

pair_style pace/extrapolation command
=====================================

Accelerator Variants: *pace/extrapolation/kk*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style pace ... keyword values ...

* one or more keyword/value pairs may be appended

  .. parsed-literal::

     keyword = *product* or *recursive* or *chunksize* or *neigh*
       *product* = use product algorithm for basis functions
       *recursive* = use recursive algorithm for basis functions
       *chunksize* value = number of atoms in each pass
       *neigh* value = *auto* or *shared* or *global*
         *auto* = select the team scratch memory level automatically (default)
         *shared* = force on-chip (level 0) shared memory scratch
         *global* = force global (level 1) memory scratch

.. code-block:: LAMMPS

   pair_style pace/extrapolation

Examples
""""""""

.. code-block:: LAMMPS

   pair_style pace
   pair_style pace product chunksize 2048
   pair_coeff * * Cu-PBE-core-rep.ace Cu

   pair_style pace
   pair_coeff * * Cu.yaml Cu

   pair_style pace/extrapolation
   pair_coeff * * Cu.yaml Cu.asi Cu

Description
"""""""""""

Pair style *pace* computes interactions using the Atomic Cluster
Expansion (ACE), which is a general expansion of the atomic energy in
multi-body basis functions. :ref:`(Drautz19) <Drautz20191>`.  The *pace*
pair style provides an efficient implementation that is described in
this paper :ref:`(Lysogorskiy21) <Lysogorskiy20211>`.

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

* ACE coefficient file (.yaml or .yace/.ace format)
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

.. versionadded:: 2Sep2026

The keyword *neigh* is only applicable when using the pair styles *pace*
and *pace/extrapolation* with the KOKKOS package on GPUs and is ignored
otherwise, so that the same input file can be used with and without the
KOKKOS package.  This keyword controls which level of Kokkos team scratch
memory is used to build the short neighbor list.  Level 0 is fast on-chip
shared memory, but it is a limited resource that can be exceeded when
atoms have many neighbors and/or when there are many atomic species, which
would otherwise abort the run with an error such as "Requested too much
scratch memory on level 0".  Level 1 is (much larger) global memory,
which avoids the limit at the cost of slower access.

With the default value *auto*, the pair style queries the amount of
shared memory available on the device (rather than assuming a fixed
value such as 48 KiB, so that larger limits available in newer versions
of the Kokkos library are used automatically) and transparently falls
back to level 1 when the request does not fit into level 0, printing a
warning the first time this happens.  The value *shared* forces the use
of level 0 (on-chip) scratch memory, and *global* forces the use of
level 1 (global) scratch memory; the latter can be used to silence the
fallback warning or to force global memory when the automatic heuristic
is too conservative.

.. versionchanged:: 2Sep2026

When *pace/kk* is used with the *product* keyword on a CPU back end
(KOKKOS built with the OpenMP or Serial back end), the calculation now
runs in the KOKKOS kernels and uses all available threads.  Previously
*pace/kk* stopped with an error when more than one thread was used on a
CPU.  With the *recursive* keyword, or for potentials whose correlation
order exceeds what the CPU kernels support, *pace/kk* prints a warning
and falls back to the non-accelerated evaluator, which is single
threaded; use the *product* keyword to get the threaded calculation.

Extrapolation grade
"""""""""""""""""""

Calculation of extrapolation grade in PACE is implemented in `pair_style
pace/extrapolation`.  It is based on the MaxVol algorithm similar to
Moment Tensor Potential (MTP) by Shapeev et al.  and is described in
:ref:`(Lysogorskiy23) <Lysogorskiy2023>`.  In order to compute
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

Core repulsion
"""""""""""""""""""
The ACE potential can be configured to initiate core-repulsion from an inner cutoff,
seamlessly transitioning from ACE to ZBL. The core repulsion factor can be accessed
as a per-atom quantity, as demonstrated in the example below:

.. code-block:: LAMMPS

    pair_style  pace
    pair_coeff  * * CuNi.yaml Cu Ni

    fix pace_corerep all pair 1 pace corerep 1

In this case, per-atom `f_pace_corerep` quantities represent the fraction of ZBL
core-repulsion for each atom.

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

**(Drautz19)** Drautz, Phys Rev B, 99, 014104 (2019).

.. _Lysogorskiy20211:

**(Lysogorskiy21)** Lysogorskiy, van der Oord, Bochkarev, Menon, Rinaldi, Hammerschmidt, Mrovec, Thompson, Csanyi, Ortner, Drautz, npj Comp Mat, 7, 97 (2021).

.. _Lysogorskiy2023:

**(Lysogorskiy23)** Lysogorskiy, Bochkarev, Mrovec, Drautz, Phys Rev Mater, 7, 043801 (2023) / arXiv:2212.08716 (2022).
