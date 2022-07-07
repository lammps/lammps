.. index:: pair_style pace
.. index:: pair_style pace/kk
.. index:: pair_style pace/extrapolation

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

   pair_style pace/extrapolation gamma_lower_bound gamma_upper_bound gamma_freq

* one or more arguments may be appended

  .. parsed-literal::
    *gamma_lower_bound* = minimal value of extrapolation grade considered as moderate extrapolation
    *gamma_upper_bound* = maximal value of extrapolation grade considered as moderate extrapolation
    *gamma_freq* value = frequency of computing extrapolation grade (in steps)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style pace
   pair_style pace product chunksize 2048
   pair_coeff * * Cu-PBE-core-rep.ace Cu

   pair_style pace/extrapolation 1.5 10 20
   pair_coeff * * Cu.yaml Cu.asi Cu

Description
"""""""""""

Pair style *pace* computes interactions using the Atomic Cluster
Expansion (ACE), which is a general expansion of the atomic energy in
multi-body basis functions. :ref:`(Drautz) <Drautz20191>`.
The *pace* pair style
provides an efficient implementation that
is described in this paper :ref:`(Lysogorskiy) <Lysogorskiy20211>`.

In ACE, the total energy is decomposed into a sum over
atomic energies. The energy of atom *i* is expressed as a
linear or non-linear function of one or more density functions.
By projecting the
density onto a local atomic base, the lowest order contributions
to the energy can be expressed as a set of scalar polynomials in
basis function contributions summed over neighbor atoms.

Only a single pair_coeff command is used with the *pace* style which
specifies an ACE coefficient file followed by N additional arguments
specifying the mapping of ACE elements to LAMMPS atom types,
where N is the number of LAMMPS atom types:

* ACE coefficient file
* N element names = mapping of ACE elements to atom types

Only a single pair_coeff command is used with the *pace* style which
specifies an ACE file that fully defines the potential.
Note that unlike for other potentials, cutoffs are
not set in the pair_style or pair_coeff command; they are specified in
the ACE file.

The pair_style *pace* command may be followed by the optional keyword
*product* or *recursive*, which determines which of two algorithms
is used for the calculation of basis functions and derivatives.
The default is *recursive*.

The keyword *chunksize* is only applicable when
using the pair style *pace* with the KOKKOS package on GPUs and is
ignored otherwise.  This keyword controls the number of atoms
in each pass used to compute the atomic cluster expansion and is used to
avoid running out of memory.  For example if there are 8192 atoms in the
simulation and the *chunksize* is set to 4096, the ACE
calculation will be broken up into two passes (running on a single GPU).

Extrapolation grade
"""""""""""""""""""

Calculation of extrapolation grade in PACE is implemented in `pair_style pace/extrapolation`.
It is based on the MaxVol algorithm similar to Moment Tensor Potential (MTP) by Shapeev et al.
and is described in :ref:`(Lysogorskiy2) <Lysogorskiy2022>`.
In order to compute extrapolation grade one needs to provide:

#. ACE potential in B-basis form (`.yaml` format) and
#. Active Set Inverted (ASI) file for corresponding potential (`.asi` format)

Calculation of extrapolation grades requires matrix-vector multiplication for each atom
and can be slower than the usual `pair_style pace recursive`,
therefore it make sense *not* to do it on every step.
Extrapolation grade calculation frequency is controlled by *gamma_freq* parameter.
On all other steps `pair_style pace recursive` is used.

The maximal value of *gamma* for all atoms in a structure determines the extrapolation grade for structure.
Both per-atom and per-structure extrapolation grades are accessible via `compute pace/extrapolation`:

.. code-block:: LAMMPS

   compute pace_gamma all pace/extrapolation

   # show maximal extrapolation grade per-structure
   thermo_style custom step etotal temp press c_pace_gamma

   # dump structure with per-atom extrapolation grades
   dump 1 all custom 100 my.dump id type mass x y z c_pace_gamma

If maximal extrapolation grade per-structure exceeds *gamma_lower_bound* but less than *gamma_upper_bound*,
the structure is considered as extrapolative and can be stored with `dump pace/extrapolation`:

.. code-block:: LAMMPS

   compute pace_gamma all pace/extrapolation
   dump pace all pace/extrapolation 1 extrapolation.dat id type mass x y z c_pace_gamma

Please note, that even if you provide dump frequency equal to one, dump will write structure
only if extrapolation grades are computed on current timestep *and* maximal extrapolation grade exceeds *gamma_lower_bound*.
If extrapolation grade exceeds *gamma_upper_bound*, simulation will be aborted.



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

This pair style does not write its information to :doc:`binary restart files <restart>`, since it is stored in potential files.  Thus, you
need to re-specify the pair_style and pair_coeff commands in an input
script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does not support the
*inner*, *middle*, *outer* keywords.

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

This pair style is part of the ML-PACE package.  It is only enabled if LAMMPS
was built with that package.
See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`pair_style snap  <pair_snap>`,
:doc:`compute pace/extrapolation  <compute_pace_extrapolation>`,
:doc:`dump pace/extrapolation  <dump_pace_extrapolation>`,

Default
"""""""

recursive, chunksize = 4096,
gamma_lower_bound = 1.5,
gamma_upper_bound = 10,
gamma_freq = 1

.. _Drautz20191:

**(Drautz)** Drautz, Phys Rev B, 99, 014104 (2019).

.. _Lysogorskiy20211:

**(Lysogorskiy)** Lysogorskiy, van der Oord, Bochkarev, Menon, Rinaldi, Hammerschmidt, Mrovec, Thompson, Csanyi, Ortner, Drautz, npj Comp Mat, 7, 97 (2021).

.. _Lysogorskiy2022:

**(Lysogorskiy2022)** Lysogorskiy, Bochkarev, Mrovec, Drautz, TBS (2022).
