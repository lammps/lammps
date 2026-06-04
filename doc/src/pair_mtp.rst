.. index:: pair_style mtp
.. index:: pair_style mtp/kk
.. index:: pair_style mtp/extrapolation
.. index:: pair_style mtp/extrapolation/kk

pair_style mtp command
=======================

Accelerator Variants: *mtp/kk*

pair_style mtp/extrapolation command
=====================================

Accelerator Variants: *mtp/extrapolation/kk*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style mtp ... keyword values ...

.. code-block:: LAMMPS

   pair_style mtp/extrapolation ... keyword values ...

* one or more keyword/value pairs may be appended

  .. parsed-literal::

     keyword = *chunksize*
       *chunksize* value = number of atoms in each pass

Examples
""""""""

.. code-block:: LAMMPS

   pair_style mtp
   pair_style mtp chunksize 2048
   pair_coeff * * SiO.mtp 0 1

   pair_style mtp/extrapolation
   pair_style mtp chunksize 2048
   pair_coeff * * SiO.almtp 0 1

Description
"""""""""""

Pair style *mtp* computes interactions using the Moment Tensor Potentials (MTP), which is a general expansion of the atomic energy in scalar contractions of moment tensors. :ref:`(Shapeev16) <Shapeev2016>`.  The *mtp* pair style provides an efficient implementation that is described in
this paper :ref:`(Meng26) <Meng2026>`. 

In the MTP, the total energy is decomposed into a sum over atomic
energies. The energy of atom *i* is expressed as a linear function of scalar contractions of moment tensors.

Only a single pair_coeff command is used with the *mtp* style which
specifies an MTP potential file followed by N additional arguments
specifying the mapping of MTP elements to LAMMPS atom types, where N is
the number of LAMMPS atom types:

* MTP potential file (.mtp or .almtp MLIP-2/MLIP-3 format)
* N element types = mapping of MTP element types (0,1,...) to atom types

Note that unlike for some other potentials, a single cutoff is used for all type pair which is not set in the pair_style or
pair_coeff command but rather, specified in the MTP file.

The keyword *chunksize* is only applicable when using the pair style
*mtp* with the KOKKOS package on GPUs and is ignored otherwise.  This
keyword controls the number of atoms in each pass and is used to avoid running out of memory. For example if there are 8192 atoms in the simulation and the *chunksize* is set to 4096, the MTP calculation will be broken up into two passes (running on a single GPU).

Extrapolation grade
"""""""""""""""""""

Calculation of extrapolation grade in MTP is implemented in `pair_style
mtp/extrapolation`.  It adapts the MaxVol algorithm and is described in
:ref:`(Podryabinkin17) <Podryabinkin2017>`.  In order to compute
extrapolation grade one needs to provide an MTP potential file in the MLIP-3 format :ref:`(Podryabinkin23) <Podryabinkin2023>` with active learning information included (.almtp).

Calculation of extrapolation grades requires additional derivatives and a matrix-vector multiplication for each atom, and is slower than the usual `pair_style mtp`, therefore it is *not* computed by default.
Extrapolation grade calculation is involved by `fix pair`, which
requests to compute `gamma`, as shown in example below:

.. code-block:: LAMMPS

    pair_style  mtp/extrapolation
    pair_coeff * * SiO.mtp 0 1

    fix mtp_gamma all pair 10 mtp/extrapolation gamma 1

    # compute max_mtp_gamma all reduce max f_mtp_gamma
    compute max_mtp_gamma all pair mtp/extrapolation
    variable dump_skip equal "c_max_mtp_gamma < 2"

    dump mtp_dump all custom 20 extrapolative_structures.dump id type x y z f_mtp_gamma
    dump_modify mtp_dump skip v_dump_skip

    variable max_mtp_gamma equal c_max_mtp_gamma
    fix extreme_extrapolation all halt 10 v_max_mtp_gamma > 10

Here extrapolation grade gamma is computed every 10 steps and is stored
in `f_mtp_gamma` per-atom variable.  The largest value of extrapolation
grade among all atoms in a structure is exposed to the `c_max_mtp_gamma`
variable.  Only if this value exceeds extrapolation threshold 2, then
the structure will be dumped into `extrapolative_structures.dump` file,
but not more often than every 20 steps. If this value exceeds 10, the simulation
is terminated.

On all other steps `pair_style mtp` will be used.

The use of pair style *mtp/extrapolation* is not recommended with 
`pair_style hybrid/overlay` since the extrapolation calculation cannot 
consider contributions from other pair styles.

----------

See the :doc:`pair_coeff <pair_coeff>` page for alternate ways
to specify the path for the MTP potential file.

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

This pair style is part of the ML-MTP package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`pair_style snap  <pair_snap>`,
:doc:`pair_style pace  <pair_pace>`,
:doc:`fix pair  <fix_pair>`
:doc:`compute pair  <compute_pair>`

Default
"""""""

chunksize = 4096,

.. _Shapeev2016:

**(Shapeev16)** Shapeev, Multiscale Model Simul, 14, 1153--1173 (2016).

.. _Meng2026:

**(Meng26)** Meng, Zongo, Torres, Maxwell, Grant, Beland, SoftwareX, 33, 102524 (2026).

.. _Podryabinkin2017:

**(Podryabinkin17)** Podryabinkin, Shapeev, Comput Mater Sci, 140, 171--180 (2017).

.. _Podryabinkin2023:

**(Podryabinkin23)** Podryabinkin, Garifullin, Shapeev, Novikov, J Chem Phys, 159, 8 (2023).
