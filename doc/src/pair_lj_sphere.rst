.. index:: pair_style lj/sphere
.. index:: pair_style lj/sphere/omp

pair_style lj/sphere command
============================

Accelerator Variant: *lj/sphere/omp*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style style args

* style = *lj/sphere*
* args = list of arguments for a particular style

.. parsed-literal::

     *lj/sphere* args = cutoff ratio
       cutoff = global cutoff ratio for Lennard Jones interactions (unitless)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style lj/sphere 2.5
   pair_coeff * * 1.0
   pair_coeff 1 1 1.1 2.8

Description
"""""""""""

The *lj/sphere* style compute the standard 12/6 Lennard-Jones potential,
given by

.. math::

   E = 4 \epsilon \left[ \left(\frac{\sigma_{ij}}{r}\right)^{12} -
       \left(\frac{\sigma_{ij}}{r}\right)^6 \right]  \qquad  r < r_c * \sigma_{ij}

:math:`r_c` is the cutoff ratio.

This is the same potential function as used by the :doc:`lj/cut
<pair_lj>` pair style, but the :math:`\sigma_{ij}` parameter is not set
as a per-type parameter via the :doc:`pair_coeff command <pair_coeff>`,
but taken from the per-atom diameter attribute of :doc:`atom_style
sphere <atom_style>`.  The individual value of :math:`\sigma_{ij}` is
computed from the diameters of the two atoms using the mixing rule for
pair coefficients as set by the :doc:`pair_modify mix <pair_modify>`
command (defaults to geometric mixing).  The cutoff is not specified as
a distance, but as ratio that is internally multiplied with
:math:`\sigma_{ij}` to obtain the actual cutoff for each pair of atoms.

Note that :math:`\sigma_{ij}` is defined in the LJ formula above as the
zero-crossing distance for the potential, *not* as the energy minimum which
is at :math:`2^{\frac{1}{6}} \sigma_{ij}`.

.. admonition:: Notes on cutoffs, neighbor lists, and efficiency
   :class: note

   Because the cutoff in this pair style depends on the diameter of the
   atoms, this influences how cutoffs are used to build neighbor lists
   and how effective those neighbor lists are at avoiding computation of
   pairwise distances between non-interacting atoms.  This pair style
   uses a conventional neighbor list construction similar to :doc:`pair
   style lj/cut <pair_lj>`, where the cutoffs are typically rather
   similar.  LAMMPS will determine the largest cutoff and use this value
   for building the neighbor lists.  This can be inefficient, if the
   differences between per-type cutoffs are large.  The command
   :doc:`neighbor multi <neighbor>` enables a modified neighbor list
   algorithm, that uses different size bins for atom types with
   different cutoffs.  It constructs adapted neighbor lists based
   on the per-type cutoffs to improve efficiency.

   If atom diameters vary largely when using pair style *lj/sphere*,
   the cutoffs computed from atom diameter and cutoff ratio with vary
   largely as well and :doc:`neighbor bin <neighbor>` based neighbor
   lists using only the largest cutoff be similarly inefficient as
   pair style *lj/cut* with largely varying per-type cutoffs.  However, the
   multi-cutoff neighbor list algorithm can only be applied when atoms
   with different cutoffs have different atom types.  Thus atoms with
   different ranges of diameters need to have different atom types to
   benefit from the  multi-cutoff neighbor lists.  LAMMPS will determine
   the largest diameter for each atom type, multiply it with the cutoff
   ratio, and use this cutoff in the same way as the per-type cutoffs in
   :doc:`pair style lj/cut <pair_lj>`

   Example input to group small and large atoms by type:

   .. code-block:: c++

      units           lj
      atom_style      sphere
      lattice         fcc 0.8442
      region          box block 0 10 0 10 0 10
      create_box      2 box
      create_atoms    1 box

      # create atoms with random diameters from bimodal distribution
      variable switch atom random(0.0,1.0,345634)
      variable diam atom (v_switch<0.75)*normal(0.4,0.075,325)+(v_switch>=0.7)*normal(1.2,0.2,453)
      set group all diameter v_diam

      # assign type 2 to atoms with diameter > 0.5
      variable large atom 2.0*radius>0.5
      group large variable large
      set group largea type 2

      pair_style      lj/sphere 2.5
      pair_coeff      * * 1.0

      neighbor 0.3 multi



Coefficients
""""""""""""

The following coefficients must be defined for each pair of atoms types via the
:doc:`pair_coeff <pair_coeff>` command as in the examples above, or in the data
file or restart files read by the :doc:`read_data <read_data>` or
:doc:`read_restart <read_restart>` commands, or by mixing as described below:

* :math:`\epsilon` (energy units)
* LJ cutoff ratio (unitless) (optional)

The last coefficient is optional.  If not specified, the global LJ
cutoff ratio specified in the :doc:`pair_style command <pair_style>` is
used.

If a repulsive only LJ interaction is desired, the coefficient for the cutoff
ratio should be set to the minimum of the LJ potential using ``$(2.0^(1.0/6.0))``

----------

.. include:: accel_styles.rst

----------

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For atom type pairs I,J and I != J, the epsilon coefficients and cutoff
ratio for the *lj/sphere* pair style can be mixed.  The default mixing
style is *geometric*.  See the :doc:`pair_modify command <pair_modify>`
for details.

The *lj/sphere* pair style supports the :doc:`pair_modify shift <pair_modify>`
option for the energy of the Lennard-Jones portion of the pair interaction.

The *lj/sphere* pair style does *not* support the :doc:`pair_modify
<pair_modify>` tail option for adding a long-range tail corrections to
the energy and pressure.

The *lj/sphere* pair style writes its information to :doc:`binary
restart files <restart>`, so pair_style and pair_coeff commands do not
need to be specified in an input script that reads a restart file.

This pair style can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  It does *not* support the
*inner*, *middle*, *outer* keywords.

----------

Restrictions
""""""""""""

The *lj/sphere* pair style is only enabled if LAMMPS was built with the
EXTRA-PAIR package.  See the :doc:`Build package <Build_package>` page
for more info.

The *lj/sphere* pair style does not support the *sixthpower* mixing rule.

----------

Related commands
""""""""""""""""

* :doc:`pair_coeff <pair_coeff>`
* :doc:`pair_style lj/cut <pair_lj>`

Default
"""""""

none
