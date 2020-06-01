.. index:: dihedral_style charmm

dihedral_style charmm command
=============================

dihedral_style charmm/intel command
===================================

dihedral_style charmm/kk command
================================

dihedral_style charmm/omp command
=================================

dihedral_style charmmfsw command
================================

Syntax
""""""

.. code-block:: LAMMPS

   dihedral_style style

* style = *charmm* or *charmmfsw*

Examples
""""""""

.. code-block:: LAMMPS

   dihedral_style charmm
   dihedral_style charmmfsw
   dihedral_coeff  1 0.2 1 180 1.0
   dihedral_coeff  2 1.8 1   0 1.0
   dihedral_coeff  1 3.1 2 180 0.5

Description
"""""""""""

The *charmm* and *charmmfsw* dihedral styles use the potential

.. math::

   E = K [ 1 + \cos (n \phi - d) ]

See :ref:`(MacKerell) <dihedral-MacKerell>` for a description of the CHARMM
force field.  This dihedral style can also be used for the AMBER force
field (see comment on weighting factors below).  See
:ref:`(Cornell) <dihedral-Cornell>` for a description of the AMBER force
field.

.. note::

   The newer *charmmfsw* style was released in March 2017.  We
   recommend it be used instead of the older *charmm* style when running
   a simulation with the CHARMM force field, either with long-range
   Coulombics or a Coulombic cutoff, via the :doc:`pair_style lj/charmmfsw/coul/long <pair_charmm>` and :doc:`pair_style lj/charmmfsw/coul/charmmfsh <pair_charmm>` commands respectively.
   Otherwise the older *charmm* style is fine to use.  See the discussion
   below and more details on the :doc:`pair_style charmm <pair_charmm>` doc
   page.

The following coefficients must be defined for each dihedral type via the
:doc:`dihedral_coeff <dihedral_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* :math:`K` (energy)
* :math:`n` (integer >= 0)
* :math:`d` (integer value of degrees)
* weighting factor (1.0, 0.5, or 0.0)

The weighting factor is required to correct for double counting
pairwise non-bonded Lennard-Jones interactions in cyclic systems or
when using the CHARMM dihedral style with non-CHARMM force fields.
With the CHARMM dihedral style, interactions between the 1st and 4th
atoms in a dihedral are skipped during the normal non-bonded force
computation and instead evaluated as part of the dihedral using
special epsilon and sigma values specified with the
:doc:`pair_coeff <pair_charmm>` command of pair styles that contain
"lj/charmm" (e.g. :doc:`pair_style lj/charmm/coul/long <pair_charmm>`)
In 6-membered rings, the same 1-4 interaction would be computed twice
(once for the clockwise 1-4 pair in dihedral 1-2-3-4 and once in the
counterclockwise dihedral 1-6-5-4) and thus the weighting factor has
to be 0.5 in this case.  In 4-membered or 5-membered rings, the 1-4
dihedral also is counted as a 1-2 or 1-3 interaction when going around
the ring in the opposite direction and thus the weighting factor is
0.0, as the 1-2 and 1-3 exclusions take precedence.

Note that this dihedral weighting factor is unrelated to the scaling
factor specified by the :doc:`special bonds <special_bonds>` command
which applies to all 1-4 interactions in the system.  For CHARMM force
fields, the special_bonds 1-4 interaction scaling factor should be set
to 0.0. Since the corresponding 1-4 non-bonded interactions are
computed with the dihedral.  This means that if any of the weighting
factors defined as dihedral coefficients (4th coeff above) are
non-zero, then you must use a pair style with "lj/charmm" and set the
special_bonds 1-4 scaling factor to 0.0 (which is the
default). Otherwise 1-4 non-bonded interactions in dihedrals will be
computed twice.

For simulations using the CHARMM force field with a Coulombic cutoff,
the difference between the *charmm* and *charmmfsw* styles is in the
computation of the 1-4 non-bond interactions, though only if the
distance between the two atoms is within the switching region of the
pairwise potential defined by the corresponding CHARMM pair style,
i.e. within the outer cutoff specified for the pair style.  The
*charmmfsw* style should only be used when using the corresponding
:doc:`pair_style lj/charmmfsw/coul/charmmfsw <pair_charmm>` or
:doc:`pair_style lj/charmmfsw/coul/long <pair_charmm>` commands.  Use
the *charmm* style with the older :doc:`pair_style <pair_charmm>`
commands that have just "charmm" in their style name.  See the
discussion on the :doc:`CHARMM pair_style <pair_charmm>` doc page for
details.

Note that for AMBER force fields, which use pair styles with "lj/cut",
the special_bonds 1-4 scaling factor should be set to the AMBER
defaults (1/2 and 5/6) and all the dihedral weighting factors (4th
coeff above) must be set to 0.0. In this case, you can use any pair
style you wish, since the dihedral does not need any Lennard-Jones
parameter information and will not compute any 1-4 non-bonded
interactions.  Likewise the *charmm* or *charmmfsw* styles are
identical in this case since no 1-4 non-bonded interactions are
computed.

----------

Styles with a *gpu*\ , *intel*\ , *kk*\ , *omp*\ , or *opt* suffix are
functionally the same as the corresponding style without the suffix.
They have been optimized to run faster, depending on your available
hardware, as discussed on the :doc:`Speed packages <Speed_packages>` doc
page.  The accelerated styles take the same arguments and should
produce the same results, except for round-off and precision issues.

These accelerated styles are part of the GPU, USER-INTEL, KOKKOS,
USER-OMP and OPT packages, respectively.  They are only enabled if
LAMMPS was built with those packages.  See the :doc:`Build package <Build_package>` doc page for more info.

You can specify the accelerated styles explicitly in your input script
by including their suffix, or you can use the :doc:`-suffix command-line switch <Run_options>` when you invoke LAMMPS, or you can use the
:doc:`suffix <suffix>` command in your input script.

See the :doc:`Speed packages <Speed_packages>` doc page for more
instructions on how to use the accelerated styles effectively.

----------

Restrictions
""""""""""""

When using run_style :doc:`respa <run_style>`, these dihedral styles
must be assigned to the same r-RESPA level as *pair* or *outer*\ .

When used in combination with CHARMM pair styles, the 1-4
:doc:`special_bonds <special_bonds>` scaling factors must be set to 0.0.
Otherwise non-bonded contributions for these 1-4 pairs will be
computed multiple times.

These dihedral styles can only be used if LAMMPS was built with the
MOLECULE package.  See the :doc:`Build package <Build_package>` doc page
for more info.

Related commands
""""""""""""""""

:doc:`dihedral_coeff <dihedral_coeff>`

**Default:** none

----------

.. _dihedral-Cornell:

**(Cornell)** Cornell, Cieplak, Bayly, Gould, Merz, Ferguson,
Spellmeyer, Fox, Caldwell, Kollman, JACS 117, 5179-5197 (1995).

.. _dihedral-MacKerell:

**(MacKerell)** MacKerell, Bashford, Bellott, Dunbrack, Evanseck, Field,
Fischer, Gao, Guo, Ha, et al, J Phys Chem B, 102, 3586 (1998).
