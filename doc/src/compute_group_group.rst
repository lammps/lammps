.. index:: compute group/group

compute group/group command
===========================

Syntax
""""""

.. parsed-literal::

   compute ID group-ID group/group group2-ID keyword value ...

* ID, group-ID are documented in :doc:`compute <compute>` command
* group/group = style name of this compute command
* group2-ID = group ID of second (or same) group
* zero or more keyword/value pairs may be appended
* keyword = *pair* or *kspace* or *boundary* or *molecule*

  .. parsed-literal::

       *pair* value = *yes* or *no*
       *kspace* value = *yes* or *no*
       *boundary* value = *yes* or *no*
       *molecule* value = *off* or *inter* or *intra*

Examples
""""""""

.. parsed-literal::

   compute 1 lower group/group upper
   compute 1 lower group/group upper kspace yes
   compute mine fluid group/group wall

Description
"""""""""""

Define a computation that calculates the total energy and force
interaction between two groups of atoms: the compute group and the
specified group2.  The two groups can be the same.

If the *pair* keyword is set to *yes*\ , which is the default, then the
the interaction energy will include a pair component which is defined
as the pairwise energy between all pairs of atoms where one atom in
the pair is in the first group and the other is in the second group.
Likewise, the interaction force calculated by this compute will
include the force on the compute group atoms due to pairwise
interactions with atoms in the specified group2.

.. note::

   The energies computed by the *pair* keyword do not include tail
   corrections, even if they are enabled via the
   :doc:`pair_modify <pair_modify>` command.

If the *molecule* keyword is set to *inter* or *intra* than an
additional check is made based on the molecule IDs of the two atoms in
each pair before including their pairwise interaction energy and
force.  For the *inter* setting, the two atoms must be in different
molecules.  For the *intra* setting, the two atoms must be in the same
molecule.

If the *kspace* keyword is set to *yes*\ , which is not the default, and
if a :doc:`kspace_style <kspace_style>` is defined, then the interaction
energy will include a Kspace component which is the long-range
Coulombic energy between all the atoms in the first group and all the
atoms in the 2nd group.  Likewise, the interaction force calculated by
this compute will include the force on the compute group atoms due to
long-range Coulombic interactions with atoms in the specified group2.

Normally the long-range Coulombic energy converges only when the net
charge of the unit cell is zero.  However, one can assume the net
charge of the system is neutralized by a uniform background plasma,
and a correction to the system energy can be applied to reduce
artifacts. For more information see :ref:`(Bogusz) <Bogusz>`.  If the
*boundary* keyword is set to *yes*\ , which is the default, and *kspace*
contributions are included, then this energy correction term will be
added to the total group-group energy.  This correction term does not
affect the force calculation and will be zero if one or both of the
groups are charge neutral.  This energy correction term is the same as
that included in the regular Ewald and PPPM routines.

.. note::

   The *molecule* setting only affects the group/group
   contributions calculated by the *pair* keyword.  It does not affect
   the group/group contributions calculated by the *kspace* keyword.

This compute does not calculate any bond or angle or dihedral or
improper interactions between atoms in the two groups.

----------

The pairwise contributions to the group-group interactions are
calculated by looping over a neighbor list.  The Kspace contribution
to the group-group interactions require essentially the same amount of
work (FFTs, Ewald summation) as computing long-range forces for the
entire system.  Thus it can be costly to invoke this compute too
frequently.

.. note::

   If you have a bonded system, then the settings of
   :doc:`special_bonds <special_bonds>` command can remove pairwise
   interactions between atoms in the same bond, angle, or dihedral.  This
   is the default setting for the :doc:`special_bonds <special_bonds>`
   command, and means those pairwise interactions do not appear in the
   neighbor list.  Because this compute uses a neighbor list, it also
   means those pairs will not be included in the group/group interaction.
   This does not apply when using long-range coulomb interactions
   (\ *coul/long*\ , *coul/msm*\ , *coul/wolf* or similar.  One way to get
   around this would be to set special\_bond scaling factors to very tiny
   numbers that are not exactly zero (e.g. 1.0e-50). Another workaround
   is to write a dump file, and use the :doc:`rerun <rerun>` command to
   compute the group/group interactions for snapshots in the dump file.
   The rerun script can use a :doc:`special_bonds <special_bonds>` command
   that includes all pairs in the neighbor list.

If you desire a breakdown of the interactions into a pairwise and
Kspace component, simply invoke the compute twice with the appropriate
yes/no settings for the *pair* and *kspace* keywords.  This is no more
costly than using a single compute with both keywords set to *yes*\ .
The individual contributions can be summed in a
:doc:`variable <variable>` if desired.

This `document <PDF/kspace.pdf>`_ describes how the long-range
group-group calculations are performed.

----------

**Output info:**

This compute calculates a global scalar (the energy) and a global
vector of length 3 (force), which can be accessed by indices 1-3.
These values can be used by any command that uses global scalar or
vector values from a compute as input.  See the :doc:`Howto output <Howto_output>` doc page for an overview of LAMMPS output
options.

Both the scalar and vector values calculated by this compute are
"extensive".  The scalar value will be in energy :doc:`units <units>`.
The vector values will be in force :doc:`units <units>`.

Restrictions
""""""""""""

Not all pair styles can be evaluated in a pairwise mode as required by
this compute.  For example, 3-body and other many-body potentials,
such as :doc:`Tersoff <pair_tersoff>` and
:doc:`Stillinger-Weber <pair_sw>` cannot be used.  :doc:`EAM <pair_eam>`
potentials will re-use previously computed embedding term contributions,
so the computed pairwise forces and energies are based on the whole
system and not valid if particles have been moved since.

Not all :doc:`Kspace styles <kspace_style>` support the calculation of
group/group interactions. The regular *ewald* and *pppm* styles do.

**Related commands:** none

Default
"""""""

The option defaults are pair = yes, kspace = no, boundary = yes,
molecule = off.

----------

.. _Bogusz:

Bogusz et al, J Chem Phys, 108, 7070 (1998)
