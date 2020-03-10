.. index:: fix poems

fix poems command
=================

Syntax:

.. parsed-literal::

   fix ID group-ID poems keyword values

* ID, group-ID are documented in :doc:`fix <fix>` command
* poems = style name of this fix command
* keyword = *group* or *file* or *molecule*

  .. parsed-literal::

       *group* values = list of group IDs
       *molecule* values = none
       *file* values = filename

Examples
""""""""

.. parsed-literal::

   fix 3 fluid poems group clump1 clump2 clump3
   fix 3 fluid poems file cluster.list

Description
"""""""""""

Treats one or more sets of atoms as coupled rigid bodies.  This means
that each timestep the total force and torque on each rigid body is
computed and the coordinates and velocities of the atoms are updated
so that the collection of bodies move as a coupled set.  This can be
useful for treating a large biomolecule as a collection of connected,
coarse-grained particles.

The coupling, associated motion constraints, and time integration is
performed by the software package `Parallelizable Open source Efficient Multibody Software (POEMS) <poems_>`_ which computes the
constrained rigid-body motion of articulated (jointed) multibody
systems :ref:`(Anderson) <Anderson>`.  POEMS was written and is distributed
by Prof Kurt Anderson, his graduate student Rudranarayan Mukherjee,
and other members of his group at Rensselaer Polytechnic Institute
(RPI).  Rudranarayan developed the LAMMPS/POEMS interface.  For
copyright information on POEMS and other details, please refer to the
documents in the poems directory distributed with LAMMPS.

.. _poems: http://www.rpi.edu/~anderk5/lab

This fix updates the positions and velocities of the rigid atoms with
a constant-energy time integration, so you should not update the same
atoms via other fixes (e.g. nve, nvt, npt, temp/rescale, langevin).

Each body must have a non-degenerate inertia tensor, which means if
must contain at least 3 non-collinear atoms.  Which atoms are in which
bodies can be defined via several options.

For option *group*\ , each of the listed groups is treated as a rigid
body.  Note that only atoms that are also in the fix group are
included in each rigid body.

For option *molecule*\ , each set of atoms in the group with a different
molecule ID is treated as a rigid body.

For option *file*\ , sets of atoms are read from the specified file and
each set is treated as a rigid body.  Each line of the file specifies
a rigid body in the following format:

ID type atom1-ID atom2-ID atom3-ID ...

ID as an integer from 1 to M (the number of rigid bodies).  Type is
any integer; it is not used by the fix poems command.  The remaining
arguments are IDs of atoms in the rigid body, each typically from 1 to
N (the number of atoms in the system).  Only atoms that are also in
the fix group are included in each rigid body.  Blank lines and lines
that begin with '#' are skipped.

A connection between a pair of rigid bodies is inferred if one atom is
common to both bodies.  The POEMS solver treats that atom as a
spherical joint with 3 degrees of freedom.  Currently, a collection of
bodies can only be connected by joints as a linear chain.  The entire
collection of rigid bodies can represent one or more chains.  Other
connection topologies (tree, ring) are not allowed, but will be added
later.  Note that if no joints exist, it is more efficient to use the
:doc:`fix rigid <fix_rigid>` command to simulate the system.

When the poems fix is defined, it will print out statistics on the
total # of clusters, bodies, joints, atoms involved.  A cluster in
this context means a set of rigid bodies connected by joints.

For computational efficiency, you should turn off pairwise and bond
interactions within each rigid body, as they no longer contribute to
the motion.  The "neigh\_modify exclude" and "delete\_bonds" commands
can be used to do this if each rigid body is a group.

For computational efficiency, you should only define one fix poems
which includes all the desired rigid bodies.  LAMMPS will allow
multiple poems fixes to be defined, but it is more expensive.

The degrees-of-freedom removed by coupled rigid bodies are accounted
for in temperature and pressure computations.  Similarly, the rigid
body contribution to the pressure virial is also accounted for.  The
latter is only correct if forces within the bodies have been turned
off, and there is only a single fix poems defined.

**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.

The :doc:`fix_modify <fix_modify>` *bodyforces* option is supported by
this fix style to set whether per-body forces and torques are computed
early or late in a timestep, i.e. at the post-force stage or at the
final-integrate stage, respectively.

No global or per-atom quantities are stored by this fix for access by
various :doc:`output commands <Howto_output>`.  No parameter of this fix
can be used with the *start/stop* keywords of the :doc:`run <run>`
command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix is part of the POEMS package.  It is only enabled if LAMMPS
was built with that package, which also requires the POEMS library be
built and linked with LAMMPS.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`fix rigid <fix_rigid>`, :doc:`delete_bonds <delete_bonds>`,
:doc:`neigh_modify <neigh_modify>` exclude

**Default:** none

----------

.. _Anderson:

**(Anderson)** Anderson, Mukherjee, Critchley, Ziegler, and Lipton
"POEMS: Parallelizable Open-source Efficient Multibody Software ",
Engineering With Computers (2006). (`link to paper <http://dx.doi.org/10.1007/s00366-006-0026-x>`_)
