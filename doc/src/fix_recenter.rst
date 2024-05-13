.. index:: fix recenter

fix recenter command
====================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID recenter x y z keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* recenter = style name of this fix command
* x,y,z = constrain center-of-mass to these coords (distance units),         any coord can also be NULL or INIT (see below)
* zero or more keyword/value pairs may be appended
* keyword = *shift* or *units*

  .. parsed-literal::

       *shift* value = group-ID
         group-ID = group of atoms whose coords are shifted
       *units* value = *box* or *lattice* or *fraction*

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all recenter 0.0 0.5 0.0
   fix 1 all recenter INIT INIT NULL
   fix 1 all recenter INIT 0.0 0.0 units box

Description
"""""""""""

Constrain the center-of-mass position of a group of atoms by adjusting
the coordinates of the atoms every timestep.  This is simply a small
shift that does not alter the dynamics of the system or change the
relative coordinates of any pair of atoms in the group.  This can be
used to ensure the entire collection of atoms (or a portion of them)
do not drift during the simulation due to random perturbations
(e.g. :doc:`fix langevin <fix_langevin>` thermostatting).

Distance units for the x,y,z values are determined by the setting of
the *units* keyword, as discussed below.  One or more x,y,z values can
also be specified as NULL, which means exclude that dimension from
this operation.  Or it can be specified as INIT which means to
constrain the center-of-mass to its initial value at the beginning of
the run.

The center-of-mass (COM) is computed for the group specified by the
fix.  If the current COM is different than the specified x,y,z, then a
group of atoms has their coordinates shifted by the difference.  By
default the shifted group is also the group specified by the fix.  A
different group can be shifted by using the *shift* keyword.  For
example, the COM could be computed on a protein to keep it in the
center of the simulation box.  But the entire system (protein + water)
could be shifted.

If the *units* keyword is set to *box*, then the distance units of
x,y,z are defined by the :doc:`units <units>` command - e.g. Angstroms
for *real* units.  A *lattice* value means the distance units are in
lattice spacings.  The :doc:`lattice <lattice>` command must have been
previously used to define the lattice spacing.  A *fraction* value
means a fractional distance between the lo/hi box boundaries, e.g. 0.5
= middle of the box.  The default is to use lattice units.

Note that the :doc:`velocity <velocity>` command can be used to create
velocities with zero aggregate linear and/or angular momentum.

.. note::

   This fix performs its operations at the same point in the
   timestep as other time integration fixes, such as :doc:`fix nve <fix_nve>`, :doc:`fix nvt <fix_nh>`, or :doc:`fix npt <fix_nh>`.
   Thus fix recenter should normally be the last such fix specified in
   the input script, since the adjustments it makes to atom coordinates
   should come after the changes made by time integration.  LAMMPS will
   warn you if your fixes are not ordered this way.

.. note::

   If you use this fix on a small group of atoms (e.g. a molecule
   in solvent) without using the *shift* keyword to adjust the positions
   of all atoms in the system, then the results can be unpredictable.
   For example, if the molecule is pushed consistently in one direction
   by a flowing solvent, its velocity will increase.  But its coordinates
   will be re-centered, meaning it is moved back towards the force.  Thus
   over time, the velocity and effective temperature of the molecule
   could become very large, though it won't actually be moving due to the
   re-centering.  If you are thermostatting the entire system, then the
   solvent would be cooled to compensate.  A better solution for this
   simulation scenario is to use the :doc:`fix spring <fix_spring>` command
   to tether the molecule in place.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.  None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.

This fix computes a global scalar which can be accessed by various
:doc:`output commands <Howto_output>`.  The scalar is the distance the
group is moved by fix recenter.

This fix also computes global 3-vector which can be accessed by
various :doc:`output commands <Howto_output>`.  The 3 quantities in the
vector are xyz components of displacement applied to the group of
atoms by the fix.

The scalar and vector values calculated by this fix are "extensive".

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix should not be used with an x,y,z setting that causes a large
shift in the system on the first timestep, due to the requested COM
being very different from the initial COM.  This could cause atoms to
be lost, especially in parallel.  Instead, use the
:doc:`displace_atoms <displace_atoms>` command, which can be used to
move atoms a large distance.

Related commands
""""""""""""""""

:doc:`fix momentum <fix_momentum>`, :doc:`velocity <velocity>`

Default
"""""""

The option defaults are shift = fix group-ID, and units = lattice.
