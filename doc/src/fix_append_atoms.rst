.. index:: fix append/atoms

fix append/atoms command
========================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID append/atoms face ... keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* append/atoms = style name of this fix command
* face = *zhi*
* zero or more keyword/value pairs may be appended
* keyword = *basis* or *size* or *freq* or *temp* or *random* or *units*

  .. parsed-literal::

       *basis* values = M itype
         M = which basis atom
         itype = atom type (1-N) to assign to this basis atom
       *size* args = Lz
         Lz = z size of lattice region appended in a single event(distance units)
       *freq* args = freq
         freq = the number of timesteps between append events
       *temp* args = target damp seed extent
         target = target temperature for the region between zhi-extent and zhi (temperature units)
         damp = damping parameter (time units)
         seed = random number seed for langevin kicks
         extent = extent of thermostatted region (distance units)
       *random* args = xmax ymax zmax seed
         *xmax*, *ymax*, *zmax* = maximum displacement in particular direction (distance units)
         *seed* = random number seed for random displacement
       *units* value = *lattice* or *box*
         *lattice* = the wall position is defined in lattice units
         *box* = the wall position is defined in simulation box units

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all append/atoms zhi size 5.0 freq 295 units lattice
   fix 4 all append/atoms zhi size 15.0 freq 5 units box
   fix A all append/atoms zhi size 1.0 freq 1000 units lattice

Description
"""""""""""

This fix creates atoms on a lattice, appended on the zhi edge of the
system box.  This can be useful when a shock or wave is propagating
from zlo.  This allows the system to grow with time to accommodate an
expanding wave.  A simulation box must already exist, which is
typically created via the :doc:`create_box <create_box>` command.
Before using this command, a lattice must also be defined using the
:doc:`lattice <lattice>` command.

This fix will automatically freeze atoms on the zhi edge of the
system, so that overlaps are avoided when new atoms are appended.

The *basis* keyword specifies an atom type that will be assigned to
specific basis atoms as they are created.  See the
:doc:`lattice <lattice>` command for specifics on how basis atoms are
defined for the unit cell of the lattice.  By default, all created
atoms are assigned type = 1 unless this keyword specifies differently.

The *size* keyword defines the size in :math:`z` of the chunk of material to
be added.

The *random* keyword will give the atoms random displacements around
their lattice points to simulate some initial temperature.

The *temp* keyword will cause a region to be thermostatted with a
Langevin thermostat on the zhi boundary.  The size of the region is
measured from zhi and is set with the *extent* argument.

The *units* keyword determines the meaning of the distance units used
to define a wall position, but only when a numeric constant is used.
A *box* value selects standard distance units as defined by the
:doc:`units <units>` command (e.g., :math:`\mathrm{\mathring A}`
for units = real or metal.
A *lattice* value means the distance units are in lattice spacings.
The :doc:`lattice <lattice>` command must have been previously used to
define the lattice spacings.

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to
:doc:`binary restart files <restart>`.  None of the
:doc:`fix_modify <fix_modify>` options
are relevant to this fix.  No global or per-atom quantities are stored
by this fix for access by various :doc:`output commands <Howto_output>`.
No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.  This fix is not invoked during
:doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This fix style is part of the SHOCK package.  It is only enabled if
LAMMPS was built with that package. See the
:doc:`Build package <Build_package>` page for more info.

The boundary on which atoms are added with append/atoms must be
shrink/minimum.  The opposite boundary may be any boundary type other
than periodic.

Related commands
""""""""""""""""

:doc:`fix wall/piston <fix_wall_piston>` command

Default
"""""""

The keyword defaults are size = 0.0, freq = 0, units = lattice.  All
added atoms are of type 1 unless the basis keyword is used.
