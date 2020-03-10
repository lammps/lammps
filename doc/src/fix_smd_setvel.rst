.. index:: fix smd/setvel

fix smd/setvel command
======================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID smd/setvel vx vy vz keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* smd/setvel = style name of this fix command
* vx,vy,vz = velocity component values
* any of vx,vy,vz can be a variable (see below)
* zero or more keyword/value pairs may be appended to args
* keyword = *region*

  .. parsed-literal::

       *region* value = region-ID
         region-ID = ID of region particles must be in to have their velocities set

Examples
""""""""

.. parsed-literal::

   fix top_velocity top_group setvel 1.0 0.0 0.0

Description
"""""""""""

Set each component of velocity on each particle in the group to the specified
values vx,vy,vz, regardless of the forces acting on the particle.  This command can
be used to impose velocity boundary conditions.

Any of the vx,vy,vz values can be specified as NULL which means do not
alter the velocity component in that dimension.

This fix is indented to be used together with a time integration fix.

Any of the 3 quantities defining the velocity components can be specified
as an equal-style or atom-style :doc:`variable <variable>`, namely *vx*\ ,
*vy*\ , *vz*\ .  If the value is a variable, it should be specified as
v\_name, where name is the variable name.  In this case, the variable
will be evaluated each timestep, and its value used to determine the
force component.

Equal-style variables can specify formulas with various mathematical
functions, and include :doc:`thermo_style <thermo_style>` command
keywords for the simulation box parameters and timestep and elapsed
time.  Thus it is easy to specify a time-dependent velocity field.

Atom-style variables can specify the same formulas as equal-style
variables but can also include per-atom values, such as atom
coordinates.  Thus it is easy to specify a spatially-dependent velocity
field with optional time-dependence as well.

If the *region* keyword is used, the particle must also be in the
specified geometric :doc:`region <region>` in order to have its velocity set by this command.

----------

**Restart, fix\_modify, output, run start/stop, minimize info:**

Currently, no part of USER-SMD supports restarting nor minimization
None of the :doc:`fix_modify <fix_modify>` options
are relevant to this fix.

This fix computes a global 3-vector of forces, which can be accessed
by various :doc:`output commands <Howto_output>`.  This is the total
force on the group of atoms.  The vector values calculated by this fix
are "extensive".

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

Restrictions
""""""""""""

This fix is part of the USER-SMD package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

**Related commands:** none

**Default:** none
