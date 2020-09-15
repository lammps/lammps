.. index:: fix cac/setvelocity

fix cac/setvelocity command
===========================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID cac/setvelocity vx vy vz keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* cac/setvelocity = style name of this fix command
* vx,vy,vz = force component values
* any of vx,vy,vz can be a variable (see below)
* zero or more keyword/value pairs may be appended to args
* keyword = *region*
  
  .. parsed-literal::
  
       *region* value = region-ID
         region-ID = ID of region atoms or element centroids must be in to have altered velocity


Examples
""""""""

.. code-block:: LAMMPS

   fix freeze indenter cac/setvelocity 0.0 0.0 0.0
   fix 2 edge cac/setvelocity NULL 0.0 0.0
   fix 2 edge cac/setvelocity NULL 0.0 v_oscillate

Description
"""""""""""

Set each component of velocity on each atom, or CAC element's nodal velocity, 
in the group to the specified values vx,vy,vz.  This erases all previously 
computed velocity on the atom/element, though additional fixes could alter velocity.  
This command can be used to displace certain atoms/elements in the simulation to
achieve some loading velocity etc.

Any of the vx,vy,vz values can be specified as NULL which means do not
alter the velocity component in that dimension.

Any of the 3 quantities defining the velocity components can be specified
as an equal-style or atom-style :doc:`variable <variable>`, namely *vx*\ ,
*vy*\ , *vz*\ .  If the value is a variable, it should be specified as
v\_name, where name is the variable name.  In this case, the variable
will be evaluated each timestep, and its value used to determine the
velocity component.

Equal-style variables can specify formulas with various mathematical
functions, and include :doc:`thermo_style <thermo_style>` command
keywords for the simulation box parameters and timestep and elapsed
time.  Thus it is easy to specify a time-dependent velocity field.

Atom-style variables can specify the same formulas as equal-style
variables but can also include per-atom values, such as atom
coordinates.  Thus it is easy to specify a spatially-dependent velocity
field with optional time-dependence as well.

.. note::

   For finite elements an atom style variable that is a function of position
   will be evaluated by the position of the element's centroid.

If the *region* keyword is used, the atom or element centroid must also 
be in the specified geometric :doc:`region <region>` in order to have its
velocity altered.

----------

**Restart, fix_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.

This fix computes a global 3-vector of velocities, which can be accessed
by various :doc:`output commands <Howto_output>`.  This is the total
velocity of the group of atom/elements before the velocities are
changed by the fix.  The vector values calculated by this fix are
"extensive".

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

Restrictions
""""""""""""

This fix requires a CAC atom style

Related commands
""""""""""""""""

:doc:`fix cac/setforce <fix_cac_setforce>`

**Default:** none
