.. index:: fix cac/addforce

fix cac/addforce command
========================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID cac/addforce fx fy fz keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* cac/addforce = style name of this fix command
* fx,fy,fz = force component values
* any of fx,fy,fz can be a variable (see below)
* zero or more keyword/value pairs may be appended to args
* keyword = *region*
  
  .. parsed-literal::
  
       *region* value = region-ID
         region-ID = ID of region atoms or element centroids must be in to have altered force

Examples
""""""""

.. code-block:: LAMMPS

   fix freeze indenter cac/addforce 0.0 0.0 0.0
   fix 2 edge cac/addforce NULL 0.0 0.0
   fix 2 edge cac/addforce NULL 0.0 v_oscillate

Description
"""""""""""

Add to each component of force on each atom, or CAC element's nodal force, 
in the group to the specified values fx,fy,fz.  This adds force to all previously 
computed forces on the atom/element, though additional fixes could add new forces.  
This command can be used to add body forces to atoms/elements in the simulation by 
either adding to the atom's force vector or adding to every nodal force of a finite 
element. It can thus only add rigid body force to finite elements.

.. warning::

   the nodal forces are not the applied force vector in FEA; it is
   for all intents and purposes the nodal acceleration times the type
   mass for the respective internal DOF.

Any of the fx,fy,fz values can be specified as NULL which means do not
alter the force component in that dimension.

Any of the 3 quantities defining the additional force components can be specified
as an equal-style or atom-style :doc:`variable <variable>`, namely *fx*\ ,
*fy*\ , *fz*\ .  If the value is a variable, it should be specified as
v_name, where name is the variable name.  In this case, the variable
will be evaluated each timestep, and its value used to determine the
addition to the force component.

Equal-style variables can specify formulas with various mathematical
functions, and include :doc:`thermo_style <thermo_style>` command
keywords for the simulation box parameters and timestep and elapsed
time.  Thus it is easy to specify a time-dependent force field.

Atom-style variables can specify the same formulas as equal-style
variables but can also include per-atom values, such as atom
coordinates.  Thus it is easy to specify a spatially-dependent force
field with optional time-dependence as well.

.. note::

   For finite elements an atom style variable that is a function of position
   will be evaluated by the position of the element's centroid.

If the *region* keyword is used, the atom or element centroid must also 
be in the specified geometric :doc:`region <region>` in order to have force added
to it.

----------

**Restart, fix_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.

This fix computes a global 3-vector of forces, which can be accessed
by various :doc:`output commands <Howto_output>`.  This is the total
force on the group of atom/elements before the forces are
changed by the fix.  The vector values calculated by this fix are
"extensive".

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

The forces due to this fix are imposed during an energy minimization,
invoked by the :doc:`minimize <minimize>` command, but you cannot set
forces to any value besides zero when performing a minimization.  Use
the :doc:`fix cac/addforce <fix_cac_addforce>` command if you want to apply a
non-zero force to atoms during a minimization.

Restrictions
""""""""""""

This fix requires a CAC atom style

Related commands
""""""""""""""""

:doc:`fix cac/addforce <fix_cac_addforce>`

**Default:** none
