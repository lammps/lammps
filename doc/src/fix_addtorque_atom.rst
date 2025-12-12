.. index:: fix addtorque/atom

fix addtorque/atom command
==========================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID addtorque/atom tx ty tz keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* addtorque/atom = style name of this fix command
* tx,ty,tz = torque component values (torque units)
  .. parsed-literal::

       any of tx,ty,tz can be a variable (see below)

* zero or more keyword/value pairs may be appended to args
* keyword = *every* or *region*

  .. parsed-literal::

       *every* value = Nevery
         Nevery = add torque every this many time steps
       *region* value = region-ID
         region-ID = ID of region atoms must be in to have added torque

Examples
""""""""

.. code-block:: LAMMPS

   fix kick flow addtorque/atom 1.0 0.0 0.0
   fix kick flow addtorque/atom 1.0 0.0 v_oscillate
   fix ff boundary addtorque/atom 0.0 0.0 v_push

Description
"""""""""""

.. versionadded:: 10Dec2025

This fix is intended to add a peratom torque of each individual
finite-sized atom in the group to the specified values. Unlike
:doc:`fix addtorque/group <fix_addtorque_group>`, it does not apply a
collective torque to a set of point particles.

Add :math:`(t_x,t_y,t_z)` to the corresponding component of the torque for each
atom in the group. Any of the three quantities defining the torque components,
namely :math:`t_x`, :math:`t_y`, and :math:`t_z`, can be specified as an
equal-style or atom-style :doc:`variable <variable>`.  If the value is a variable,
it should be specified as v_name, where name is the variable name.  In this case,
the variable will be evaluated each time step, and its value(s) will be used to
determine the torque component(s).

Equal-style variables can specify formulas with various mathematical
functions and include :doc:`thermo_style <thermo_style>` command
keywords for the simulation box parameters, time step, and elapsed time.
Thus, it is easy to specify a time-dependent torque field.

Atom-style variables can specify the same formulas as equal-style
variables but can also include per-atom values, such as atom
coordinates.  Thus, it is easy to specify a spatially-dependent torque
field with optional time-dependence as well.

If the *every* keyword is used, the *Nevery* setting determines how
often the torques are applied.  The default value is 1, for every
time step.

If the *region* keyword is used, the atom must also be in the
specified geometric :doc:`region <region>` in order to have torque added
to it.

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.

The :doc:`fix_modify <fix_modify>` *respa* option is supported by this
fix. This allows to set at which level of the :doc:`r-RESPA
<run_style>` integrator the fix is adding its torques. Default is the
outermost level.

This fix computes a global three-vector of torques which can be accessed
by various :doc:`output commands <Howto_output>`. The vector is the total
torque on the group of atoms before the torques on individual atoms are
changed by the fix.  The vector values calculated by this fix are "extensive".

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

The torques due to this fix are imposed during an energy minimization,
invoked by the :doc:`minimize <minimize>` command.  You should not
specify torque components with a variable that has time-dependence for
use with a minimizer, since the minimizer increments the time step as
the iteration count during the minimization.

.. note::

   This fix is not (currently) designed to be used with rigid fixes.
   While it will apply additional torques to all of the atoms in a
   rigid body as described above, there is not always an easy mapping
   between these peratom torques and the torque experienced by the body.

Restrictions
""""""""""""

Fix *addtorque/atom* is part of the GRANULAR package.  It is only
enabled if LAMMPS was built with that package.  See the :doc:`Build
package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`fix settorque/atom <fix_settorque_atom>`, :doc:`fix addforce <fix_addforce>`

Default
"""""""

The option default for the every keyword is every = 1.
