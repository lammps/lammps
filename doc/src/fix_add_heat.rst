.. index:: fix add/heat

fix add/heat command
====================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID add/heat style args keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* add/heat = style name of this fix command
* style = *constant* or *linear* or *quartic*

  .. parsed-literal::

       *constant* args = *rate*
         *rate* = rate of heat flow (energy/time units)
       *linear* args = *Ttarget* *k*
         *Ttarget* = target temperature (temperature units)
         *k* = prefactor (energy/(time*temperature) units)
       *quartic* args = *Ttarget* *k*
         *Ttarget* = target temperature (temperature units)
         *k* = prefactor (energy/(time*temperature^4) units)

* zero or more keyword/value pairs may be appended to args
* keyword = *overwrite*

  .. parsed-literal::

       *overwrite* value = *yes* or *no*
         *yes* = sets current heat flow of particle
         *no* = adds to current heat flow of particle

Examples
""""""""

.. code-block:: LAMMPS

   fix 1 all add/heat constant v_heat
   fix 1 all add/heat linear 10.0 1.0 overwrite yes

Description
"""""""""""

This fix adds heat to particles with the temperature attribute every timestep
at a given rate. Note that this is an internal temperature of a particle intended
for use with non-atomistic models like the discrete element method.

For the *constant* style, heat is added at the specified rate. For the *linear* style,
heat is added at a rate of :math:`k (T_{target} - T)` where :math:`k` is the
specified prefactor, :math:`T_{target}` is the specified target temperature, and
:math:`T` is the temperature of the atom. This may be more representative of a
conductive process. For the *quartic* style, heat is added at a rate of
:math:`k (T_{target}^4 - T^4)`, akin to radiative heat transfer.

The rate or temperature can be can be specified as an equal-style or atom-style
:doc:`variable <variable>`.  If the value is a variable, it should be
specified as v_name, where name is the variable name.  In this case, the
variable will be evaluated each time step, and its value will be used to
determine the rate of heat added.

Equal-style variables can specify formulas with various mathematical
functions and include :doc:`thermo_style <thermo_style>` command
keywords for the simulation box parameters, time step, and elapsed time
to specify time-dependent heating.

Atom-style variables can specify the same formulas as equal-style
variables but can also include per-atom values, such as atom
coordinates to specify spatially-dependent heating.

If the *overwrite* keyword is set to *yes*, this fix will set the total
heat flow on a particle every timestep, overwriting contributions from pair
styles or other fixes. If *overwrite* is *no*, this fix will add heat on
top of other contributions.

----------

Restart, fix_modify, output, run start/stop, minimize info
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files <restart>`.
None of the :doc:`fix_modify <fix_modify>` options are relevant to this fix.
No global or per-atom quantities are stored by this fix for access by various
:doc:`output commands <Howto_output>`. No parameter of this fix can be used
with the *start/stop* keywords of the :doc:`run <run>` command.  This fix is
not invoked during :doc:`energy minimization <minimize>`.

Restrictions
""""""""""""

This pair style is part of the GRANULAR package.  It is
only enabled if LAMMPS was built with that package.
See the :doc:`Build package <Build_package>` page for more info.

This fix requires that atoms store temperature and heat flow
as defined by the :doc:`fix property/atom <fix_property_atom>` command or
included in certain atom styles, such as atom_style rheo/thermal.

Related commands
""""""""""""""""

:doc:`fix heat/flow <fix_heat_flow>`,
:doc:`fix property/atom <fix_property_atom>`,
:doc:`fix rheo/thermal <fix_rheo_thermal>`

Default
"""""""

The default for the *overwrite* keyword is *no*
