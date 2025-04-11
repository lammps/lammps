.. index:: fix spring/self
.. index:: fix spring/self/kk

fix spring/self command
=======================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID spring/self K dir

* ID, group-ID are documented in :doc:`fix <fix>` command
* spring/self = style name of this fix command
* K = spring constant (force/distance units), can be a variable (see below)
* dir = xyz, xy, xz, yz, x, y, or z (optional, default: xyz)

Examples
""""""""

.. code-block:: LAMMPS

   fix tether boundary-atoms spring/self 10.0
   fix var all spring/self v_kvar
   fix zrest  move spring/self 10.0 z

Description
"""""""""""

Apply a spring force independently to each atom in the group to tether
it to its initial position.  The initial position for each atom is its
location at the time the fix command was issued.  At each timestep,
the magnitude of the force on each atom is -Kr, where r is the
displacement of the atom from its current position to its initial
position.  The distance r correctly takes into account any crossings
of periodic boundary by the atom since it was in its initial
position.

With the (optional) dir flag, one can select in which direction the
spring force is applied. By default, the restraint is applied in all
directions, but it can be limited to the xy-, xz-, yz-plane and the
x-, y-, or z-direction, thus restraining the atoms to a line or a
plane, respectively.

The force constant *k* can be specified as an equal-style or atom-style
:doc:`variable <variable>`.  If the value is a variable, it should be specified
as v_name, where name is the variable name.  In this case, the variable
will be evaluated each time step, and its value(s) will be used as
force constant for the spring force.

Equal-style variables can specify formulas with various mathematical
functions and include :doc:`thermo_style <thermo_style>` command
keywords for the simulation box parameters, time step, and elapsed time.
Thus, it is easy to specify a time-dependent force field.

Atom-style variables can specify the same formulas as equal-style
variables but can also include per-atom values, such as atom
coordinates.  Thus, it is easy to specify a spatially-dependent force
field with optional time-dependence as well.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This fix writes the original coordinates of tethered atoms to
:doc:`binary restart files <restart>`, so that the spring effect will
be the same in a restarted simulation.  See the :doc:`read_restart
<read_restart>` command for info on how to re-specify a fix in an
input script that reads a restart file, so that the operation of the
fix continues in an uninterrupted fashion.

The :doc:`fix_modify <fix_modify>` *energy* option is supported by
this fix to add the energy stored in the per-atom springs to the
global potential energy of the system as part of :doc:`thermodynamic
output <thermo_style>`.  The default setting for this fix is
:doc:`fix_modify energy no <fix_modify>`.

The :doc:`fix_modify <fix_modify>` *respa* option is supported by
this fix. This allows to set at which level of the :doc:`r-RESPA <run_style>`
integrator the fix is adding its forces. Default is the outermost level.

This fix computes a global scalar which can be accessed by various
:doc:`output commands <Howto_output>`.  The scalar is an energy which is
the sum of the spring energy for each atom, where the per-atom energy
is 0.5 \* K \* r\^2.  The scalar value calculated by this fix is
"extensive".

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

The forces due to this fix are imposed during an energy minimization,
invoked by the :doc:`minimize <minimize>` command.

.. note::

   If you want the per-atom spring energy to be included in the
   total potential energy of the system (the quantity being minimized),
   you MUST enable the :doc:`fix_modify <fix_modify>` *energy* option for
   this fix.

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

The KOKKOS version, *fix spring/self/kk* may only be used with a constant
value of K, not a variable.

Related commands
""""""""""""""""

:doc:`fix drag <fix_drag>`, :doc:`fix spring <fix_spring>`,
:doc:`fix smd <fix_smd>`, :doc:`fix spring/rg <fix_spring_rg>`

Default
"""""""

none
