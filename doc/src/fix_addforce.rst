.. index:: fix addforce

fix addforce command
====================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID addforce fx fy fz keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* addforce = style name of this fix command
* fx,fy,fz = force component values (force units)

  .. parsed-literal::

       any of fx,fy,fz can be a variable (see below)

* zero or more keyword/value pairs may be appended to args
* keyword = *every* or *region* or *energy*

  .. parsed-literal::

       *every* value = Nevery
         Nevery = add force every this many timesteps
       *region* value = region-ID
         region-ID = ID of region atoms must be in to have added force
       *energy* value = v_name
         v_name = variable with name that calculates the potential energy of each atom in the added force field

Examples
""""""""

.. parsed-literal::

   fix kick flow addforce 1.0 0.0 0.0
   fix kick flow addforce 1.0 0.0 v_oscillate
   fix ff boundary addforce 0.0 0.0 v_push energy v_espace

Description
"""""""""""

Add fx,fy,fz to the corresponding component of force for each atom in
the group.  This command can be used to give an additional push to
atoms in a simulation, such as for a simulation of Poiseuille flow in
a channel.

Any of the 3 quantities defining the force components can be specified
as an equal-style or atom-style :doc:`variable <variable>`, namely *fx*\ ,
*fy*\ , *fz*\ .  If the value is a variable, it should be specified as
v\_name, where name is the variable name.  In this case, the variable
will be evaluated each timestep, and its value(s) used to determine
the force component.

Equal-style variables can specify formulas with various mathematical
functions, and include :doc:`thermo_style <thermo_style>` command
keywords for the simulation box parameters and timestep and elapsed
time.  Thus it is easy to specify a time-dependent force field.

Atom-style variables can specify the same formulas as equal-style
variables but can also include per-atom values, such as atom
coordinates.  Thus it is easy to specify a spatially-dependent force
field with optional time-dependence as well.

If the *every* keyword is used, the *Nevery* setting determines how
often the forces are applied.  The default value is 1, for every
timestep.

If the *region* keyword is used, the atom must also be in the
specified geometric :doc:`region <region>` in order to have force added
to it.

----------

Adding a force to atoms implies a change in their potential energy as
they move due to the applied force field.  For dynamics via the "run"
command, this energy can be optionally added to the system's potential
energy for thermodynamic output (see below).  For energy minimization
via the "minimize" command, this energy must be added to the system's
potential energy to formulate a self-consistent minimization problem
(see below).

The *energy* keyword is not allowed if the added force is a constant
vector F = (fx,fy,fz), with all components defined as numeric
constants and not as variables.  This is because LAMMPS can compute
the energy for each atom directly as E = -x dot F = -(x\*fx + y\*fy +
z\*fz), so that -Grad(E) = F.

The *energy* keyword is optional if the added force is defined with
one or more variables, and if you are performing dynamics via the
:doc:`run <run>` command.  If the keyword is not used, LAMMPS will set
the energy to 0.0, which is typically fine for dynamics.

The *energy* keyword is required if the added force is defined with
one or more variables, and you are performing energy minimization via
the "minimize" command.  The keyword specifies the name of an
atom-style :doc:`variable <variable>` which is used to compute the
energy of each atom as function of its position.  Like variables used
for *fx*\ , *fy*\ , *fz*\ , the energy variable is specified as v\_name,
where name is the variable name.

Note that when the *energy* keyword is used during an energy
minimization, you must insure that the formula defined for the
atom-style :doc:`variable <variable>` is consistent with the force
variable formulas, i.e. that -Grad(E) = F.  For example, if the force
were a spring-like F = kx, then the energy formula should be E =
-0.5kx\^2.  If you don't do this correctly, the minimization will not
converge properly.

----------

Styles with a *gpu*\ , *intel*\ , *kk*\ , *omp*\ , or *opt* suffix are
functionally the same as the corresponding style without the suffix.
They have been optimized to run faster, depending on your available
hardware, as discussed on the :doc:`Speed packages <Speed_packages>` doc
page.  The accelerated styles take the same arguments and should
produce the same results, except for round-off and precision issues.

These accelerated styles are part of the GPU, USER-INTEL, KOKKOS,
USER-OMP and OPT packages, respectively.  They are only enabled if
LAMMPS was built with those packages.  See the :doc:`Build package <Build_package>` doc page for more info.

You can specify the accelerated styles explicitly in your input script
by including their suffix, or you can use the :doc:`-suffix command-line switch <Run_options>` when you invoke LAMMPS, or you can use the
:doc:`suffix <suffix>` command in your input script.

See the :doc:`Speed packages <Speed_packages>` doc page for more
instructions on how to use the accelerated styles effectively.

----------

**Restart, fix\_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.

The :doc:`fix_modify <fix_modify>` *energy* option is supported by this
fix to add the potential "energy" inferred by the added force to the
system's potential energy as part of :doc:`thermodynamic output <thermo_style>`.  This is a fictitious quantity but is
needed so that the :doc:`minimize <minimize>` command can include the
forces added by this fix in a consistent manner.  I.e. there is a
decrease in potential energy when atoms move in the direction of the
added force.

The :doc:`fix_modify <fix_modify>` *virial* option is supported by this
fix to add the contribution due to the added forces on atoms to the
system's virial as part of :doc:`thermodynamic output <thermo_style>`.
The default is *virial no*

The :doc:`fix_modify <fix_modify>` *respa* option is supported by this
fix. This allows to set at which level of the :doc:`r-RESPA <run_style>`
integrator the fix is adding its forces. Default is the outermost
level.

This fix computes a global scalar and a global 3-vector of forces,
which can be accessed by various :doc:`output commands <Howto_output>`.
The scalar is the potential energy discussed above.  The vector is the
total force on the group of atoms before the forces on individual
atoms are changed by the fix.  The scalar and vector values calculated
by this fix are "extensive".

No parameter of this fix can be used with the *start/stop* keywords of
the :doc:`run <run>` command.

The forces due to this fix are imposed during an energy minimization,
invoked by the :doc:`minimize <minimize>` command.  You should not
specify force components with a variable that has time-dependence for
use with a minimizer, since the minimizer increments the timestep as
the iteration count during the minimization.

.. note::

   If you want the fictitious potential energy associated with the
   added forces to be included in the total potential energy of the
   system (the quantity being minimized), you MUST enable the
   :doc:`fix_modify <fix_modify>` *energy* option for this fix.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`fix setforce <fix_setforce>`, :doc:`fix aveforce <fix_aveforce>`

Default
"""""""

The option default for the every keyword is every = 1.
