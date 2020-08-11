.. index:: fix efield

fix efield command
==================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID efield ex ey ez keyword value ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* efield = style name of this fix command
* ex,ey,ez = E-field component values (electric field units)
* any of ex,ey,ez can be a variable (see below)
* zero or more keyword/value pairs may be appended to args
* keyword = *region* or *energy*

  .. parsed-literal::

       *region* value = region-ID
         region-ID = ID of region atoms must be in to have added force
       *energy* value = v_name
         v_name = variable with name that calculates the potential energy of each atom in the added E-field

Examples
""""""""

.. code-block:: LAMMPS

   fix kick external-field efield 1.0 0.0 0.0
   fix kick external-field efield 0.0 0.0 v_oscillate

Description
"""""""""""

Add a force F = qE to each charged atom in the group due to an
external electric field being applied to the system.  If the system
contains point-dipoles, also add a torque on the dipoles due to the
external electric field.

For charges, any of the 3 quantities defining the E-field components
can be specified as an equal-style or atom-style
:doc:`variable <variable>`, namely *ex*\ , *ey*\ , *ez*\ .  If the value is a
variable, it should be specified as v_name, where name is the variable
name.  In this case, the variable will be evaluated each timestep, and
its value used to determine the E-field component.

For point-dipoles, equal-style variables can be used, but atom-style
variables are not currently supported, since they imply a spatial
gradient in the electric field which means additional terms with
gradients of the field are required for the force and torque on
dipoles.

Equal-style variables can specify formulas with various mathematical
functions, and include :doc:`thermo_style <thermo_style>` command
keywords for the simulation box parameters and timestep and elapsed
time.  Thus it is easy to specify a time-dependent E-field.

Atom-style variables can specify the same formulas as equal-style
variables but can also include per-atom values, such as atom
coordinates.  Thus it is easy to specify a spatially-dependent E-field
with optional time-dependence as well.

If the *region* keyword is used, the atom must also be in the
specified geometric :doc:`region <region>` in order to have force added
to it.

----------

Adding a force or torque to atoms implies a change in their potential
energy as they move or rotate due to the applied E-field.

For dynamics via the "run" command, this energy can be optionally
added to the system's potential energy for thermodynamic output (see
below).  For energy minimization via the "minimize" command, this
energy must be added to the system's potential energy to formulate a
self-consistent minimization problem (see below).

The *energy* keyword is not allowed if the added field is a constant
vector (ex,ey,ez), with all components defined as numeric constants
and not as variables.  This is because LAMMPS can compute the energy
for each charged particle directly as E = -x dot qE = -q (x\*ex + y\*ey
+ z\*ez), so that -Grad(E) = F.  Similarly for point-dipole particles
the energy can be computed as E = -mu dot E = -(mux\*ex + muy\*ey +
muz\*ez).

The *energy* keyword is optional if the added force is defined with
one or more variables, and if you are performing dynamics via the
:doc:`run <run>` command.  If the keyword is not used, LAMMPS will set
the energy to 0.0, which is typically fine for dynamics.

The *energy* keyword is required if the added force is defined with
one or more variables, and you are performing energy minimization via
the "minimize" command for charged particles.  It is not required for
point-dipoles, but a warning is issued since the minimizer in LAMMPS
does not rotate dipoles, so you should not expect to be able to
minimize the orientation of dipoles in an applied electric field.

The *energy* keyword specifies the name of an atom-style
:doc:`variable <variable>` which is used to compute the energy of each
atom as function of its position.  Like variables used for *ex*\ , *ey*\ ,
*ez*\ , the energy variable is specified as v_name, where name is the
variable name.

Note that when the *energy* keyword is used during an energy
minimization, you must insure that the formula defined for the
atom-style :doc:`variable <variable>` is consistent with the force
variable formulas, i.e. that -Grad(E) = F.  For example, if the force
due to the electric field were a spring-like F = kx, then the energy
formula should be E = -0.5kx\^2.  If you don't do this correctly, the
minimization will not converge properly.

----------

**Restart, fix_modify, output, run start/stop, minimize info:**

No information about this fix is written to :doc:`binary restart files <restart>`.

The :doc:`fix_modify <fix_modify>` *energy* option is supported by this
fix to add the potential "energy" inferred by the added force due to
the electric field to the system's potential energy as part of
:doc:`thermodynamic output <thermo_style>`.  This is a fictitious
quantity but is needed so that the :doc:`minimize <minimize>` command
can include the forces added by this fix in a consistent manner.
I.e. there is a decrease in potential energy when atoms move in the
direction of the added force due to the electric field.

The :doc:`fix_modify <fix_modify>` *virial* option is supported by this
fix to add the contribution due to the added forces on atoms to the
system's virial as part of :doc:`thermodynamic output <thermo_style>`.
The default is *virial no*

The :doc:`fix_modify <fix_modify>` *respa* option is supported by this
fix. This allows to set at which level of the :doc:`r-RESPA <run_style>`
integrator the fix adding its forces. Default is the outermost level.

This fix computes a global scalar and a global 3-vector of forces,
which can be accessed by various :doc:`output commands <Howto_output>`.
The scalar is the potential energy discussed above.  The vector is the
total force added to the group of atoms.  The scalar and vector values
calculated by this fix are "extensive".

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

This fix is part of the MISC package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package <Build_package>` doc page for more info.

Related commands
""""""""""""""""

:doc:`fix addforce <fix_addforce>`

**Default:** none
