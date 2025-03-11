.. index:: fix efield/lepton

fix efield/lepton command
=========================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID efield/lepton V ...

* ID, group-ID are documented in the :doc:`fix <fix>` command
* style = *efield/lepton*
* V = electric potential (electric field * distance units)
* V must be a Lepton expression (see below)
* zero or more keyword/value pairs may be appended to args
* keyword = *region* or *step*

  .. parsed-literal::

       *region* value = region-ID
         region-ID = ID of region atoms must be in to have effect
       *step* value = h
         h = step size for numerical differentiation (distance units)

Examples
""""""""

.. code-block:: LAMMPS

   fix ex all efield/lepton "-E*x; E=1"
   fix dexx all efield/lepton "-0.5*x^2" step 1
   fix yukawa all efield/lepton "A*exp(-B*r)/r; r=abs(sqrt(x^2+y^2+z^2)); A=1; B=1" step 1e-6
   fix infp all efield/lepton "-abs(x)" step 1

   variable th equal 2*PI*ramp(0,1)
   fix erot all efield/lepton "-(x*cos(v_th)+y*sin(v_th))"

Description
"""""""""""

.. versionadded:: 4Feb2025

Add an electric potential :math:`V` that applies to a group of charged atoms a force :math:`\vec{F} = q \vec{E}`,
and to dipoles a force :math:`\vec{F} = (\vec{p} \cdot \nabla) \vec{E}` and torque :math:`\vec{T} = \vec{p} \times \vec{E}`,
where :math:`\vec{E} = - \nabla V`. The fix also evaluates the electrostatic energy (:math:`U_{q} = q V` and :math:`U_{p} = - \vec{p} \cdot \vec{E}`)
due to this potential when the :doc:`fix_modify energy yes <fix_modify>` command is specified (see below).

.. note::

   This command should be used instead of :doc:`fix efield <fix_efield>` if you want to impose a non-uniform electric field on a system with dipoles
   since the latter does not include the dipole force term. If you only have charges or if the electric field gradient is negligible,
   :doc:`fix efield <fix_efield>` should be used since it is faster.

The `Lepton library <https://simtk.org/projects/lepton>`_, that the *efield/lepton* fix style interfaces with, evaluates
the expression string at run time to compute the energy, forces, and torques. It creates an analytical representation
of :math:`V` and :math:`\vec{E}`, while the gradient force is computed using a central difference scheme

.. math::

   \vec{F} = \frac{|\vec{p}|}{2h} \left[ \vec{E}(\vec{x} + h \hat{p}) - \vec{E}(\vec{x} - h \hat{p}) \right] .

The Lepton expression must be either enclosed in quotes or must not contain any whitespace so that LAMMPS
recognizes it as a single keyword. More on valid Lepton expressions below. The final Lepton expression must
be a function of only :math:`x, y, z`, which refer to the current *unwrapped* coordinates of the atoms to ensure continuity.
Special care must be taken when using this fix with periodic boundary conditions or box-changing commands.

----------

.. include:: lepton_expression.rst

----------

If the *region* keyword is used, the atom must also be in the specified
geometric :doc:`region <region>` in order to be affected by the potential.

The *step* keyword is required when :doc:`atom_style dipole <atom_style>` is used and the electric field is non-uniform.

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.

The :doc:`fix_modify <fix_modify>` *energy* option is supported by this
fix to add the potential energy defined above to the global potential energy
of the system as part of :doc:`thermodynamic output <thermo_style>`.
The default setting for this fix is :doc:`fix_modify energy no <fix_modify>`.

The :doc:`fix_modify <fix_modify>` *virial* option is supported by this
fix to add the contribution due to the added ***forces*** on charges and dipoles
to both the global pressure and per-atom stress of the system via the
:doc:`compute pressure <compute_pressure>` and :doc:`compute stress/atom
<compute_stress_atom>` commands. The former can be accessed by
:doc:`thermodynamic output <thermo_style>`. The default setting for
this fix is :doc:`fix_modify virial no <fix_modify>`.

The :doc:`fix_modify <fix_modify>` *respa* option is supported by this
fix. This allows to set at which level of the :doc:`r-RESPA <run_style>`
integrator the fix adding its forces. Default is the outermost level.

This fix computes a global scalar and a global 3-vector of forces,
which can be accessed by various :doc:`output commands <Howto_output>`.
The scalar is the potential energy discussed above.
The vector is the total force added to the group of atoms.
The scalar and vector values calculated by this fix are "extensive".

This fix cannot be used with the *start/stop* keywords of
the :doc:`run <run>` command.

The forces due to this fix are imposed during an energy minimization,
invoked by the :doc:`minimize <minimize>` command. You should not
specify force components with a variable that has time-dependence for
use with a minimizer, since the minimizer increments the timestep as
the iteration count during the minimization.

.. note::

   If you want the electric potential energy to be included in the
   total potential energy of the system (the quantity being minimized),
   you MUST enable the :doc:`fix_modify <fix_modify>` *energy* option for this fix.

----------

Restrictions
""""""""""""

Fix style *efield/lepton* is part of the LEPTON package. It is only enabled if LAMMPS was built with that package.
See the :doc:`Build package <Build_package>` page for more info.


Related commands
""""""""""""""""

:doc:`fix efield <fix_efield>`

Default
"""""""

none
