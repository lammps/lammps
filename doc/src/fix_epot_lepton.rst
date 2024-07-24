.. index:: fix epot/lepton

fix epot command
==================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID epot/lepton V ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* style = *epot/lepton*
* V = electric potential (energy/charge units)
* V must be a Lepton expression (see below)
* zero or one keyword/value pair may be appended to args
* keyword = *region*

  .. parsed-literal::

       *region* value = region-ID
         region-ID = ID of region atoms must be in to have effect

Examples
""""""""

.. code-block:: LAMMPS

   fix ex all epot/lepton "-E*x; E = 1"
   fix dexx all epot/lepton "-0.5*x^2"
   fix yukawa all epot/lepton "A*exp(-B*r)/r; r = abs(sqrt(x^2+y^2+z^2)); A = 1; B = 1"
   fix infp all epot/lepton "-abs(x)"

Description
"""""""""""

.. versionadded:: 27Jun2024

Add an electric potential :math:`V` that applies to a group of charged atoms a force
:math:`\vec{F} = q \vec{E}`, and to dipoles a force :math:`\vec{F} = (\vec{p} \cdot) \nabla \vec{E}`
and torque :math:`\vec{T} = \vec{p} \times \vec{E}`, where :math:`\vec{E} = - \nabla V`.

The `Lepton library <https://simtk.org/projects/lepton>`_, that the
*epot/lepton* fix style interfaces with, evaluates this expression string at
run time to compute the energy, forces, and torques. It creates an analytical
representation of :math:`V` and its first- and second-order derivatives.

The Lepton expression must be either enclosed in quotes or must not
contain any whitespace so that LAMMPS recognizes it as a single keyword.
More on valid Lepton expressions below. 

----------

.. include:: lepton_expression.rst

----------

If the *region* keyword is used, the atom must also be in the specified 
geometric :doc:`region <region>` in order to be affected by the potential.

----------

For dynamics via the "run" command and energy minimization via the "minimize" command, 
the energies :math:`U_{q} = q V` and :math:`U_{p} = - \vec{p} \cdot{E}`  can be optionally 
added to the system's potential energy for thermodynamic output (see below).  

----------

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.

The :doc:`fix_modify <fix_modify>` *energy* option is supported by this
fix to add the potential energy inferred by the added force due to the
electric field to the global potential energy of the system as part of
:doc:`thermodynamic output <thermo_style>`.  The default setting for
this fix is :doc:`fix_modify energy no <fix_modify>`.  

The :doc:`fix_modify <fix_modify>` *virial* option is supported by this
fix to add the contribution due to the added forces on atoms to both the
global pressure and per-atom stress of the system via the :doc:`compute
pressure <compute_pressure>` and :doc:`compute stress/atom
<compute_stress_atom>` commands.  The former can be accessed by
:doc:`thermodynamic output <thermo_style>`.  The default setting for
this fix is :doc:`fix_modify virial no <fix_modify>`.

The :doc:`fix_modify <fix_modify>` *respa* option is supported by this
fix. This allows to set at which level of the :doc:`r-RESPA <run_style>`
integrator the fix adding its forces. Default is the outermost level.

This fix computes a global scalar and a global 3-vector of forces,
which can be accessed by various :doc:`output commands
<Howto_output>`. The scalar is the potential energy discussed above.
The vector is the total force added to the group of atoms. The scalar
and vector values calculated by this fix are "extensive".

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

.. note::

   The potential :math:`V` can only depend on the spatial coordinates :math:`x, y, z`.
   The potential will not be updated with periodic boundary conditions and the 
   :doc:`fix deform <fix_deform>` command. 
   
----------

Restrictions
""""""""""""

Fix style *epot/lepton* is part of the LEPTON package. It is only
enabled if LAMMPS was built with that package.  See the :doc:`Build
package <Build_package>` page for more info.


Related commands
""""""""""""""""

:doc:`fix efield <fix_efield>`

Default
"""""""

none
