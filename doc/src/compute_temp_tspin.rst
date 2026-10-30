.. index:: compute temp/tspin

compute temp/tspin command
==========================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID temp/tspin

* ID, group-ID are documented in :doc:`compute <compute>` command
* temp/tspin = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute ts all temp/tspin
   thermo_style custom step temp c_ts c_ts[1]

Description
"""""""""""

.. versionadded:: TBD

Define a computation that calculates the instantaneous kinetic
temperature and the kinetic energy of the auxiliary spin degrees of
freedom of inertial spin dynamics (TSPIN) for all atoms in the group
that have a positive spin mass.

The kinetic energy is

.. math::

   K_s = \sum_i \frac{1}{2} m_s \left| \vec{v}^{s}_i \right|^2

where :math:`m_s` is the spin mass of atom i and :math:`\vec{v}^{s}_i`
its auxiliary spin velocity.  The temperature is

.. math::

   T_s = \frac{2K_s}{3N_s k_B},

where :math:`N_s` is the number of these atoms.

Each spin contributes three Cartesian degrees of freedom.  No degrees
of freedom are subtracted for spin momentum removal by
:doc:`velocity/tspin <velocity_tspin>`.

This temperature characterizes the auxiliary TSPIN degrees of freedom.
It differs from the torque-based estimator returned by
:doc:`compute spin <compute_spin>` and from the usual atomic
temperature *temp*.  The kinetic energy belongs to the sampling dynamics
and is distinct from the magnetic potential energy.

The kinetic energy is not included in the thermodynamic keyword *ke*,
which counts only the translational kinetic energy of the atoms, or in
*etotal* or *econserve*.  Reference it as *c_ID[1]* in
:doc:`thermo_style custom <thermo_style>` to output it.  For a
microcanonical run with :doc:`fix nve/tspin <fix_nve_tspin>` and a moving
lattice, with no thermostat or barostat active and with all relevant
potential-energy contributions included in *pe*, the conserved energy
includes the potential energy, the lattice kinetic energy, and this
auxiliary spin kinetic energy:

.. code-block:: LAMMPS

   compute      ts all temp/tspin
   variable     H equal etotal+c_ts[1]
   thermo_style custom step pe ke c_ts[1] v_H

Include all potential energy contributions in *pe*.  For example, the
energy of :doc:`fix spring/tspin <fix_spring_tspin>` requires
:doc:`fix_modify energy yes <fix_modify>`.  With *lattice frozen* and
zero atomic velocities, the lattice kinetic energy vanishes and the
expression reduces to ``pe+c_ts[1]``.

Output info
"""""""""""

This compute calculates a global scalar (the temperature) and a global
vector of length 1, whose only element, accessed with index 1, is the
kinetic energy.  These values can be used by any command that uses
global scalar or vector values from a compute as input.  See the
:doc:`Howto output <Howto_output>` page for an overview of LAMMPS output
options.

The scalar value is "intensive" and in temperature :doc:`units <units>`
(K in metal units).  The vector value is "extensive" and in energy
:doc:`units <units>`.

Restrictions
""""""""""""

This compute requires :doc:`atom_style spin <atom_style>` and the SPIN
package.  A TSPIN integrator must be defined first.

This compute is for spin-temperature output only.  It cannot be used
with :doc:`thermo_modify temp <thermo_modify>` or as a temperature
compute for atomic thermostats or :doc:`pressure calculations
<compute_pressure>`.

Related commands
""""""""""""""""

:doc:`compute ke <compute_ke>`,
:doc:`compute spin <compute_spin>`,
:doc:`velocity/tspin <velocity_tspin>`,
:doc:`fix nve/tspin <fix_nve_tspin>`,
:doc:`fix nvt/tspin <fix_nvt_tspin>`,
:doc:`fix langevin/tspin <fix_langevin_tspin>`

Default
"""""""

none
