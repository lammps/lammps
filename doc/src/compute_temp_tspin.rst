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

   compute 1 all temp/tspin
   thermo_style custom step temp c_1

Description
"""""""""""

.. versionadded:: TBD

Compute the instantaneous kinetic temperature of the auxiliary spin
degrees of freedom used by TSPIN,

.. math::

   T_s = \frac{2K_s}{3N_s k_B},

where :math:`K_s` is the auxiliary spin kinetic energy (see
:doc:`compute ke/tspin <compute_ke_tspin>`) and :math:`N_s` is the
number of magnetic atoms in the group.

Each spin contributes three Cartesian degrees of freedom.  No degrees
of freedom are subtracted for spin momentum removal by
:doc:`velocity/tspin <velocity_tspin>`.

This temperature characterizes the auxiliary TSPIN degrees of freedom.
It differs from the torque-based estimator returned by
:doc:`compute spin <compute_spin>` and from the usual atomic
temperature *temp*.

Output info
"""""""""""

This compute returns a global scalar that is "intensive" and in
temperature :doc:`units <units>` (K in metal units).

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

:doc:`compute ke/tspin <compute_ke_tspin>`,
:doc:`velocity/tspin <velocity_tspin>`,
:doc:`fix nvt/tspin <fix_nvt_tspin>`,
:doc:`fix langevin/tspin <fix_langevin_tspin>`

Default
"""""""

none
