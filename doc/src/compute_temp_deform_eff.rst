.. index:: compute temp/deform/eff

compute temp/deform/eff command
===============================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID temp/deform/eff

* ID, group-ID are documented in :doc:`compute <compute>` command
* temp/deform/eff = style name of this compute command

Examples
""""""""

.. code-block:: LAMMPS

   compute myTemp all temp/deform/eff

Description
"""""""""""

Define a computation that calculates the temperature of a group of
nuclei and electrons in the :doc:`electron force field <pair_eff>`
model, after subtracting out a streaming velocity induced by the
simulation box changing size and/or shape, for example in a
non-equilibrium MD (NEMD) simulation.  The size/shape change is
induced by use of the :doc:`fix deform <fix_deform>` command.  A
compute of this style is created by the
:doc:`fix nvt/sllod/eff <fix_nvt_sllod_eff>` command to compute the thermal
temperature of atoms for thermostatting purposes.  A compute of this
style can also be used by any command that computes a temperature
(e.g., :doc:`thermo_modify <thermo_modify>`, :doc:`fix npt/eff <fix_nh_eff>`).

The calculation performed by this compute is exactly like that
described by the :doc:`compute temp/deform <compute_temp_deform>`
command, except that the formula for the temperature includes the
radial electron velocity contributions, as discussed by the :doc:`compute temp/eff <compute_temp_eff>` command.  Note that only the
translational degrees of freedom for each nuclei or electron are
affected by the streaming velocity adjustment.  The radial velocity
component of the electrons is not affected.

Output info
"""""""""""

This compute calculates a global scalar (the temperature) and a global
vector of length 6 (KE tensor), which can be accessed by indices 1--6.
These values can be used by any command that uses global scalar or
vector values from a compute as input.  See the
:doc:`Howto output <Howto_output>` page for an overview of LAMMPS output
options.

The scalar value calculated by this compute is "intensive".  The
vector values are "extensive".

The scalar value will be in temperature :doc:`units <units>`.  The
vector values will be in energy :doc:`units <units>`.

Restrictions
""""""""""""

This compute is part of the EFF package.  It is only enabled if
LAMMPS was built with that package.  See the
:doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`compute temp/ramp <compute_temp_ramp>`, :doc:`fix deform <fix_deform>`,
:doc:`fix nvt/sllod/eff <fix_nvt_sllod_eff>`

Default
"""""""

none
