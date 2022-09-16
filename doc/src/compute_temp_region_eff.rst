.. index:: compute temp/region/eff

compute temp/region/eff command
===============================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID temp/region/eff region-ID

* ID, group-ID are documented in :doc:`compute <compute>` command
* temp/region/eff = style name of this compute command
* region-ID = ID of region to use for choosing atoms

Examples
""""""""

.. code-block:: LAMMPS

   compute mine flow temp/region/eff boundary

Description
"""""""""""

Define a computation that calculates the temperature of a group of
nuclei and electrons in the :doc:`electron force field <pair_eff>`
model, within a geometric region using the electron force field.
A compute of this style can be used by commands that compute a
temperature (e.g., :doc:`thermo_modify <thermo_modify>`).

The operation of this compute is exactly like that described by the
:doc:`compute temp/region <compute_temp_region>` command, except that
the formula for the temperature itself includes the radial electron
velocity contributions, as discussed by the
:doc:`compute temp/eff <compute_temp_eff>` command.

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
LAMMPS was built with that package.
See the :doc:`Build package <Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`compute temp/region <compute_temp_region>`,
:doc:`compute temp/eff <compute_temp_eff>`,
:doc:`compute pressure <compute_pressure>`

Default
"""""""

none
