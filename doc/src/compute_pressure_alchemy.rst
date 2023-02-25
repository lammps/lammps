.. index:: compute pressure/alchemy

compute pressure/alchemy command
================================

Syntax
""""""

.. code-block:: LAMMPS

   compute ID group-ID pressure/alchemy fix-ID

* ID, group-ID are documented in :doc:`compute <compute>` command
* pressure/alchemy = style name of this compute command
* fix-ID = ID of :doc:`fix alchemy <fix_alchemy>` command

Examples
""""""""

.. code-block:: LAMMPS

   fix trans all alchemy
   compute mixed all pressure/alchemy trans
   thermo_modify press mixed

Description
"""""""""""

.. versionadded:: TBD

Define a compute style that makes the "mixed" system pressure available
for a system that uses the :doc:`fix alchemy <fix_alchemy>` command to
transform one topology to another.  This can be used in combination with
either :doc:`thermo_modify press <thermo_modify>` or :doc:`fix_modify
press <fix_modify>` to output and access a pressure consistent with the
simulated combined two topology system.

The actual pressure is determined with :doc:`compute pressure
<compute_pressure>` commands that are internally used by :doc:`fix
alchemy <fix_alchemy>` for each topology individually and then combined.
This command just extracts the information from the fix.

The ``examples/PACKAGES/alchemy`` folder contains an example input for this command.

----------

Output info
"""""""""""

This compute calculates a global scalar (the pressure) and a global
vector of length 6 (the pressure tensor), which can be accessed by
indices 1--6.  These values can be used by any command that uses global
scalar or vector values from a compute as input.  See the :doc:`Howto
output <Howto_output>` page for an overview of LAMMPS output options.

The ordering of values in the symmetric pressure tensor is as follows:
:math:`p_{xx},` :math:`p_{yy},` :math:`p_{zz},` :math:`p_{xy},`
:math:`p_{xz},` :math:`p_{yz}.`

The scalar and vector values calculated by this compute are "intensive".
The scalar and vector values will be in pressure :doc:`units <units>`.

Restrictions
""""""""""""

This compute is part of the REPLICA package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.


Related commands
""""""""""""""""

:doc:`fix alchemy <fix_alchemy>`, :doc:`compute pressure <compute_pressure>`,
:doc:`thermo_modify <thermo_modify>`, :doc:`fix_modify <fix_modify>`

Default
"""""""

none
