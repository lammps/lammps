.. index:: bond_style none

bond_style none command
=======================

Syntax
""""""


.. code-block:: LAMMPS

   bond_style none

Examples
""""""""


.. code-blocK:: LAMMPS

   bond_style none

Description
"""""""""""

Using a bond style of none means bond forces and energies are not
computed, even if pairs of bonded atoms were listed in the data file
read by the :doc:`read_data <read_data>` command.

See the :doc:`bond_style zero <bond_zero>` command for a way to
calculate bond statistics, but compute no bond interactions.

Restrictions
""""""""""""
 none

**Related commands:** none

:doc:`bond_style zero <bond_zero>`

**Default:** none
