.. index:: bond_style zero

bond_style zero command
=======================

Syntax
""""""

.. code-block:: LAMMPS

   bond_style zero [nocoeff]

Examples
""""""""

.. code-block:: LAMMPS

   bond_style zero
   bond_style zero nocoeff
   bond_coeff *
   bond_coeff * 2.14

Description
"""""""""""

Using an bond style of zero means bond forces and energies are not
computed, but the geometry of bond pairs is still accessible to other
commands.

As an example, the :doc:`compute bond/local <compute_bond_local>`
command can be used to compute distances for the list of pairs of bond
atoms listed in the data file read by the :doc:`read_data <read_data>`
command.  If no bond style is defined, this command cannot be used.

The optional *nocoeff* flag allows to read data files with a BondCoeff
section for any bond style. Similarly, any bond\_coeff commands
will only be checked for the bond type number and the rest ignored.

Note that the :doc:`bond_coeff <bond_coeff>` command must be used for
all bond types. If specified, there can be only one value, which is
going to be used to assign an equilibrium distance, e.g. for use with
:doc:`fix shake <fix_shake>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`bond_style none <bond_none>`

**Default:** none
