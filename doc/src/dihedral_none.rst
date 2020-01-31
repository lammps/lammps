.. index:: dihedral\_style none

dihedral\_style none command
============================

Syntax
""""""


.. parsed-literal::

   dihedral_style none

Examples
""""""""


.. parsed-literal::

   dihedral_style none

Description
"""""""""""

Using a dihedral style of none means dihedral forces and energies are
not computed, even if quadruplets of dihedral atoms were listed in the
data file read by the :doc:`read_data <read_data>` command.

See the :doc:`dihedral_style zero <dihedral_zero>` command for a way to
calculate dihedral statistics, but compute no dihedral interactions.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`dihedral_style zero <dihedral_zero>`

**Default:** none


