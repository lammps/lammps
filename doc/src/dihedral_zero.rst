.. index:: dihedral\_style zero

dihedral\_style zero command
============================

Syntax
""""""


.. parsed-literal::

   dihedral_style zero *nocoeff*

Examples
""""""""


.. parsed-literal::

   dihedral_style zero
   dihedral_style zero nocoeff
   dihedral_coeff \*

Description
"""""""""""

Using a dihedral style of zero means dihedral forces and energies are
not computed, but the geometry of dihedral quadruplets is still
accessible to other commands.

As an example, the :doc:`compute dihedral/local <compute_dihedral_local>` command can be used to
compute the theta values for the list of quadruplets of dihedral atoms
listed in the data file read by the :doc:`read_data <read_data>`
command.  If no dihedral style is defined, this command cannot be
used.

The optional *nocoeff* flag allows to read data files with a DihedralCoeff
section for any dihedral style. Similarly, any dihedral\_coeff commands
will only be checked for the dihedral type number and the rest ignored.

Note that the :doc:`dihedral_coeff <dihedral_coeff>` command must be
used for all dihedral types, though no additional values are
specified.

Restrictions
""""""""""""
 none

**Related commands:** none

:doc:`dihedral_style none <dihedral_none>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
