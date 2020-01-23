.. index:: dihedral\_coeff

dihedral\_coeff command
=======================

Syntax
""""""


.. parsed-literal::

   dihedral_coeff N args

* N = dihedral type (see asterisk form below)
* args = coefficients for one or more dihedral types

Examples
""""""""


.. parsed-literal::

   dihedral_coeff 1 80.0 1 3
   dihedral_coeff \* 80.0 1 3 0.5
   dihedral_coeff 2\* 80.0 1 3 0.5

Description
"""""""""""

Specify the dihedral force field coefficients for one or more dihedral types.
The number and meaning of the coefficients depends on the dihedral style.
Dihedral coefficients can also be set in the data file read by the
:doc:`read_data <read_data>` command or in a restart file.

N can be specified in one of two ways.  An explicit numeric value can
be used, as in the 1st example above.  Or a wild-card asterisk can be
used to set the coefficients for multiple dihedral types.  This takes the
form "\*" or "\*n" or "n\*" or "m\*n".  If N = the number of dihedral types,
then an asterisk with no numeric values means all types from 1 to N.  A
leading asterisk means all types from 1 to n (inclusive).  A trailing
asterisk means all types from n to N (inclusive).  A middle asterisk
means all types from m to n (inclusive).

Note that using a dihedral\_coeff command can override a previous setting
for the same dihedral type.  For example, these commands set the coeffs
for all dihedral types, then overwrite the coeffs for just dihedral type 2:


.. parsed-literal::

   dihedral_coeff \* 80.0 1 3
   dihedral_coeff 2 200.0 1 3

A line in a data file that specifies dihedral coefficients uses the exact
same format as the arguments of the dihedral\_coeff command in an input
script, except that wild-card asterisks should not be used since
coefficients for all N types must be listed in the file.  For example,
under the "Dihedral Coeffs" section of a data file, the line that
corresponds to the 1st example above would be listed as


.. parsed-literal::

   1 80.0 1 3

The :doc:`dihedral_style class2 <dihedral_class2>` is an exception to
this rule, in that an additional argument is used in the input script
to allow specification of the cross-term coefficients.  See its doc
page for details.

.. note::

   When comparing the formulas and coefficients for various LAMMPS
   dihedral styles with dihedral equations defined by other force fields,
   note that some force field implementations divide/multiply the energy
   prefactor *K* by the multiple number of torsions that contain the J-K
   bond in an I-J-K-L torsion.  LAMMPS does not do this, i.e. the listed
   dihedral equation applies to each individual dihedral.  Thus you need
   to define *K* appropriately to account for this difference if
   necessary.


----------


The list of all dihedral styles defined in LAMMPS is given on the
:doc:`dihedral_style <dihedral_style>` doc page.  They are also listed
in more compact form on the :ref:`Commands dihedral <dihedral>` doc page.

On either of those pages, click on the style to display the formula it
computes and its coefficients as specified by the associated
dihedral\_coeff command.


----------


Restrictions
""""""""""""


This command must come after the simulation box is defined by a
:doc:`read_data <read_data>`, :doc:`read_restart <read_restart>`, or
:doc:`create_box <create_box>` command.

A dihedral style must be defined before any dihedral coefficients are
set, either in the input script or in a data file.

Related commands
""""""""""""""""

:doc:`dihedral_style <dihedral_style>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
