.. index:: bond_coeff

bond_coeff command
==================

Syntax
""""""

.. code-block:: LAMMPS

   bond_coeff N args

* N = bond type (see asterisk form below)
* args = coefficients for one or more bond types

Examples
""""""""

.. code-block:: LAMMPS

   bond_coeff 5 80.0 1.2
   bond_coeff * 30.0 1.5 1.0 1.0
   bond_coeff 1*4 30.0 1.5 1.0 1.0
   bond_coeff 1 harmonic 200.0 1.0

Description
"""""""""""

Specify the bond force field coefficients for one or more bond types.
The number and meaning of the coefficients depends on the bond style.
Bond coefficients can also be set in the data file read by the
:doc:`read_data <read_data>` command or in a restart file.

N can be specified in one of two ways.  An explicit numeric value can
be used, as in the first example above.  Or a wild-card asterisk can be
used to set the coefficients for multiple bond types.  This takes the
form "\*" or "\*n" or "n\*" or "m\*n".  If N = the number of bond types,
then an asterisk with no numeric values means all types from 1 to N.  A
leading asterisk means all types from 1 to n (inclusive).  A trailing
asterisk means all types from n to N (inclusive).  A middle asterisk
means all types from m to n (inclusive).

Note that using a bond_coeff command can override a previous setting
for the same bond type.  For example, these commands set the coeffs
for all bond types, then overwrite the coeffs for just bond type 2:

.. code-block:: LAMMPS

   bond_coeff * 100.0 1.2
   bond_coeff 2 200.0 1.2

A line in a data file that specifies bond coefficients uses the exact
same format as the arguments of the bond_coeff command in an input
script, except that wild-card asterisks should not be used since
coefficients for all N types must be listed in the file.  For example,
under the "Bond Coeffs" section of a data file, the line that
corresponds to the first example above would be listed as

.. parsed-literal::

   5 80.0 1.2

----------

The list of all bond styles defined in LAMMPS is given on the
:doc:`bond_style <bond_style>` doc page.  They are also listed in more
compact form on the :doc:`Commands bond <Commands_bond>` doc page.

On either of those pages, click on the style to display the formula it
computes and its coefficients as specified by the associated
bond_coeff command.

----------

Restrictions
""""""""""""

This command must come after the simulation box is defined by a
:doc:`read_data <read_data>`, :doc:`read_restart <read_restart>`, or
:doc:`create_box <create_box>` command.

A bond style must be defined before any bond coefficients are set,
either in the input script or in a data file.

Related commands
""""""""""""""""

:doc:`bond_style <bond_style>`

Default
"""""""

none
