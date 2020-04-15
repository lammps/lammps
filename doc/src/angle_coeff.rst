.. index:: angle_coeff

angle_coeff command
===================

Syntax
""""""

.. code-block:: LAMMPS

   angle_coeff N args

* N = angle type (see asterisk form below)
* args = coefficients for one or more angle types

Examples
""""""""

.. code-block:: LAMMPS

   angle_coeff 1 300.0 107.0
   angle_coeff * 5.0
   angle_coeff 2*10 5.0

Description
"""""""""""

Specify the angle force field coefficients for one or more angle types.
The number and meaning of the coefficients depends on the angle style.
Angle coefficients can also be set in the data file read by the
:doc:`read_data <read_data>` command or in a restart file.

N can be specified in one of two ways.  An explicit numeric value can
be used, as in the 1st example above.  Or a wild-card asterisk can be
used to set the coefficients for multiple angle types.  This takes the
form "\*" or "\*n" or "n\*" or "m\*n".  If N = the number of angle types,
then an asterisk with no numeric values means all types from 1 to N.  A
leading asterisk means all types from 1 to n (inclusive).  A trailing
asterisk means all types from n to N (inclusive).  A middle asterisk
means all types from m to n (inclusive).

Note that using an :doc:`angle_coeff <angle_coeff>` command can override a previous setting
for the same angle type.  For example, these commands set the coeffs
for all angle types, then overwrite the coeffs for just angle type 2:

.. code-block:: LAMMPS

   angle_coeff * 200.0 107.0 1.2
   angle_coeff 2 50.0 107.0

A line in a data file that specifies angle coefficients uses the exact
same format as the arguments of the :doc:`angle_coeff <angle_coeff>` command in an input
script, except that wild-card asterisks should not be used since
coefficients for all N types must be listed in the file.  For example,
under the "Angle Coeffs" section of a data file, the line that
corresponds to the 1st example above would be listed as

.. parsed-literal::

   1 300.0 107.0

The :doc:`angle_style class2 <angle_class2>` is an exception to this
rule, in that an additional argument is used in the input script to
allow specification of the cross-term coefficients.   See its
doc page for details.

----------

The list of all angle styles defined in LAMMPS is given on the
:doc:`angle_style <angle_style>` doc page.  They are also listed in more
compact form on the :ref:`Commands angle <angle>` doc
page.

On either of those pages, click on the style to display the formula it
computes and its coefficients as specified by the associated
:doc:`angle_coeff <angle_coeff>` command.

----------

Restrictions
""""""""""""

This command must come after the simulation box is defined by a
:doc:`read_data <read_data>`, :doc:`read_restart <read_restart>`, or
:doc:`create_box <create_box>` command.

An angle style must be defined before any angle coefficients are
set, either in the input script or in a data file.

Related commands
""""""""""""""""

:doc:`angle_style <angle_style>`

**Default:** none
