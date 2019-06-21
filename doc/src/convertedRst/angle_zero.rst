.. index:: angle\_style zero

angle\_style zero command
=========================

Syntax
""""""


.. parsed-literal::

   angle_style zero *nocoeff*

Examples
""""""""


.. parsed-literal::

   angle_style zero
   angle_style zero nocoeff
   angle_coeff \*
   angle_coeff \* 120.0

Description
"""""""""""

Using an angle style of zero means angle forces and energies are not
computed, but the geometry of angle triplets is still accessible to
other commands.

As an example, the :doc:`compute angle/local <compute_angle_local>`
command can be used to compute the theta values for the list of
triplets of angle atoms listed in the data file read by the
:doc:`read\_data <read_data>` command.  If no angle style is defined,
this command cannot be used.

The optional *nocoeff* flag allows to read data files with AngleCoeff
section for any angle style. Similarly, any angle\_coeff commands
will only be checked for the angle type number and the rest ignored.

Note that the :doc:`angle\_coeff <angle_coeff>` command must be used for
all angle types. If specified, there can be only one value, which is
going to be used to assign an equilibrium angle, e.g. for use with
:doc:`fix shake <fix_shake>`.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`angle\_style none <angle_none>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
