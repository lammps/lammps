.. index:: angle\_style none

angle\_style none command
=========================

Syntax
""""""


.. parsed-literal::

   angle_style none

Examples
""""""""


.. parsed-literal::

   angle_style none

Description
"""""""""""

Using an angle style of none means angle forces and energies are not
computed, even if triplets of angle atoms were listed in the data file
read by the :doc:`read\_data <read_data>` command.

See the :doc:`angle\_style zero <angle_zero>` command for a way to
calculate angle statistics, but compute no angle interactions.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`angle\_style zero <angle_zero>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
