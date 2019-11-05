.. index:: bond\_style none

bond\_style none command
========================

Syntax
""""""


.. parsed-literal::

   bond_style none

Examples
""""""""


.. parsed-literal::

   bond_style none

Description
"""""""""""

Using a bond style of none means bond forces and energies are not
computed, even if pairs of bonded atoms were listed in the data file
read by the :doc:`read\_data <read_data>` command.

See the :doc:`bond\_style zero <bond_zero>` command for a way to
calculate bond statistics, but compute no bond interactions.

Restrictions
""""""""""""
 none

**Related commands:** none

:doc:`bond\_style zero <bond_zero>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
