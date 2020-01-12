.. index:: dielectric

dielectric command
==================

Syntax
""""""


.. parsed-literal::

   dielectric value

* value = dielectric constant

Examples
""""""""


.. parsed-literal::

   dielectric 2.0

Description
"""""""""""

Set the dielectric constant for Coulombic interactions (pairwise and
long-range) to this value.  The constant is unitless, since it is used
to reduce the strength of the interactions.  The value is used in the
denominator of the formulas for Coulombic interactions - e.g. a value
of 4.0 reduces the Coulombic interactions to 25% of their default
strength.  See the :doc:`pair\_style <pair_style>` command for more
details.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`pair\_style <pair_style>`

Default
"""""""


.. parsed-literal::

   dielectric 1.0


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
