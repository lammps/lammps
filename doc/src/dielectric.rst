.. index:: dielectric

dielectric command
==================

Syntax
""""""

.. code-block:: LAMMPS

   dielectric value

* value = dielectric constant

Examples
""""""""

.. code-block:: LAMMPS

   dielectric 2.0

Description
"""""""""""

Set the dielectric constant for Coulombic interactions (pairwise and
long-range) to this value.  The constant is unitless, since it is used
to reduce the strength of the interactions.  The value is used in the
denominator of the formulas for Coulombic interactions (e.g., a value
of 4.0 reduces the Coulombic interactions to 25% of their default
strength).  See the :doc:`pair_style <pair_style>` command for more
details.

Restrictions
""""""""""""
 none

Related commands
""""""""""""""""

:doc:`pair_style <pair_style>`

Default
"""""""

.. code-block:: LAMMPS

   dielectric 1.0
