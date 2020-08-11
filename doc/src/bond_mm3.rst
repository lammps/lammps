.. index:: bond_style mm3

bond_style mm3 command
======================

Syntax
""""""

.. code-block:: LAMMPS

   bond_style mm3

Examples
""""""""

.. code-block:: LAMMPS

   bond_style mm3
   bond_coeff 1 100.0 107.0

Description
"""""""""""

The *mm3* bond style uses the potential that is anharmonic in the bond
as defined in :ref:`(Allinger) <mm3-allinger1989>`

.. math::

   E = K (r - r_0)^2 \left[ 1 - 2.55(r-r_0) + (7/12) 2.55^2(r-r_0)^2 \right]

where :math:`r_0` is the equilibrium value of the bond, and :math:`K` is a
prefactor. The anharmonic prefactors have units angstrom\^(-n):
-2.55 angstrom\^(-1) and (7/12)2.55\^2 angstrom\^(-2). The code takes
care of the necessary unit conversion for these factors internally.
Note that the MM3 papers contains an error in Eq (1):
(7/12)2.55 should be replaced with (7/12)2.55\^2

The following coefficients must be defined for each bond type via the
:doc:`bond_coeff <bond_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* :math:`K` (energy/distance\^2)
* :math:`r_0` (distance)

Restrictions
""""""""""""

This bond style can only be used if LAMMPS was built with the
USER_YAFF package.  See the :doc:`Build package <Build_package>` doc
page for more info.

Related commands
""""""""""""""""

:doc:`bond_coeff <bond_coeff>`

**Default:** none

----------

.. _mm3-allinger1989:

**(Allinger)** Allinger, Yuh, Lii, JACS, 111(23), 8551-8566
(1989),
