.. index:: bond\_style mm3

bond\_style mm3 command
=======================

Syntax
""""""


.. parsed-literal::

   bond_style mm3

Examples
""""""""


.. parsed-literal::

   bond_style mm3
   bond_coeff 1 100.0 107.0

Description
"""""""""""

The *mm3* bond style uses the potential that is anharmonic in the bond
as defined in :ref:`(Allinger) <mm3-allinger1989>`

.. math::

  E = K (r - r_0)^2 \left[ 1 - 2.55(r-r_0) + (7/12) 2.55^2(r-r_0)^2 \right]


where r0 is the equilibrium value of the bond, and K is a
prefactor. The anharmonic prefactors have units angstrom\^(-n):
-2.55 angstrom\^(-1) and (7/12)2.55\^2 angstrom\^(-2). The code takes
care of the necessary unit conversion for these factors internally.
Note that the MM3 papers contains an error in Eq (1):
(7/12)2.55 should be replaced with (7/12)2.55\^2

The following coefficients must be defined for each bond type via the
:doc:`bond\_coeff <bond_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read\_data <read_data>`
or :doc:`read\_restart <read_restart>` commands:

* K (energy/distance\^2)
* r0 (distance)

Restrictions
""""""""""""


This bond style can only be used if LAMMPS was built with the
USER\_YAFF package.  See the :doc:`Build package <Build_package>` doc
page for more info.

Related commands
""""""""""""""""

:doc:`bond\_coeff <bond_coeff>`

**Default:** none


----------


.. _mm3-allinger1989:



**(Allinger)** Allinger, Yuh, Lii, JACS, 111(23), 8551-8566
(1989),


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
