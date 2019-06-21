.. index:: bond\_style class2

bond\_style class2 command
==========================

bond\_style class2/omp command
==============================

bond\_style class2/kk command
=============================

Syntax
""""""


.. parsed-literal::

   bond_style class2

Examples
""""""""


.. parsed-literal::

   bond_style class2
   bond_coeff 1 1.0 100.0 80.0 80.0

Description
"""""""""""

The *class2* bond style uses the potential

.. math::

   E = K_2 (r - r_0)^2 + K_3 (r - r_0)^3 + K_4 (r - r_0)^4 


where r0 is the equilibrium bond distance.

See :ref:`(Sun) <bond-Sun>` for a description of the COMPASS class2 force field.

The following coefficients must be defined for each bond type via the
:doc:`bond\_coeff <bond_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read\_data <read_data>`
or :doc:`read\_restart <read_restart>` commands:

* R0 (distance)
* K2 (energy/distance\^2)
* K3 (energy/distance\^3)
* K4 (energy/distance\^4)


----------


Styles with a *gpu*\ , *intel*\ , *kk*\ , *omp*\ , or *opt* suffix are
functionally the same as the corresponding style without the suffix.
They have been optimized to run faster, depending on your available
hardware, as discussed on the :doc:`Speed packages <Speed_packages>` doc
page.  The accelerated styles take the same arguments and should
produce the same results, except for round-off and precision issues.

These accelerated styles are part of the GPU, USER-INTEL, KOKKOS,
USER-OMP and OPT packages, respectively.  They are only enabled if
LAMMPS was built with those packages.  See the :doc:`Build package <Build_package>` doc page for more info.

You can specify the accelerated styles explicitly in your input script
by including their suffix, or you can use the :doc:`-suffix command-line switch <Run_options>` when you invoke LAMMPS, or you can use the
:doc:`suffix <suffix>` command in your input script.

See the :doc:`Speed packages <Speed_packages>` doc page for more
instructions on how to use the accelerated styles effectively.


----------


Restrictions
""""""""""""


This bond style can only be used if LAMMPS was built with the CLASS2
package.  See the :doc:`Build package <Build_package>` doc page for more
info.

Related commands
""""""""""""""""

:doc:`bond\_coeff <bond_coeff>`, :doc:`delete\_bonds <delete_bonds>`

**Default:** none


----------


.. _bond-Sun:



**(Sun)** Sun, J Phys Chem B 102, 7338-7364 (1998).


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
