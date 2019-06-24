.. index:: dihedral\_style fourier

dihedral\_style fourier command
===============================

dihedral\_style fourier/intel command
=====================================

dihedral\_style fourier/omp command
===================================

Syntax
""""""


.. parsed-literal::

   dihedral_style fourier

Examples
""""""""


.. parsed-literal::

   dihedral_style fourier
   dihedral_coeff 1 3 -0.846200 3 0.0 7.578800 1 0 0.138000 2 -180.0

Description
"""""""""""

The *fourier* dihedral style uses the potential:

.. math::

   E = \sum_{i=1,m} K_i  [ 1.0 + \cos ( n_i \phi - d_i ) ]


The following coefficients must be defined for each dihedral type via the
:doc:`dihedral\_coeff <dihedral_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read\_data <read_data>`
or :doc:`read\_restart <read_restart>` commands:

* m (integer >=1)
* K1 (energy)
* n1 (integer >= 0)
* d1 (degrees)
* [...]
* Km (energy)
* nm (integer >= 0)
* dm (degrees)


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


This angle style can only be used if LAMMPS was built with the
USER\_MISC package.  See the :doc:`Build package <Build_package>` doc
page for more info.

Related commands
""""""""""""""""

:doc:`dihedral\_coeff <dihedral_coeff>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
