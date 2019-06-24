.. index:: dihedral\_style nharmonic

dihedral\_style nharmonic command
=================================

dihedral\_style nharmonic/omp command
=====================================

Syntax
""""""


.. parsed-literal::

   dihedral_style nharmonic

Examples
""""""""


.. parsed-literal::

   dihedral_style nharmonic
   dihedral_coeff \* 3 10.0 20.0 30.0

Description
"""""""""""

The *nharmonic* dihedral style uses the potential:

.. math::

   E = \sum_{n=1,n} A_n  \cos^{n-1}(\phi)


The following coefficients must be defined for each dihedral type via the
:doc:`dihedral\_coeff <dihedral_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read\_data <read_data>`
or :doc:`read\_restart <read_restart>` commands:

* n (integer >=1)
* A1 (energy)
* A2 (energy)
* ...
* An (energy)


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
