.. index:: improper\_style fourier

improper\_style fourier command
===============================

improper\_style fourier/omp command
===================================

Syntax
""""""


.. parsed-literal::

   improper_style fourier

Examples
""""""""


.. parsed-literal::

   improper_style fourier
   improper_coeff 1 100.0 180.0

Description
"""""""""""

The *fourier* improper style uses the following potential:

.. math::

   E = K [C_0 + C_1 \cos ( \omega) + C_2 \cos( 2 \omega) ] 


where K is the force constant and omega is the angle between the IL
axis and the IJK plane:

.. image:: Eqs/umbrella.jpg
   :align: center

If all parameter (see bellow) is not zero, the all the three possible angles will taken in account.

The following coefficients must be defined for each improper type via
the :doc:`improper\_coeff <improper_coeff>` command as in the example
above, or in the data file or restart files read by the
:doc:`read\_data <read_data>` or :doc:`read\_restart <read_restart>`
commands:

* K (energy)
* C0 (real)
* C1 (real)
* C2 (real)
* all  (integer >= 0)


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

:doc:`improper\_coeff <improper_coeff>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
