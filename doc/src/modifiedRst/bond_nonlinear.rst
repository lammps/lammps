.. index:: bond\_style nonlinear

bond\_style nonlinear command
=============================

bond\_style nonlinear/omp command
=================================

Syntax
""""""


.. parsed-literal::

   bond_style nonlinear

Examples
""""""""


.. parsed-literal::

   bond_style nonlinear
   bond_coeff 2 100.0 1.1 1.4

Description
"""""""""""

The *nonlinear* bond style uses the potential

.. math source doc: src/Eqs/bond_nonlinear.tex
.. math::

   E = \frac{\epsilon (r - r_0)^2}{ [ \lambda^2 - (r - r_0)^2 ]}


to define an anharmonic spring :ref:`(Rector) <Rector>` of equilibrium
length r0 and maximum extension lamda.

The following coefficients must be defined for each bond type via the
:doc:`bond\_coeff <bond_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read\_data <read_data>`
or :doc:`read\_restart <read_restart>` commands:

* epsilon (energy)
* r0 (distance)
* lamda (distance)


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


This bond style can only be used if LAMMPS was built with the MOLECULE
package.  See the :doc:`Build package <Build_package>` doc page for more
info.

Related commands
""""""""""""""""

:doc:`bond\_coeff <bond_coeff>`, :doc:`delete\_bonds <delete_bonds>`

**Default:** none


----------


.. _Rector:



**(Rector)** Rector, Van Swol, Henderson, Molecular Physics, 82, 1009 (1994).


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
