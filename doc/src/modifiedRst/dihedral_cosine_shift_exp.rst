.. index:: dihedral\_style cosine/shift/exp

dihedral\_style cosine/shift/exp command
========================================

dihedral\_style cosine/shift/exp/omp command
============================================

Syntax
""""""


.. parsed-literal::

   dihedral_style cosine/shift/exp

Examples
""""""""


.. parsed-literal::

   dihedral_style cosine/shift/exp
   dihedral_coeff 1 10.0 45.0 2.0

Description
"""""""""""

The *cosine/shift/exp* dihedral style uses the potential

.. math::

E=-U_{min} 
\frac{e^{-a U(\theta,\theta_0)}-1}{e^a-1}
\quad\mbox{with}\quad
U(\theta,\theta_0)
=-0.5 \left(1+\cos(\theta-\theta_0) \right)


where Umin, theta, and a are defined for each dihedral type.

The potential is bounded between [-Umin:0] and the minimum is located
at the angle theta0. The a parameter can be both positive or negative
and is used to control the spring constant at the equilibrium.

The spring constant is given by k=a exp(a) Umin/ [2 (Exp(a)-1)].
For a>3 k/Umin = a/2 to better than 5% relative error. For negative
values of the a parameter, the spring constant is essentially zero,
and anharmonic terms takes over. The potential is furthermore well
behaved in the limit a->0, where it has been implemented to linear
order in a for a < 0.001.

The following coefficients must be defined for each dihedral type via
the :doc:`dihedral\_coeff <dihedral_coeff>` command as in the example
above, or in the data file or restart files read by the
:doc:`read\_data <read_data>` or :doc:`read\_restart <read_restart>`
commands:

* umin (energy)
* theta (angle)
* A (real number)


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


This dihedral style can only be used if LAMMPS was built with the
USER-MISC package.  See the :doc:`Build package <Build_package>` doc
page for more info.

Related commands
""""""""""""""""

:doc:`dihedral\_coeff <dihedral_coeff>`,
:doc:`angle\_cosine\_shift\_exp <angle_cosine_shift_exp>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
