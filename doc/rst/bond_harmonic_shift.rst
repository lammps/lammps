.. index:: bond\_style harmonic/shift

bond\_style harmonic/shift command
==================================

bond\_style harmonic/shift/omp command
======================================

Syntax
""""""


.. parsed-literal::

   bond_style harmonic/shift

Examples
""""""""


.. parsed-literal::

   bond_style harmonic/shift
   bond_coeff 5 10.0 0.5 1.0

Description
"""""""""""

The *harmonic/shift* bond style is a shifted harmonic bond that uses
the potential

.. image:: Eqs/bond_harmonic_shift.jpg
   :align: center

where r0 is the equilibrium bond distance, and rc the critical distance.
The potential is -Umin at r0 and zero at rc. The spring constant is
k = Umin / [ 2 (r0-rc)\^2].

The following coefficients must be defined for each bond type via the
:doc:`bond\_coeff <bond_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read\_data <read_data>`
or :doc:`read\_restart <read_restart>` commands:

* Umin (energy)

* r0 (distance)

* rc (distance)


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


This bond style can only be used if LAMMPS was built with the
USER-MISC package.  See the :doc:`Build package <Build_package>` doc
page for more info.

Related commands
""""""""""""""""

:doc:`bond\_coeff <bond_coeff>`, :doc:`delete\_bonds <delete_bonds>`,
:doc:`bond\_harmonic <bond_harmonic>`

**Default:** none


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
