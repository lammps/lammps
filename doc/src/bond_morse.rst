.. index:: bond_style morse

bond_style morse command
========================

bond_style morse/omp command
============================

Syntax
""""""


.. code-block:: LAMMPS

   bond_style morse

Examples
""""""""


.. code-block:: LAMMPS

   bond_style morse
   bond_coeff 5 1.0 2.0 1.2

Description
"""""""""""

The *morse* bond style uses the potential

.. math::

   E = D \left[ 1 - e^{-\alpha (r - r_0)} \right]^2


where :math:`r_0` is the equilibrium bond distance, :math:`\alpha` is a stiffness
parameter, and :math:`D` determines the depth of the potential well.

The following coefficients must be defined for each bond type via the
:doc:`bond\_coeff <bond_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read\_data <read_data>`
or :doc:`read\_restart <read_restart>` commands:

* :math:`D` (energy)
* :math:`\alpha` (inverse distance)
* :math:`r_0` (distance)


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
