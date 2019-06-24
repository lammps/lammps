.. index:: angle\_style charmm

angle\_style charmm command
===========================

angle\_style charmm/intel command
=================================

angle\_style charmm/kk command
==============================

angle\_style charmm/omp command
===============================

Syntax
""""""


.. parsed-literal::

   angle_style charmm

Examples
""""""""


.. parsed-literal::

   angle_style charmm
   angle_coeff 1 300.0 107.0 50.0 3.0

Description
"""""""""""

The *charmm* angle style uses the potential

.. math::

   E = K (\theta - \theta_0)^2 + K_{UB} (r - r_{UB})^2 


with an additional Urey\_Bradley term based on the distance *r* between
the 1st and 3rd atoms in the angle.  K, theta0, Kub, and Rub are
coefficients defined for each angle type.

See :ref:`(MacKerell) <angle-MacKerell>` for a description of the CHARMM force
field.

The following coefficients must be defined for each angle type via the
:doc:`angle\_coeff <angle_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read\_data <read_data>`
or :doc:`read\_restart <read_restart>` commands:

* K (energy/radian\^2)
* theta0 (degrees)
* K\_ub (energy/distance\^2)
* r\_ub (distance)

Theta0 is specified in degrees, but LAMMPS converts it to radians
internally; hence the units of K are in energy/radian\^2.


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

See :doc:`Speed packages <Speed_packages>` doc page for more
instructions on how to use the accelerated styles effectively.


----------


Restrictions
""""""""""""


This angle style can only be used if LAMMPS was built with the
MOLECULE package.  See the :doc:`Build package <Build_package>` doc page
for more info.

Related commands
""""""""""""""""

:doc:`angle\_coeff <angle_coeff>`

**Default:** none


----------


.. _angle-MacKerell:



**(MacKerell)** MacKerell, Bashford, Bellott, Dunbrack, Evanseck, Field,
Fischer, Gao, Guo, Ha, et al, J Phys Chem, 102, 3586 (1998).


.. _lws: http://lammps.sandia.gov
.. _ld: Manual.html
.. _lc: Commands_all.html
