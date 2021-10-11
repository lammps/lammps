.. index:: angle_style charmm
.. index:: angle_style charmm/intel
.. index:: angle_style charmm/kk
.. index:: angle_style charmm/omp

angle_style charmm command
==========================

Accelerator Variants: *charmm/intel*, *charmm/kk*, *charmm/omp*

Syntax
""""""

.. code-block:: LAMMPS

   angle_style charmm

Examples
""""""""

.. code-block:: LAMMPS

   angle_style charmm
   angle_coeff 1 300.0 107.0 50.0 3.0

Description
"""""""""""

The *charmm* angle style uses the potential

.. math::

   E = K (\theta - \theta_0)^2 + K_{ub} (r - r_{ub})^2

with an additional Urey_Bradley term based on the distance :math:`r` between
the first and third atoms in the angle.  :math:`K`, :math:`\theta_0`,
:math:`K_{ub}`, and :math:`R_{ub}` are coefficients defined for each angle
type.

See :ref:`(MacKerell) <angle-MacKerell>` for a description of the CHARMM force
field.

The following coefficients must be defined for each angle type via the
:doc:`angle_coeff <angle_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* :math:`K` (energy)
* :math:`\theta_0` (degrees)
* :math:`K_{ub}` (energy/distance\^2)
* :math:`r_{ub}` (distance)

:math:`\theta_0` is specified in degrees, but LAMMPS converts it to
radians internally; hence :math:`K` is effectively energy per
radian\^2.

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

This angle style can only be used if LAMMPS was built with the
MOLECULE package.  See the :doc:`Build package <Build_package>` doc page
for more info.

Related commands
""""""""""""""""

:doc:`angle_coeff <angle_coeff>`

Default
"""""""

none

----------

.. _angle-MacKerell:

**(MacKerell)** MacKerell, Bashford, Bellott, Dunbrack, Evanseck, Field,
Fischer, Gao, Guo, Ha, et al, J Phys Chem, 102, 3586 (1998).
