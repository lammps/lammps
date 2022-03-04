.. index:: angle_style cosine/buck6d

angle_style cosine/buck6d command
=================================

Syntax
""""""

.. code-block:: LAMMPS

   angle_style cosine/buck6d

Examples
""""""""

.. code-block:: LAMMPS

   angle_style cosine/buck6d
   angle_coeff 1  cosine/buck6d  1.978350  4  180.000000

Description
"""""""""""

The *cosine/buck6d* angle style uses the potential

.. math::

   E = K \left[ 1 + \cos(n\theta - \theta_0)\right]

where :math:`K` is the energy constant, :math:`n` is the periodic multiplicity and
:math:`\theta_0` is the equilibrium angle.

The coefficients must be defined for each angle type via the
:doc:`angle_coeff <angle_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands in the following order:

* :math:`K` (energy)
* :math:`n`
* :math:`\theta_0` (degrees)

:math:`\theta_0` is specified in degrees, but LAMMPS converts it to radians
internally.

Additional to the cosine term the *cosine/buck6d* angle style computes
the short range (vdW) interaction belonging to the
:doc:`pair_style buck6d <pair_buck6d_coul_gauss>` between the end atoms of the
angle.  For this reason this angle style only works in combination
with the :doc:`pair_style buck6d <pair_buck6d_coul_gauss>` styles and needs
the :doc:`special_bonds <special_bonds>` 1-3 interactions to be weighted
0.0 to prevent double counting.

----------

Restrictions
""""""""""""

*cosine/buck6d* can only be used in combination with the
:doc:`pair_style buck6d <pair_buck6d_coul_gauss>` style and with a
:doc:`special_bonds <special_bonds>` 0.0 weighting of 1-3 interactions.

This angle style can only be used if LAMMPS was built with the
MOFFF package.  See the :doc:`Build package <Build_package>` doc
page for more info.

Related commands
""""""""""""""""

:doc:`angle_coeff <angle_coeff>`

Default
"""""""

none
