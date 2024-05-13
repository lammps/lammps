.. index:: angle_style spica
.. index:: angle_style spica/omp

angle_style spica command
=========================

Accelerator Variants: *spica/omp*

Syntax
""""""

.. code-block:: LAMMPS

   angle_style spica

   angle_style spica/omp

Examples
""""""""

.. code-block:: LAMMPS

   angle_style spica
   angle_coeff 1 300.0 107.0

Description
"""""""""""

The *spica* angle style is a combination of the harmonic angle potential,

.. math::

   E = K (\theta - \theta_0)^2

where :math:`\theta_0` is the equilibrium value of the angle and
:math:`K` a prefactor, with the *repulsive* part of the non-bonded
*lj/spica* pair style between the atoms 1 and 3.  This angle potential is
intended for coarse grained MD simulations with the SPICA (formerly called SDK) parameterization
using the :doc:`pair_style lj/spica <pair_spica>`.  Relative to the
pair_style *lj/spica*, however, the energy is shifted by
:math:`\epsilon`, to avoid sudden jumps.  Note that the usual 1/2 factor
is included in :math:`K`.

The following coefficients must be defined for each angle type via the
:doc:`angle_coeff <angle_coeff>` command as in the example above:

* :math:`K` (energy)
* :math:`\theta_0` (degrees)

:math:`\theta_0` is specified in degrees, but LAMMPS converts it to
radians internally; hence :math:`K` is effectively energy per
radian\^2.

The required *lj/spica* parameters are extracted automatically from the
pair_style.

Style *sdk*, the original implementation of style *spica*, is available
for backward compatibility.

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

This angle style can only be used if LAMMPS was built with the
CG-SPICA package.  See the :doc:`Build package <Build_package>` doc
page for more info.

Related commands
""""""""""""""""

:doc:`angle_coeff <angle_coeff>`, :doc:`angle_style harmonic <angle_harmonic>`, :doc:`pair_style lj/spica <pair_spica>`,
:doc:`pair_style lj/spica/coul/long <pair_spica>`

Default
"""""""

none
