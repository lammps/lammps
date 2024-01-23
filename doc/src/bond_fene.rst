.. index:: bond_style fene
.. index:: bond_style fene/nm
.. index:: bond_style fene/intel
.. index:: bond_style fene/kk
.. index:: bond_style fene/omp

bond_style fene command
=======================

Accelerator Variants: *fene/intel*, *fene/kk*, *fene/omp*

bond_style fene/nm command
==========================

Syntax
""""""

.. code-block:: LAMMPS

   bond_style fene
   bond_style fene/nm

Examples
""""""""

.. code-block:: LAMMPS

   bond_style fene
   bond_coeff 1 30.0 1.5 1.0 1.0

   bond_style fene/nm
   bond_coeff 1 2.25344 1.5 1.0 1.12246 2 6

Description
"""""""""""

The *fene* bond style uses the potential

.. math::

   E = -0.5 K R_0^2  \ln \left[ 1 - \left(\frac{r}{R_0}\right)^2\right] + 4 \epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^6 \right] + \epsilon

to define a finite extensible nonlinear elastic (FENE) potential
:ref:`(Kremer) <fene-Kremer>`, used for bead-spring polymer models.  The first
term is attractive, the second Lennard-Jones term is repulsive.  The
first term extends to :math:`R_0`, the maximum extent of the bond.  The second
term is cutoff at :math:`2^\frac{1}{6} \sigma`, the minimum of the LJ potential.

The *fene/nm* bond style substitutes the standard LJ potential with the generalized LJ potential
in the same form as in pair style :doc:`nm/cut <pair_nm>`. The bond energy is then given by

.. math::

  E = -0.5 K R_0^2  \ln \left[ 1 - \left(\frac{r}{R_0}\right)^2\right] + \frac{E_0}{(n-m)} \left[ m \left(\frac{r_0}{r}\right)^n - n \left(\frac{r_0}{r}\right)^m \right]

Similar to the *fene* style, the generalized Lennard-Jones is cut off at
the potential minimum, :math:`r_0`, to be repulsive only.  The following
coefficients must be defined for each bond type via the :doc:`bond_coeff
<bond_coeff>` command as in the example above, or in the data file or
restart files read by the :doc:`read_data <read_data>` or
:doc:`read_restart <read_restart>` commands:

* :math:`K` (energy/distance\^2)
* :math:`R_0` (distance)
* :math:`\epsilon` (energy)
* :math:`\sigma` (distance)

For the *fene/nm* style, the following coefficients are used.  Please
note, that the standard LJ potential and thus the regular FENE potential
is recovered for (n=12 m=6) and :math:`r_0 = 2^\frac{1}{6} \sigma`.

* :math:`K` (energy/distance\^2)
* :math:`R_0` (distance)
* :math:`E_0` (energy)
* :math:`r_0` (distance)
* :math:`n` (unitless)
* :math:`m` (unitless)

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

The *fene* bond style can only be used if LAMMPS was built with the MOLECULE
package; the *fene/nm* bond style can only be used if LAMMPS was built
with the EXTRA-MOLECULE package. See the :doc:`Build package <Build_package>`
page for more info.

You typically should specify :doc:`special_bonds fene <special_bonds>`
or :doc:`special_bonds lj/coul 0 1 1 <special_bonds>` to use this bond
style.  LAMMPS will issue a warning it that's not the case.

Related commands
""""""""""""""""

:doc:`bond_coeff <bond_coeff>`, :doc:`delete_bonds <delete_bonds>`,
:doc:`pair style lj/cut <pair_lj>`, :doc:`pair style nm/cut <pair_nm>`.

Default
"""""""

none

----------

.. _fene-Kremer:

**(Kremer)** Kremer, Grest, J Chem Phys, 92, 5057 (1990).
