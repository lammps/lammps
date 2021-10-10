.. index:: bond_style fene
.. index:: bond_style fene/intel
.. index:: bond_style fene/kk
.. index:: bond_style fene/omp

bond_style fene command
=======================

Accelerator Variants: *fene/intel*, *fene/kk*, *fene/omp*

Syntax
""""""

.. code-block:: LAMMPS

   bond_style fene

Examples
""""""""

.. code-block:: LAMMPS

   bond_style fene
   bond_coeff 1 30.0 1.5 1.0 1.0

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

The following coefficients must be defined for each bond type via the
:doc:`bond_coeff <bond_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* :math:`K` (energy/distance\^2)
* :math:`R_0` (distance)
* :math:`\epsilon` (energy)
* :math:`\sigma` (distance)

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

This bond style can only be used if LAMMPS was built with the MOLECULE
package.  See the :doc:`Build package <Build_package>` page for more
info.

You typically should specify :doc:`special_bonds fene <special_bonds>`
or :doc:`special_bonds lj/coul 0 1 1 <special_bonds>` to use this bond
style.  LAMMPS will issue a warning it that's not the case.

Related commands
""""""""""""""""

:doc:`bond_coeff <bond_coeff>`, :doc:`delete_bonds <delete_bonds>`

Default
"""""""

none

----------

.. _fene-Kremer:

**(Kremer)** Kremer, Grest, J Chem Phys, 92, 5057 (1990).
