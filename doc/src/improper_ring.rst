.. index:: improper_style ring
.. index:: improper_style ring/omp

improper_style ring command
===========================

Accelerator Variants: *ring/omp*

Syntax
""""""

.. code-block:: LAMMPS

   improper_style ring

Examples
""""""""

.. code-block:: LAMMPS

   improper_style ring
   improper_coeff 1 8000 70.5

Description
"""""""""""

The *ring* improper style uses the potential

.. math::

   E = &\frac{1}{6} K \left(\Delta_{ijl} + \Delta_{ijk} + \Delta_{kjl} \right)^6 \\
   \Delta_{ijl} = & \cos{\theta_{ijl} - \cos{\theta_0}} \\
   \Delta_{ijk} = & \cos{\theta_{ijk} - \cos{\theta_0}} \\
   \Delta_{kjl} = & \cos{\theta_{kjl} - \cos{\theta_0}}

where :math:`K` is a prefactor, :math:`\theta` is the angle formed by
the atoms specified by (i,j,k,l) indices and :math:`\theta_0` its
equilibrium value.

If the 4 atoms in an improper quadruplet (listed in the data file read
by the :doc:`read_data <read_data>` command) are ordered i,j,k,l then
:math:`\theta_{ijl}` is the angle between atoms i,j and l,
:math:`\theta_{ijk}` is the angle between atoms i,j and k,
:math:`\theta_{kjl}` is the angle between atoms j,k, and l.

The "ring" improper style implements the improper potential introduced
by Destree et al., in Equation (9) of :ref:`(Destree) <Destree>`.  This
potential does not affect small amplitude vibrations but is used in an
ad-hoc way to prevent the onset of accidentally large amplitude
fluctuations leading to the occurrence of a planar conformation of the
three bonds i-j, j-k and j-l, an intermediate conformation toward the
chiral inversion of a methine carbon.  In the "Impropers" section of
data file four atoms: i, j, k and l are specified with i,j and l lying
on the backbone of the chain and k specifying the chirality of j.

The following coefficients must be defined for each improper type via
the :doc:`improper_coeff <improper_coeff>` command as in the example
above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands:

* :math:`K` (energy)
* :math:`\theta_0` (degrees)

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

This improper style can only be used if LAMMPS was built with the
MOLECULE package.  See the :doc:`Build package <Build_package>` doc
page for more info.

Related commands
""""""""""""""""

:doc:`improper_coeff <improper_coeff>`

.. _Destree:

**(Destree)** M. Destree, F. Laupretre, A. Lyulin, and J.-P.  Ryckaert,
J Chem Phys, 112, 9632 (2000).
