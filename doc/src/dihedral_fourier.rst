.. index:: dihedral_style fourier
.. index:: dihedral_style fourier/intel
.. index:: dihedral_style fourier/omp

dihedral_style fourier command
==============================

Accelerator Variants: *fourier/intel*, *fourier/omp*

Syntax
""""""

.. code-block:: LAMMPS

   dihedral_style fourier

Examples
""""""""

.. code-block:: LAMMPS

   dihedral_style fourier
   dihedral_coeff 1 3 -0.846200 3 0.0 7.578800 1 0 0.138000 2 -180.0

Description
"""""""""""

The *fourier* dihedral style uses the potential:

.. math::

   E = \sum_{i=1,m} K_i  [ 1.0 + \cos ( n_i \phi - d_i ) ]

The following coefficients must be defined for each dihedral type via the
:doc:`dihedral_coeff <dihedral_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* :math:`m` (integer >=1)
* :math:`K_1` (energy)
* :math:`n_1` (integer >= 0)
* :math:`d_1` (degrees)
* [...]
* :math:`K_m` (energy)
* :math:`n_m` (integer >= 0)
* :math:`d_m` (degrees)

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

This dihedral style can only be used if LAMMPS was built with the
EXTRA-MOLECULE package.  See the :doc:`Build package <Build_package>` doc
page for more info.

Related commands
""""""""""""""""

:doc:`dihedral_coeff <dihedral_coeff>`

Default
"""""""

none
