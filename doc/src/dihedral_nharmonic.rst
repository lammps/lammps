.. index:: dihedral_style nharmonic
.. index:: dihedral_style nharmonic/omp

dihedral_style nharmonic command
================================

Accelerator Variants: *nharmonic/omp*

Syntax
""""""

.. code-block:: LAMMPS

   dihedral_style nharmonic

Examples
""""""""

.. code-block:: LAMMPS

   dihedral_style nharmonic
   dihedral_coeff * 3 10.0 20.0 30.0

Description
"""""""""""

The *nharmonic* dihedral style uses the potential:

.. math::

   E = \sum_{i=1,n} A_i  \cos^{i-1}(\phi)

The following coefficients must be defined for each dihedral type via the
:doc:`dihedral_coeff <dihedral_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* :math:`n` (integer >=1)
* :math:`A_1` (energy)
* :math:`A_2` (energy)
* ...
* :math:`A_n` (energy)

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
