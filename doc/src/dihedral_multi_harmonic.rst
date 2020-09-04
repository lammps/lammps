.. index:: dihedral_style multi/harmonic
.. index:: dihedral_style multi/harmonic/omp

dihedral_style multi/harmonic command
=====================================

Accelerator Variants: *multi/harmonic/omp*

Syntax
""""""

.. code-block:: LAMMPS

   dihedral_style multi/harmonic

Examples
""""""""

.. code-block:: LAMMPS

   dihedral_style multi/harmonic
   dihedral_coeff 1 20 20 20 20 20

Description
"""""""""""

The *multi/harmonic* dihedral style uses the potential

.. math::

   E = \sum_{n=1,5} A_n  \cos^{n-1}(\phi)

The following coefficients must be defined for each dihedral type via the
:doc:`dihedral_coeff <dihedral_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* :math:`A_1` (energy)
* :math:`A_2` (energy)
* :math:`A_3` (energy)
* :math:`A_4` (energy)
* :math:`A_5` (energy)

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

This dihedral style can only be used if LAMMPS was built with the
MOLECULE package.  See the :doc:`Build package <Build_package>` doc page
for more info.

Related commands
""""""""""""""""

:doc:`dihedral_coeff <dihedral_coeff>`

Default
"""""""

none
