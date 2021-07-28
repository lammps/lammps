.. index:: dihedral_style quadratic
.. index:: dihedral_style quadratic/omp

dihedral_style quadratic command
================================

Accelerator Variants: *quadratic/omp*

Syntax
""""""

.. code-block:: LAMMPS

   dihedral_style quadratic

Examples
""""""""

.. code-block:: LAMMPS

   dihedral_style quadratic
   dihedral_coeff 100.0 80.0

Description
"""""""""""

The *quadratic* dihedral style uses the potential:

.. math::

   E = K (\phi - \phi_0)^2

This dihedral potential can be used to keep a dihedral in a predefined
value (cis=zero, right-hand convention is used).

The following coefficients must be defined for each dihedral type via
the :doc:`dihedral_coeff <dihedral_coeff>` command as in the example
above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands:

* :math:`K` (energy)
* :math:`\phi_0` (degrees)

:math:`\phi_0` is specified in degrees, but LAMMPS converts it to
radians internally; hence :math:`K` is effectively energy per
radian\^2.

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

This angle style can only be used if LAMMPS was built with the
MOLECULE package.  See the :doc:`Build package <Build_package>` doc
page for more info.

Related commands
""""""""""""""""

:doc:`dihedral_coeff <dihedral_coeff>`

Default
"""""""

none
