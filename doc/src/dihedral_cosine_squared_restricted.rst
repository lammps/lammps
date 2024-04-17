.. index:: dihedral_style cosine/squared/restricted

dihedral_style cosine/squared/restricted command
================================================


Syntax
""""""

.. code-block:: LAMMPS

   dihedral_style cosine/squared/restricted

Examples
""""""""

.. code-block:: LAMMPS

   dihedral_style cosine/squared/restricted
   dihedral_coeff 1 10.0 120

Description
"""""""""""

.. versionadded:: 17Apr2024

The *cosine/squared/restricted* dihedral style uses the potential

.. math::

   E = K [\cos(\phi) - \cos(\phi_0)]^2 / \sin^2(\phi)

, which is commonly used in the MARTINI force field.

See :ref:`(Bulacu) <restricted-Bul>` for a description of the restricted dihedral for the MARTINI force field.

The following coefficients must be defined for each dihedral type via the
:doc:`dihedral_coeff <dihedral_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* :math:`K` (energy)
* :math:`\phi_0` (degrees)

:math:`\phi_0` is specified in degrees, but LAMMPS converts it to radians internally.

----------

Restrictions
""""""""""""

This dihedral style can only be used if LAMMPS was built with the
EXTRA-MOLECULE package.  See the :doc:`Build package <Build_package>` doc page
for more info.

Related commands
""""""""""""""""

:doc:`dihedral_coeff <dihedral_coeff>`

Default
"""""""

none

----------

.. _restricted-Bul:

**(Bulacu)** Bulacu, Goga, Zhao, Rossi, Monticelli, Periole, Tieleman, Marrink, J Chem Theory Comput, 9, 3282-3292
(2013).
