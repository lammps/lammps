.. index:: dihedral_style harmonic

dihedral_style harmonic command
===============================

dihedral_style harmonic/intel command
=====================================

dihedral_style harmonic/kk command
==================================

dihedral_style harmonic/omp command
===================================

Syntax
""""""

.. code-block:: LAMMPS

   dihedral_style harmonic

Examples
""""""""

.. code-block:: LAMMPS

   dihedral_style harmonic
   dihedral_coeff 1 80.0 1 2

Description
"""""""""""

The *harmonic* dihedral style uses the potential

.. math::

   E = K [ 1 + d  \cos (n \phi) ]

The following coefficients must be defined for each dihedral type via the
:doc:`dihedral_coeff <dihedral_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* :math:`K` (energy)
* :math:`d` (+1 or -1)
* :math:`n` (integer >= 0)

.. note::

   Here are important points to take note of when defining LAMMPS
   dihedral coefficients for the harmonic style, so that they are
   compatible with how harmonic dihedrals are defined by other force
   fields:

* The LAMMPS convention is that the trans position = 180 degrees, while
  in some force fields trans = 0 degrees.
* Some force fields reverse the sign convention on :math:`d`.
* Some force fields let :math:`n` be positive or negative which corresponds to
  :math:`d = 1` or :math:`d = -1` for the harmonic style.

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

**Default:** none
