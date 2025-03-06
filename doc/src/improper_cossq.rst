.. index:: improper_style cossq
.. index:: improper_style cossq/omp

improper_style cossq command
============================

Accelerator Variants: *cossq/omp*

Syntax
""""""

.. code-block:: LAMMPS

   improper_style cossq

Examples
""""""""

.. code-block:: LAMMPS

   improper_style cossq
   improper_coeff 1 4.0 0.0

Description
"""""""""""

The *cossq* improper style uses the potential

.. math::

   E = \frac{1}{2} K \cos^2{\left(\chi - \chi_0\right)}

where :math:`\chi` is the improper angle, :math:`\chi_0` is its
equilibrium value, and :math:`K` is a prefactor.

If the 4 atoms in an improper quadruplet (listed in the data file read
by the :doc:`read_data <read_data>` command) are ordered I,J,K,L then
:math:`\chi` is the angle between the plane of I,J,K and the plane of J,K,L.
Alternatively, you can think of atoms J,K,L as being in a plane, and
atom I above the plane, and :math:`\chi` as a measure of how far
out-of-plane I is with respect to the other 3 atoms.

Note that defining 4 atoms to interact in this way, does not mean that
bonds necessarily exist between I-J, J-K, or K-L, as they would in a
linear dihedral.  Normally, the bonds I-J, I-K, I-L would exist for an
improper to be defined between the 4 atoms.

The following coefficients must be defined for each improper type via
the :doc:`improper_coeff <improper_coeff>` command as in the example
above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands:

* :math:`K` (energy)
* :math:`\chi_0` (degrees)

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

This improper style can only be used if LAMMPS was built with the
EXTRA-MOLECULE package.  See the :doc:`Build package <Build_package>`
doc page for more info.

Related commands
""""""""""""""""

:doc:`improper_coeff <improper_coeff>`

Default
"""""""

none
