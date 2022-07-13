.. index:: improper_style amoeba

improper_style harmonic command
===============================

Syntax
""""""

.. code-block:: LAMMPS

   improper_style amoeba

Examples
""""""""

.. code-block:: LAMMPS

   improper_style amoeba
   improper_coeff 1 49.6

Description
"""""""""""

The *amoeba* improper style uses the potential

.. math::

   E = K (\chi)^2

where :math:`\chi` is the improper angle and :math:`K` is a prefactor.
Note that the usual 1/2 factor is included in :math:`K`.

This formula seems like a simplified version of the formula for the
:doc:`improper_style harmonic <improper_harmonic>` command with
:math:`\chi_0` = 0.0.  However the computation of the angle
:math:`\chi` is done differently to match how the Tinker MD code
computes its out-of-plane improper for the AMOEBA and HIPPO force
fields.  See the :doc:`Howto amoeba <Howto_amoeba>` doc page for more
information about the implementation of AMOEBA and HIPPO in LAMMPS.

If the 4 atoms in an improper quadruplet (listed in the data file read
by the :doc:`read_data <read_data>` command are ordered I,J,K,L then
atoms I,K,L are considered to lie in a plane and atom J is
out-of-place.  The angle :math:`\chi_0` is computed as the Allinger
angle which is defined as the angle between the plane of I,K,L, and
the vector from atom I to atom J.

The following coefficient must be defined for each improper type via
the :doc:`improper_coeff <improper_coeff>` command as in the example
above, or in the data file or restart files read by the
:doc:`read_data <read_data>` or :doc:`read_restart <read_restart>`
commands:

* :math:`K` (energy)

Note that the angle :math:`\chi` is computed in radians; hence
:math:`K` is effectively energy per radian\^2.

----------

Restrictions
""""""""""""

This improper style can only be used if LAMMPS was built with the
AMOEBA package.  See the :doc:`Build package <Build_package>` doc page
for more info.

Related commands
""""""""""""""""

:doc:`improper_coeff <improper_coeff>`, `improper_harmonic
:doc:<improper_harmonic>`

Default
"""""""

none
