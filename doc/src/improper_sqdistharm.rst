.. index:: improper_style sqdistharm

improper_style sqdistharm command
=================================

Syntax
""""""

.. code-block:: LAMMPS

   improper_style sqdistharm

Examples
""""""""

.. code-block:: LAMMPS

   improper_style sqdistharm
   improper_coeff 1 50.0 0.1

Description
"""""""""""

The *sqdistharm* improper style uses the potential

.. math::

   E = K (d^2 - {d_0}^2)^2

where :math:`d` is the distance between the central atom and the plane formed
by the other three atoms.  If the 4 atoms in an improper quadruplet
(listed in the data file read by the :doc:`read_data <read_data>`
command) are ordered I,J,K,L then the L-atom is assumed to be the
central atom. Note that this is different from the convention used
in the improper_style distance.

The following coefficients must be defined for each improper type via
the improper_coeff command as in the example above, or in the data
file or restart files read by the read_data or read_restart commands:

* :math:`K` (energy/distance\^4)
* :math:`{d_0}^2` (distance\^2)

Note that :math:`{d_0}^2` (in units distance\^2) has be provided and not :math:`d_0`.

----------

Restrictions
""""""""""""

This improper style can only be used if LAMMPS was built with the
MOLECULE package.  See the :doc:`Build package <Build_package>` doc
page for more info.

Related commands
""""""""""""""""

:doc:`improper_coeff <improper_coeff>`

Default
"""""""

none
