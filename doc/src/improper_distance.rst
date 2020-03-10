.. index:: improper_style distance

improper_style distance command
===============================

Syntax
""""""

.. code-block:: LAMMPS

   improper_style distance

Examples
""""""""

.. code-block:: LAMMPS

   improper_style distance
   improper_coeff 1 80.0 100.0

Description
"""""""""""

The *distance* improper style uses the potential

.. math::

   E = K_2 d^2 + K_4 d^4

where :math:`d` is the distance between the central atom and the plane formed
by the other three atoms.  If the 4 atoms in an improper quadruplet
(listed in the data file read by the :doc:`read_data <read_data>`
command) are ordered I,J,K,L then the I-atom is assumed to be the
central atom.

.. image:: JPG/improper_distance.jpg
   :align: center

Note that defining 4 atoms to interact in this way, does not mean that
bonds necessarily exist between I-J, J-K, or K-L, as they would in a
linear dihedral. Normally, the bonds I-J, I-K, I-L would exist for an
improper to be defined between the 4 atoms.

The following coefficients must be defined for each improper type via
the improper\_coeff command as in the example above, or in the data
file or restart files read by the read\_data or read\_restart commands:

* :math:`K_2` (energy/distance\^2)
* :math:`K_4` (energy/distance\^4)

----------

Restrictions
""""""""""""

This improper style can only be used if LAMMPS was built with the
USER-MISC package.  See the :doc:`Build package <Build_package>` doc
page for more info.

Related commands
""""""""""""""""

:doc:`improper_coeff <improper_coeff>`

**Default:** none
