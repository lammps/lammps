.. index:: improper_style periodic

improper_style periodic command
===============================

Syntax
""""""

.. code-block:: LAMMPS

   improper_style periodic

Examples
""""""""

.. code-block:: LAMMPS

   improper_style periodic
   improper_coeff 1 10.5 2 180.0

Description
"""""""""""

The *periodic* improper style uses the potential

.. math::

   E = K [1 + \cos (n \phi - \phi_s) ]

where :math:`\phi` is the improper dihedral angle.

This improper style is the LAMMPS equivalent to the GROMACS dihedral 
type 4, the `periodic improper dihedral <https://manual.gromacs.org/documentation/current/reference-manual/functions/bonded-interactions.html#improper-dihedrals-periodic-type>`_.
It enables the correct bookkeeping of topological information and energy
contributions on impropers when moving forth and back between the GROMACS
and LAMMPS ecosystems.

If the 4 atoms in an improper quadruplet (listed in the data file read
by the :doc:`read_data <read_data>` command) are ordered I,J,K,L then
the improper dihedral angle is between the plane of I,J,K and the
plane of J,K,L.  Note that because this is effectively a dihedral
angle, the formula for this improper style is similar to 
:doc:`dihedral_style harmonic <dihedral_harmonic>`, augmented by phase 
shift :math:`\phi_s`.

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
* :math:`n` (integer >=1)
* :math:`\phi_s` (degrees)

----------

Restrictions
""""""""""""

This improper style can only be used if LAMMPS was built with the
MOLECULE-EXTRA package.  See the :doc:`Build package <Build_package>` doc page
for more info.

Related commands
""""""""""""""""

:doc:`improper_coeff <improper_coeff>`

Default
"""""""

none
