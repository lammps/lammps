.. index:: bond_style quartic
.. index:: bond_style quartic/omp

bond_style quartic command
==========================

Accelerator Variants: *quartic/omp*

Syntax
""""""

.. code-block:: LAMMPS

   bond_style quartic

Examples
""""""""

.. code-block:: LAMMPS

   bond_style quartic
   bond_coeff 2 1200 -0.55 0.25 1.3 34.6878

Description
"""""""""""

The *quartic* bond style uses the potential

.. math::

   E      & = E_q + E_{LJ} \\
   E_q    & = K (r - R_c)^ 2 (r - R_c - B_1) (r - R_c - B_2) + U_0 \\
   E_{LJ} & = \left\{ \begin{array} {l@{\quad:\quad}l}
   4 \epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^6 \right] + \epsilon & r < 2^{\frac{1}{6}}, \epsilon = 1, \sigma = 1 \\
                                                  0 & r >= 2^{\frac{1}{6}}
                         \end{array} \right.

to define a bond that can be broken as the simulation proceeds (e.g.
due to a polymer being stretched).  The :math:`\sigma` and
:math:`\epsilon` used in the LJ portion of the formula are both set
equal to 1.0 by LAMMPS and the LJ portion is cut off at its minimum,
i.e. at :math:`r_c = 2^{\frac{1}{6}}`.

The following coefficients must be defined for each bond type via the
:doc:`bond_coeff <bond_coeff>` command as in the example above, or in
the data file or restart files read by the :doc:`read_data <read_data>`
or :doc:`read_restart <read_restart>` commands:

* :math:`K` (energy/distance\^4)
* :math:`B_1` (distance)
* :math:`B_2` (distance)
* :math:`R_c` (distance)
* :math:`U_0` (energy)

This potential was constructed to mimic the FENE bond potential for
coarse-grained polymer chains.  When monomers with :math:`\sigma =
\epsilon = 1.0` are used, the following choice of parameters gives a
quartic potential that looks nearly like the FENE potential:

.. math::

   K &= 1200 \\
   B_1 &= -0.55 \\
   B_2 &= 0.25 \\
   R_c &= 1.3 \\
   U_0 &= 34.6878

Different parameters can be specified using the :doc:`bond_coeff <bond_coeff>`
command, but you will need to choose them carefully so they form a suitable
bond potential.

:math:`R_c` is the cutoff length at which the bond potential goes smoothly to a
local maximum.  If a bond length ever becomes :math:`> R_c`, LAMMPS "breaks"
the bond, which means two things.  First, the bond potential is turned
off by setting its type to 0, and is no longer computed.  Second, a
pairwise interaction between the two atoms is turned on, since they
are no longer bonded. See the :doc:`Howto <Howto_broken_bonds>` page
on broken bonds for more information.

LAMMPS does the second task via a computational sleight-of-hand.  It
subtracts the pairwise interaction as part of the bond computation.
When the bond breaks, the subtraction stops.  For this to work, the
pairwise interaction must always be computed by the
:doc:`pair_style <pair_style>` command, whether the bond is broken or
not.  This means that :doc:`special_bonds <special_bonds>` must be set
to 1,1,1, as indicated as a restriction below.

Note that when bonds are dumped to a file via the :doc:`dump local <dump>` command, bonds with type 0 are not included.  The
:doc:`delete_bonds <delete_bonds>` command can also be used to query the
status of broken bonds or permanently delete them, e.g.:

.. code-block:: LAMMPS

   delete_bonds all stats
   delete_bonds all bond 0 remove

----------

.. include:: accel_styles.rst

----------

Restrictions
""""""""""""

This bond style can only be used if LAMMPS was built with the MOLECULE
package.  See the :doc:`Build package <Build_package>` page for more
info.

The *quartic* style requires that :doc:`special_bonds <special_bonds>`
parameters be set to 1,1,1.  Three- and four-body interactions (angle,
dihedral, etc) cannot be used with *quartic* bonds.

Related commands
""""""""""""""""

:doc:`bond_coeff <bond_coeff>`, :doc:`delete_bonds <delete_bonds>`

Default
"""""""

none
