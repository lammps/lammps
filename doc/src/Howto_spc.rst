SPC water model
===============

The SPC water model specifies a 3-site rigid water molecule with
charges and Lennard-Jones parameters assigned to each of the 3 atoms.
In LAMMPS the :doc:`fix shake <fix_shake>` command can be used to hold
the two O-H bonds and the H-O-H angle rigid.  A bond style of
*harmonic* and an angle style of *harmonic* or *charmm* should also be
used.

These are the additional parameters (in real units) to set for O and H
atoms and the water molecule to run a rigid SPC model.

| O mass = 15.9994
| H mass = 1.008
| O charge = -0.820
| H charge = 0.410
| LJ :math:`\epsilon` of OO = 0.1553
| LJ :math:`\sigma` of OO = 3.166
| LJ :math:`\epsilon`, :math:`\sigma` of OH, HH = 0.0
| :math:`r_0` of OH bond = 1.0
| :math:`\theta_0` of HOH angle = 109.47\ :math:`^{\circ}`

Note that as originally proposed, the SPC model was run with a 9
Angstrom cutoff for both LJ and Coulomb terms.  It can also be used
with long-range electrostatic solvers (e.g. Ewald or PPPM in LAMMPS)
without changing any of the parameters above, although it becomes
a different model in that mode of usage.

The SPC/E (extended) water model is the same, except
the partial charge assignments change:

| O charge = -0.8476
| H charge = 0.4238

See the :ref:`(Berendsen) <howto-Berendsen>` reference for more details on both
the SPC and SPC/E models.

Below is the code for a LAMMPS input file and a molecule file
(``spce.mol``) of SPC/E water for use with the :doc:`molecule command
<molecule>` demonstrating how to set up a small bulk water system for
SPC/E with rigid bonds.

.. code-block:: LAMMPS

   units real
   atom_style full
   region box block -5 5 -5 5 -5 5
   create_box 2 box  bond/types 1 angle/types 1 &
                   extra/bond/per/atom 2 extra/angle/per/atom 1 extra/special/per/atom 2

   mass 1 15.9994
   mass 2 1.008

   pair_style lj/cut/coul/cut 10.0
   pair_coeff 1 1 0.1553 3.166
   pair_coeff 1 2 0.0    1.0
   pair_coeff 2 2 0.0    1.0

   bond_style zero
   bond_coeff 1 1.0

   angle_style zero
   angle_coeff 1 109.47

   molecule water spce.mol
   create_atoms 0 random 33 34564 NULL mol water 25367 overlap 1.33

   timestep 1.0
   fix rigid     all shake 0.0001 10 10000 b 1 a 1
   minimize 0.0 0.0 1000 10000
   velocity all create 300.0 5463576
   fix integrate all nvt temp 300.0 300.0 100.0

   thermo_style custom step temp press etotal density pe ke
   thermo 1000
   run 20000 upto
   write_data spce.data nocoeff

.. _spce_molecule:
.. code-block::

   # Water molecule. SPC/E geometry

   3 atoms
   2 bonds
   1 angles

   Coords

   1    0.00000  -0.06461   0.00000
   2    0.81649   0.51275   0.00000
   3   -0.81649   0.51275   0.00000

   Types

   1        1   # O
   2        2   # H
   3        2   # H

   Charges

   1       -0.8476
   2        0.4238
   3        0.4238

   Bonds

   1   1      1      2
   2   1      1      3

   Angles

   1   1      2      1      3

   Shake Flags

   1 1
   2 1
   3 1

   Shake Atoms

   1 1 2 3
   2 1 2 3
   3 1 2 3

   Shake Bond Types

   1 1 1 1
   2 1 1 1
   3 1 1 1

   Special Bond Counts

   1 2 0 0
   2 1 1 0
   3 1 1 0

   Special Bonds

   1 2 3
   2 1 3
   3 1 2

Wikipedia also has a nice article on `water models <https://en.wikipedia.org/wiki/Water_model>`_.

----------

.. _howto-Berendsen:

**(Berendsen)** Berendsen, Grigera, Straatsma, J Phys Chem, 91, 6269-6271 (1987).
