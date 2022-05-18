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
|

Note that as originally proposed, the SPC model was run with a 9
Angstrom cutoff for both LJ and Coulomb terms.  It can also be used
with long-range electrostatic solvers (e.g. Ewald or PPPM in LAMMPS)
without changing any of the parameters above, although it becomes
a different model in that mode of usage.

The SPC/E (extended) water model is the same, except
the partial charge assignments change:

| O charge = -0.8476
| H charge = 0.4238
|

See the :ref:`(Berendsen) <howto-Berendsen>` reference for more details on both
the SPC and SPC/E models.

Wikipedia also has a nice article on `water models <http://en.wikipedia.org/wiki/Water_model>`_.

----------

.. _howto-Berendsen:

**(Berendsen)** Berendsen, Grigera, Straatsma, J Phys Chem, 91,
6269-6271 (1987).
