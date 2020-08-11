TIP3P water model
=================

The TIP3P water model as implemented in CHARMM
:ref:`(MacKerell) <howto-tip3p>` specifies a 3-site rigid water molecule with
charges and Lennard-Jones parameters assigned to each of the 3 atoms.
In LAMMPS the :doc:`fix shake <fix_shake>` command can be used to hold
the two O-H bonds and the H-O-H angle rigid.  A bond style of
*harmonic* and an angle style of *harmonic* or *charmm* should also be
used.

These are the additional parameters (in real units) to set for O and H
atoms and the water molecule to run a rigid TIP3P-CHARMM model with a
cutoff.  The K values can be used if a flexible TIP3P model (without
fix shake) is desired.  If the LJ epsilon and sigma for HH and OH are
set to 0.0, it corresponds to the original 1983 TIP3P model
:ref:`(Jorgensen) <Jorgensen1>`.

| O mass = 15.9994
| H mass = 1.008
| O charge = -0.834
| H charge = 0.417
| LJ :math:`\epsilon` of OO = 0.1521
| LJ :math:`\sigma` of OO = 3.1507
| LJ :math:`\epsilon` of HH = 0.0460
| LJ :math:`\sigma` of HH = 0.4000
| LJ :math:`\epsilon` of OH = 0.0836
| LJ :math:`\sigma` of OH = 1.7753
| K of OH bond = 450
| :math:`r_0` of OH bond = 0.9572
| K of HOH angle = 55
| :math:`\theta` of HOH angle = 104.52\ :math:`^{\circ}`
|

These are the parameters to use for TIP3P with a long-range Coulomb
solver (e.g. Ewald or PPPM in LAMMPS), see :ref:`(Price) <Price1>` for
details:

| O mass = 15.9994
| H mass = 1.008
| O charge = -0.830
| H charge = 0.415
| LJ :math:`\epsilon` of OO = 0.102
| LJ :math:`\sigma` of OO = 3.188
| LJ :math:`\epsilon`, :math:`\sigma` of OH, HH = 0.0
| K of OH bond = 450
| :math:`r_0` of OH bond = 0.9572
| K of HOH angle = 55
| :math:`\theta` of HOH angle = 104.52\ :math:`^{\circ}`
|

Wikipedia also has a nice article on `water models <http://en.wikipedia.org/wiki/Water_model>`_.

----------

.. _howto-tip3p:

**(MacKerell)** MacKerell, Bashford, Bellott, Dunbrack, Evanseck, Field,
Fischer, Gao, Guo, Ha, et al, J Phys Chem, 102, 3586 (1998).

.. _Jorgensen1:

**(Jorgensen)** Jorgensen, Chandrasekhar, Madura, Impey, Klein, J Chem
Phys, 79, 926 (1983).

.. _Price1:

**(Price)** Price and Brooks, J Chem Phys, 121, 10096 (2004).
