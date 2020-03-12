TIP4P water model
=================

The four-point TIP4P rigid water model extends the traditional
three-point TIP3P model by adding an additional site, usually
massless, where the charge associated with the oxygen atom is placed.
This site M is located at a fixed distance away from the oxygen along
the bisector of the HOH bond angle.  A bond style of *harmonic* and an
angle style of *harmonic* or *charmm* should also be used.

A TIP4P model is run with LAMMPS using either this command
for a cutoff model:

:doc:`pair_style lj/cut/tip4p/cut <pair_lj>`

or these two commands for a long-range model:

* :doc:`pair_style lj/cut/tip4p/long <pair_lj>`
* :doc:`kspace_style pppm/tip4p <kspace_style>`

For both models, the bond lengths and bond angles should be held fixed
using the :doc:`fix shake <fix_shake>` command.

These are the additional parameters (in real units) to set for O and H
atoms and the water molecule to run a rigid TIP4P model with a cutoff
:ref:`(Jorgensen) <Jorgensen5>`.  Note that the OM distance is specified in
the :doc:`pair_style <pair_style>` command, not as part of the pair
coefficients.

| O mass = 15.9994
| H mass = 1.008
| O charge = -1.040
| H charge = 0.520
| r0 of OH bond = 0.9572
| theta of HOH angle = 104.52
| OM distance = 0.15
| LJ epsilon of O-O = 0.1550
| LJ sigma of O-O = 3.1536
| LJ epsilon, sigma of OH, HH = 0.0
| Coulombic cutoff = 8.5
| 

For the TIP4/Ice model (J Chem Phys, 122, 234511 (2005);
https://doi.org/10.1063/1.1931662) these values can be used:

| O mass = 15.9994
| H mass =  1.008
| O charge = -1.1794
| H charge =  0.5897
| r0 of OH bond = 0.9572
| theta of HOH angle = 104.52
| OM distance = 0.1577
| LJ epsilon of O-O = 0.21084
| LJ sigma of O-O = 3.1668
| LJ epsilon, sigma of OH, HH = 0.0
| Coulombic cutoff = 8.5
| 

For the TIP4P/2005 model (J Chem Phys, 123, 234505 (2005);
https://doi.org/10.1063/1.2121687), these values can be used:

| O mass = 15.9994
| H mass =  1.008
| O charge = -1.1128
| H charge = 0.5564
| r0 of OH bond = 0.9572
| theta of HOH angle = 104.52
| OM distance = 0.1546
| LJ epsilon of O-O = 0.1852
| LJ sigma of O-O = 3.1589
| LJ epsilon, sigma of OH, HH = 0.0
| Coulombic cutoff = 8.5
| 

These are the parameters to use for TIP4P with a long-range Coulombic
solver (e.g. Ewald or PPPM in LAMMPS):

| O mass = 15.9994
| H mass = 1.008
| O charge = -1.0484
| H charge = 0.5242
| r0 of OH bond = 0.9572
| theta of HOH angle = 104.52
| OM distance = 0.1250
| LJ epsilon of O-O = 0.16275
| LJ sigma of O-O = 3.16435
| LJ epsilon, sigma of OH, HH = 0.0
| 

Note that the when using the TIP4P pair style, the neighbor list
cutoff for Coulomb interactions is effectively extended by a distance
2 \* (OM distance), to account for the offset distance of the
fictitious charges on O atoms in water molecules.  Thus it is
typically best in an efficiency sense to use a LJ cutoff >= Coulomb
cutoff + 2\*(OM distance), to shrink the size of the neighbor list.
This leads to slightly larger cost for the long-range calculation, so
you can test the trade-off for your model.  The OM distance and the LJ
and Coulombic cutoffs are set in the :doc:`pair_style lj/cut/tip4p/long <pair_lj>` command.

Wikipedia also has a nice article on `water models <http://en.wikipedia.org/wiki/Water_model>`_.


----------


.. _Jorgensen5:



**(Jorgensen)** Jorgensen, Chandrasekhar, Madura, Impey, Klein, J Chem
Phys, 79, 926 (1983).
