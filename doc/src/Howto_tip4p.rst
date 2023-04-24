TIP4P water model
=================

The four-point TIP4P rigid water model extends the traditional
:doc:`three-point TIP3P <Howto_tip3p>` model by adding an additional
site M, usually massless, where the charge associated with the oxygen
atom is placed.  This site M is located at a fixed distance away from
the oxygen along the bisector of the HOH bond angle.  A bond style of
:doc:`harmonic <bond_harmonic>` and an angle style of :doc:`harmonic
<angle_harmonic>` or :doc:`charmm <angle_charmm>` should also be used.
In case of rigid bonds also bond style :doc:`zero <bond_zero>` and angle
style :doc:`zero <angle_zero>` can be used.

There are two ways to implement TIP4P water in LAMMPS:

#. Use a specially written pair style that uses the :ref:`TIP3P geometry
   <tip3p_molecule>` without the point M. The point M location is then
   implicitly derived from the other atoms or each water molecule and
   used during the force computation.  The forces on M are then
   projected on the oxygen and the two hydrogen atoms.  This is
   computationally very efficient, but the charge distribution in space
   is only correct within the tip4p labeled styles.  So all other
   computations using charges will "see" the negative charge incorrectly
   on the oxygen atom.

   This can be done with the following pair styles for Coulomb with a cutoff:

   * :doc:`pair_style tip4p/cut <pair_lj_cut_tip4p>`
   * :doc:`pair_style lj/cut/tip4p/cut <pair_lj_cut_tip4p>`

   or these commands for a long-range Coulomb treatment:

   * :doc:`pair_style tip4p/long <pair_coul>`
   * :doc:`pair_style lj/cut/tip4p/long <pair_lj_cut_tip4p>`
   * :doc:`pair_style lj/long/tip4p/long <pair_lj_long>`
   * :doc:`pair_style tip4p/long/soft <pair_fep_soft>`
   * :doc:`pair_style lj/cut/tip4p/long/soft <pair_fep_soft>`
   * :doc:`kspace_style pppm/tip4p <kspace_style>`
   * :doc:`kspace_style pppm/disp/tip4p <kspace_style>`

   The bond lengths and bond angles should be held fixed using the
   :doc:`fix shake <fix_shake>` or :doc:`fix rattle <fix_shake>` command,
   unless a parameterization for a flexible TIP4P model is used.  The
   parameter sets listed below are all for rigid TIP4P model variants and
   thus the bond and angle force constants are not used and can be set to
   any legal value; only equilibrium length and angle are used.

#. Use an :ref:`explicit 4 point TIP4P geometry <tip4p_molecule>` where
   the oxygen atom carries no charge and the M point no Lennard-Jones
   interactions.  Since :doc:`fix shake <fix_shake>` or :doc:`fix rattle
   <fix_shake>` may not be applied to this kind of geometry, :doc:`fix
   rigid or fix rigid/small <fix_rigid>` or its thermostatted variants
   are required to maintain a rigid geometry.  This avoids some of the
   issues with respect to analysis and non-tip4p styles, but it is a
   more costly force computation (more atoms in the same volume and thus
   more neighbors in the neighbor lists) and requires a much shorter
   timestep for stable integration of the rigid body motion.  Since no
   bonds or angles are required, they do not need to be defined and atom
   style charge would be sufficient for a bulk TIP4P water system.  In
   order to avoid that LAMMPS produces an error due to the massless M
   site a tiny non-zero mass needs to be assigned.

The table below lists the force field parameters (in real :doc:`units
<units>`) to for a selection of popular variants of the TIP4P model.
There is the rigid TIP4P model with a cutoff :ref:`(Jorgensen)
<Jorgensen5>`, the TIP4/Ice model :ref:`(Abascal1) <Abascal1>`, the
TIP4P/2005 model :ref:`(Abascal2) <Abascal2>` and a version of TIP4P
parameters adjusted for use with a long-range Coulombic solver
(e.g. Ewald or PPPM in LAMMPS).  Note that for implicit TIP4P models the
OM distance is specified in the :doc:`pair_style <pair_style>` command,
not as part of the pair coefficients.

   .. list-table::
      :header-rows: 1
      :widths: auto

      * - Parameter
        - TIP4P (original)
        - TIP4P/Ice
        - TIP4P/2005
        - TIP4P (Ewald)
      * - O mass (amu)
        - 15.9994
        - 15.9994
        - 15.9994
        - 15.9994
      * - H mass (amu)
        - 1.008
        - 1.008
        - 1.008
        - 1.008
      * - O or M charge (:math:`e`)
        - -1.040
        - -1.1794
        - -1.1128
        - -1.04844
      * - H charge (:math:`e`)
        - 0.520
        - 0.5897
        - 0.5564
        - 0.52422
      * - LJ :math:`\epsilon` of OO (kcal/mole)
        - 0.1550
        - 0.21084
        - 0.1852
        - 0.16275
      * - LJ :math:`\sigma` of OO (:math:`\AA`)
        - 3.1536
        - 3.1668
        - 3.1589
        - 3.16435
      * - LJ :math:`\epsilon` of HH, MM, OH, OM, HM (kcal/mole)
        - 0.0
        - 0.0
        - 0.0
        - 0.0
      * - LJ :math:`\sigma` of HH, MM, OH, OM, HM (:math:`\AA`)
        - 1.0
        - 1.0
        - 1.0
        - 1.0
      * - :math:`r_0` of OH bond (:math:`\AA`)
        - 0.9572
        - 0.9572
        - 0.9572
        - 0.9572
      * - :math:`\theta_0` of HOH angle
        - 104.52\ :math:`^{\circ}`
        - 104.52\ :math:`^{\circ}`
        - 104.52\ :math:`^{\circ}`
        - 104.52\ :math:`^{\circ}`
      * - OM distance (:math:`\AA`)
        - 0.15
        - 0.1577
        - 0.1546
        - 0.1250

Note that the when using the TIP4P pair style, the neighbor list cutoff
for Coulomb interactions is effectively extended by a distance 2 \* (OM
distance), to account for the offset distance of the fictitious charges
on O atoms in water molecules.  Thus it is typically best in an
efficiency sense to use a LJ cutoff >= Coulomb cutoff + 2\*(OM
distance), to shrink the size of the neighbor list.  This leads to
slightly larger cost for the long-range calculation, so you can test the
trade-off for your model.  The OM distance and the LJ and Coulombic
cutoffs are set in the :doc:`pair_style lj/cut/tip4p/long
<pair_lj_cut_tip4p>` command.

Below is the code for a LAMMPS input file using the implicit method and
the :ref:`TIP3P molecule file <tip3p_molecule>`.  Because the TIP4P
charges are different from TIP3P they need to be reset (or the molecule
file changed):

.. code-block:: LAMMPS

    units real
    atom_style full
    region box block -5 5 -5 5 -5 5
    create_box 2 box bond/types 1 angle/types 1 &
                extra/bond/per/atom 2 extra/angle/per/atom 1 extra/special/per/atom 2

    mass 1 15.9994
    mass 2 1.008

    pair_style lj/cut/tip4p/cut 1 2 1 1 0.15 8.0
    pair_coeff 1 1 0.1550 3.1536
    pair_coeff 2 2 0.0    1.0

    bond_style zero
    bond_coeff 1 0.9574

    angle_style zero
    angle_coeff 1 104.52

    molecule water tip3p.mol  # this uses the TIP3P geometry
    create_atoms 0 random 33 34564 NULL mol water 25367 overlap 1.33
    # must change charges for TIP4P
    set type 1 charge -1.040
    set type 2 charge  0.520

    fix rigid all shake 0.001 10 10000 b 1 a 1
    minimize 0.0 0.0 1000 10000
    run 0 post no

    reset_timestep 0
    velocity all create 300.0 5463576
    fix integrate all nvt temp 300 300 1.0

    thermo_style custom step temp press etotal pe

    thermo 1000
    run 20000
    write_data tip3p.data nocoeff

Below is the code for a LAMMPS input file using the explicit method and
a TIP4P molecule file.  Because of using :doc:`fix rigid/nvt/small
<fix_rigid>` no bonds need to be defined and thus no extra storage needs
to be reserved for them, but we need to switch to atom style full or use
:doc:`fix property/atom mol <fix_property_atom>` so that fix
rigid/nvt/small can identify rigid bodies by their molecule ID:

.. code-block:: LAMMPS

    units real
    atom_style charge
    region box block -5 5 -5 5 -5 5
    create_box 3 box

    mass 1 15.9994
    mass 2 1.008
    mass 3 1.0e-100

    pair_style lj/cut/coul/cut 8.0
    pair_coeff 1 1 0.1550 3.1536
    pair_coeff 2 2 0.0    1.0
    pair_coeff 3 3 0.0    1.0

    fix mol all property/atom mol
    molecule water tip4p.mol
    create_atoms 0 random 33 34564 NULL mol water 25367 overlap 1.33

    timestep 0.1
    fix integrate all rigid/nvt/small molecule temp 300.0 300.0 1.0
    velocity all create 300.0 5463576

    thermo_style custom step temp press etotal density pe ke
    thermo 1000
    run 20000
    write_data tip4p.data nocoeff

.. _tip4p_molecule:
.. code-block::

   # Water molecule. Explicit TIP4P geometry for use with fix rigid

   4 atoms

   Coords

   1    0.00000  -0.06556   0.00000
   2    0.75695   0.52032   0.00000
   3   -0.75695   0.52032   0.00000
   4    0.00000   0.08444   0.00000

   Types

   1        1   # O
   2        2   # H
   3        2   # H
   4        3   # M

   Charges

   1        0.000
   2        0.520
   3        0.520
   4       -1.040


Wikipedia also has a nice article on `water models <https://en.wikipedia.org/wiki/Water_model>`_.

----------

.. _Jorgensen5:

**(Jorgensen)** Jorgensen, Chandrasekhar, Madura, Impey, Klein, J Chem
Phys, 79, 926 (1983).

.. _Abascal1:

**(Abascal1)** Abascal, Sanz, Fernandez, Vega, J Chem Phys, 122, 234511 (2005)
   https://doi.org/10.1063/1.1931662

.. _Abascal2:

**(Abascal2)** Abascal, J Chem Phys, 123, 234505 (2005)
   https://doi.org/10.1063/1.2121687
