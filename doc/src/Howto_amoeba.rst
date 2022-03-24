AMOEBA and HIPPO force fields
=============================

The AMOEBA and HIPPO polarizable force fields were developed by Jay
Ponder's group at U Washington at St Louis.  Their implementation in
LAMMPS was done using F90 code provided by the Ponder group from their
`Tinker MD code <https://dasher.wustl.edu/tinker/>`_

NOTE: what version of AMOEBA and HIPPO does LAMMPS implement?

These force fields can be used when polarization effects are desired
in simulations of water, organic molecules, and biomolecules including
proteins, provided that parameterizations (force field files) are
available for the systems you are interested in.  Files in the LAMMPS
potentials directory with a "amoeba" or "hippo" suffix can be used.
The Tinker distribution and website may have other force field files.

Note that currently, HIPPO can only be used for water systems, but
HIPPO files for a variety of small organic and biomolecules are in
preparation by the Ponder group.  Those force field files will be
included in the LAMMPS distribution when available.

The :doc:`pair_style amoeba <pair_amoeba>` doc page gives a brief
description of the AMOEBA and HIPPO force fields.  Further details for
AMOEBA are in these papers: :ref:`(Ren) <amoeba-Ren>`, :ref:`(Shi)
<amoeba-Shi>`.  Further details for HIPPO are in this paper:
:ref:`(Rackers) <amoeba-Rackers>`.

----------

To use the AMOEBA force field in LAMMPS you should use command like
these appropriately in your input script.  The only change needed for
a HIPPO simulation is for the pair_style and pair_coeff commands.  See
examples/amoeba for example input scripts for both AMOEBA and HIPPO.

.. code-block:: LAMMPS

   units              real                           # required
   atom_style         amoeba
   bond_style         class2                         # CLASS2 package
   angle_style        amoeba
   dihedral_style     fourier                        # EXTRA-MOLECULE package
   improper_style     amoeba
                                                     # required per-atom data
   fix                amtype all property/atom i_amtype ghost yes
   fix                extra all property/atom &                      
                      i_amgroup i_ired i_xaxis i_yaxis i_zaxis d_pval ghost yes
   fix                polaxe all property/atom i_polaxe

   fix                pit all amoeba/pitorsion       # PiTorsion terms in FF
   fix_modify         pit energy yes
                                                     # Bitorsion terms in FF
   fix                bit all amoeba/bitorsion bitorsion.ubiquitin.data  
   fix_modify         bit energy yes

   read_data          data.ubiquitin fix amtype NULL "Tinker Types" &
                      fix pit "pitorsion types" "PiTorsion Coeffs" &
                      fix pit pitorsions PiTorsions &
                      fix bit bitorsions BiTorsions

   pair_style         amoeba                          # AMOEBA FF
   pair_coeff         * * amoeba_ubiquitin.prm amoeba_ubiquitin.key

   pair_style         hippo                           # HIPPO FF
   pair_coeff         * * hippo_water.prm hippo_water.key

   special_bonds      lj/coul 0.5 0.5 0.5 one/five yes     # 1-5 neighbors


NOTE: that some options not needed for simpler systems

NOTE: tinker2lmp.py tool

NOTE: are some commands not needed if system is simple and not ubi2 ?

Both AMOEBA and HIPPO use amoeba for bond/angle/dihedral styles,
assuming the molecular system has bonds, angles, or dihedrals.  These
style commands must be be used before the data file is read, as the
data file defines the coefficients for the various bond/angle/dihedral
types.  The pair_style command can come after the read_data command,
as no pair coefficients are defined in the data file.

The :doc:`fix property/atom <fix_property_atom>` commands are required
to create per-atom data used by the AMOEBA/HIPPO force fields, some of
which is defined in the LAMMPS data file.  The LAMMPS data file should
be produced as a pre-processing step by the tinker2lmp.py Python
script in tools/amoeba with a command like one of the following.
The tools/amoeba/README file has more details on using this tool.

.. code-block:: bash

   % python tinker2lmp.py -xyz tinker.xyz -amoeba protein.prm.amoeba -data lmp.data  # AMOEBA for non-periodic systems
   % python tinker2lmp.py -xyz tinker.xyz -amoeba protein.prm.amoeba -data lmp.data -pbc 10.0 20.0 15.0  # AMOEBA for periodic systems
   % python tinker2lmp.py -xyz tinker.xyz -hippo water.prm.hippo -data lmp.data  # HIPPO for non-periodic systems
   % python tinker2lmp.py -xyz tinker.xyz -hippo water.prm.hippo -data lmp.data -pbc 10.0 20.0 15.0  # HIPPO for periodic systems

Note that two input files are needed and a LAMMPS data file (lmp.data)
is produced. The data file will have information on Tinker atom types
and AMOEBA/HIPPO force field parameters for bonds, angles, and
dihedrals.

The first input is an XYZ file listing all the atoms in the
system.

NOTE: is this a Tinker-augmented-XYZ format or standard?  In either
case, how do we suggest LAMMPS users come up with these files?

The second input is an AMOEBA or HIPPO PRM (force field) file.  The
format of these files is defined by Tinker.  A few such files are
provided in the LAMMPS potentials directory.  Others may be available
in the Tinker distribution or from the Ponder group.

The pair_coeff command should specify the same PRM file, and
optionally a Tinker-format KEY file.  See the :doc:`pair_style amoeba
<pair_amoeba>` doc page for more information about Tinker PRM and KEY
files.

Finally, the :doc:`special_bonds <special_bonds>` command is used to
set all LJ and Coulombic 1-2, 1-3, 1-4 weighting factors to non-zero
and non-unity values, and to generate a per-atom list of 1-5 neighbors
as well.  This is to insure all bond-topology neighbors are included
in the neighbor lists used by AMOEBA/HIPPO.  These force fields apply
their own custom weighting factors to all these terms, including the
1-5 neighbors.

----------

These command doc pages have additional details:

* :doc:`pair_style amoeba or hippo <pair_ameoba>`
* :doc:`bond_style amoeba <bond_amoeba>`
* :doc:`angle_style amoeba <angle_charmm>`
* :doc:`dihedral_style amoeba <dihedral_amoeba>`
* :doc:`fix property/atom <fix_property_atom>`
* :doc:`special_bonds <special_bonds>`

----------

.. _howto-Ren:

**(Ren)** Ren and Ponder, J Phys Chem B, 107, 5933 (2003).

.. _howto-Shi:

**(Shi)** Shi, Xiz, Znahg, Best, Wu, Ponder, Ren, J Chem Theory Comp,
 9, 4046, 2013.

.. _howto-Rackers:

**(Rackers)** Rackers and Ponder, J Chem Phys, 150, 084104 (2010).
