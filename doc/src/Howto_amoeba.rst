AMOEBA and HIPPO force fields
=============================

The AMOEBA and HIPPO polarizable force fields were developed by Jay
Ponder's group at the U Washington at St Louis.  The LAMMPS
implementation is based on Fortran 90 code provided by the Ponder
group in their `Tinker MD software <https://dasher.wustl.edu/tinker/>`_.

The current implementation (July 2022) of AMOEBA in LAMMPS matches the
version discussed in :ref:`(Ponder) <amoeba-Ponder>`, :ref:`(Ren)
<amoeba-Ren>`, and :ref:`(Shi) <amoeba-Shi>`.  Likewise the current
implementation of HIPPO in LAMMPS matches the version discussed in
:ref:`(Rackers) <amoeba-Rackers>`.

These force fields can be used when polarization effects are desired
in simulations of water, organic molecules, and biomolecules including
proteins, provided that parameterizations (Tinker PRM force field
files) are available for the systems you are interested in.  Files in
the LAMMPS potentials directory with a "amoeba" or "hippo" suffix can
be used.  The Tinker distribution and website have additional force
field files as well:
`https://github.com/TinkerTools/tinker/tree/release/params
<https://github.com/TinkerTools/tinker/tree/release/params>`_.

Note that currently, HIPPO can only be used for water systems, but
HIPPO files for a variety of small organic and biomolecules are in
preparation by the Ponder group.  Those force field files will be
included in the LAMMPS distribution when available.

To use the AMOEBA or HIPPO force fields, a simulation must be 3d, and
fully periodic or fully non-periodic, and use an orthogonal (not
triclinic) simulation box.

----------

The AMOEBA and HIPPO force fields contain the following terms in their
energy (U) computation.  Further details for AMOEBA equations are in
:ref:`(Ponder) <amoeba-Ponder>`, further details for the HIPPO
equations are in :ref:`(Rackers) <amoeba-Rackers>`.

.. math::

  U & = U_{intermolecular} + U_{intramolecular} \\
  U_{intermolecular} & = U_{hal} + U_{repulsion} + U_{dispersion} + U_{multipole} + U_{polar} + U_{qxfer} \\
  U_{intramolecular} & = U_{bond} + U_{angle} + U_{torsion} + U_{oop} + U_{b\theta} + U_{UB} + U_{pitorsion} + U_{bitorsion}

For intermolecular terms, the AMOEBA force field includes only the
:math:`U_{hal}`, :math:`U_{multipole}`, :math:`U_{polar}` terms.  The
HIPPO force field includes all but the :math:`U_{hal}` term.  In
LAMMPS, these are all computed by the :doc:`pair_style amoeba or hippo
<pair_style>` command.  Note that the :math:`U_{multipole}` and
:math:`U_{polar}` terms in this formula are not the same for the
AMOEBA and HIPPO force fields.

For intramolecular terms, the :math:`U_{bond}`, :math:`U_{angle}`,
:math:`U_{torsion}`, :math:`U_{oop}` terms are computed by the
:doc:`bond_style class2 <bond_class2>` :doc:`angle_style amoeba
<angle_amoeba>`, :doc:`dihedral_style fourier <dihedral_fourier>`, and
:doc:`improper_style amoeba <improper_amoeba>` commands respectively.
The :doc:`angle_style amoeba <angle_amoeba>` command includes the
:math:`U_{b\theta}` bond-angle cross term, and the :math:`U_{UB}` term
for a Urey-Bradley bond contribution between the I,K atoms in the IJK
angle.

The :math:`U_{pitorsion}` term is computed by the :doc:`fix
amoeba/pitorsion <fix_amoeba_pitorsion>` command.  It computes 6-body
interaction between a pair of bonded atoms which each have 2
additional bond partners.

The :math:`U_{bitorsion}` term is computed by the :doc:`fix
amoeba/bitorsion <fix_amoeba_bitorsion>` command.  It computes 5-body
interaction between two 4-body torsions (dihedrals) which overlap,
having 3 atoms in common.

These command doc pages have additional details on the terms they
compute:

* :doc:`pair_style amoeba or hippo <pair_amoeba>`
* :doc:`bond_style class2 <bond_class2>`
* :doc:`angle_style amoeba <angle_amoeba>`
* :doc:`dihedral_style fourier <dihedral_fourier>`
* :doc:`improper_style amoeba <improper_amoeba>`
* :doc:`fix amoeba/pitorsion <fix_amoeba_pitorsion>`
* :doc:`fix amoeba/bitorsion <fix_amoeba_bitorsion>`

----------

To use the AMOEBA or HIPPO force fields in LAMMPS, use commands like
the following appropriately in your input script.  The only change
needed for AMOEBA vs HIPPO simulation is for the :doc:`pair_style
<pair_style>` and :doc:`pair_coeff <pair_coeff>` commands, as shown
below.  See examples/amoeba for example input scripts for both AMOEBA
and HIPPO.

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

The data file read by the :doc:`read_data <read_data>` command should
be created by the tools/tinker/tinker2lmp.py conversion program
described below.  It will create a section in the data file with the
header "Tinker Types".  A :doc:`fix property/atom <fix_property_atom>`
command for the data must be specified before the read_data command.
In the example above the fix ID is *amtype*.

Similarly, if the system you are simulating defines AMOEBA/HIPPO
pitorsion or bitorsion interactions, there will be entries in the data
file for those interactions.  They require a :doc:`fix
amoeba/pitortion <fix_amoeba_pitorsion>` and :doc:`fix
amoeba/bitorsion <fix_amoeba_bitorsion>` command be defined.  In the
example above, the IDs for these two fixes are *pit* and *bit*.

Of course, if the system being modeled does not have one or more of
the following -- bond, angle, dihedral, improper, pitorsion,
bitorsion interactions -- then the corresponding style and fix
commands above do not need to be used.  See the example scripts in
examples/amoeba for water systems as examples; they are simpler than
what is listed above.

The two :doc:`fix property/atom <fix_property_atom>` commands with IDs
(in the example above) *extra* and *polaxe* are also needed to define
internal per-atom quantities used by the AMOEBA and HIPPO force
fields.

The :doc:`pair_coeff <pair_coeff>` command used for either the AMOEBA
or HIPPO force field takes two arguments for Tinker force field files,
namely a PRM and KEY file.  The keyfile can be specified as NULL and
default values for a various settings will be used.  Note that these 2
files are meant to allow use of native Tinker files as-is.  However
LAMMPS does not support all the options which can be included
in a Tinker PRM or KEY file.  See specifics below.

A :doc:`special_bonds <special_bonds>` command with the *one/five*
option is required, since the AMOEBA/HIPPO force fields define
weighting factors for not only 1-2, 1-3, 1-4 interactions, but also
1-5 interactions.  This command will trigger a per-atom list of 1-5
neighbors to be generated.  The AMOEBA and HIPPO force fields define
their own custom weighting factors for all the 1-2, 1-3, 1-4, 1-5
terms which in the Tinker PRM and KEY files; they can be different for
different terms in the force field.

In addition to the list above, these command doc pages have additional
details:

* :doc:`atom_style amoeba <atom_style>`
* :doc:`fix property/atom <fix_property_atom>`
* :doc:`special_bonds <special_bonds>`

----------

Tinker PRM and KEY files

A Tinker PRM file is composed of sections, each of which has multiple
lines.  This is the list of PRM sections LAMMPS knows how to parse and
use.  Any other sections are skipped:

* Angle Bending Parameters
* Atom Type Definitions
* Atomic Multipole Parameters
* Bond Stretching Parameters
* Charge Penetration Parameters
* Charge Transfer Parameters
* Dipole Polarizability Parameters
* Dispersion Parameters
* Force Field Definition
* Literature References
* Out-of-Plane Bend Parameters
* Pauli Repulsion Parameters
* Pi-Torsion Parameters
* Stretch-Bend Parameters
* Torsion-Torsion Parameters
* Torsional Parameters
* Urey-Bradley Parameters
* Van der Waals Pair Parameters
* Van der Waals Parameters

A Tinker KEY file is composed of lines, each of which has a keyword
followed by zero or more parameters.  This is the list of keywords
LAMMPS knows how to parse and use in the same manner Tinker does.  Any
other keywords are skipped.  The value in parenthesis is the default
value for the keyword if it is not specified, or if the keyfile in the
:doc:`pair_coeff <pair_coeff>` command is specified as NULL:

* a-axis (0.0)
* b-axis (0.0)
* c-axis (0.0)
* ctrn-cutoff (6.0)
* ctrn-taper (0.9 * ctrn-cutoff)
* cutoff
* delta-halgren (0.07)
* dewald (no long-range dispersion unless specified)
* dewald-alpha (0.4)
* dewald-cutoff (7.0)
* dispersion-cutoff (9.0)
* dispersion-taper (9.0 * dispersion-cutoff)
* dpme-grid
* dpme-order (4)
* ewald (no long-range electrostatics unless specified)
* ewald-alpha (0.4)
* ewald-cutoff (7.0)
* gamma-halgren (0.12)
* mpole-cutoff (9.0)
* mpole-taper (0.65 * mpole-cutoff)
* pcg-guess (enabled by default)
* pcg-noguess (disable pcg-guess if specified)
* pcg-noprecond (disable pcg-precond if specified)
* pcg-peek (1.0)
* pcg-precond (enabled by default)
* pewald-alpha (0.4)
* pme-grid
* pme-order (5)
* polar-eps (1.0e-6)
* polar-iter (100)
* polar-predict (no prediction operation unless specified)
* ppme-order (5)
* repulsion-cutoff (6.0)
* repulsion-taper (0.9 * repulsion-cutoff)
* taper
* usolve-cutoff (4.5)
* usolve-diag (2.0)
* vdw-cutoff (9.0)
* vdw-taper (0.9 * vdw-cutoff)

----------

Tinker2lmp.py tool

This conversion tool is found in the tools/tinker directory.
As shown in examples/amoeba/README, these commands produce
the data files found in examples/amoeba, and also illustrate
all the options available to use with the tinker2lmp.py script:

.. code-block:: bash

   % python tinker2lmp.py -xyz water_dimer.xyz -amoeba amoeba_water.prm -data data.water_dimer.amoeba                # AMOEBA non-periodic system
   % python tinker2lmp.py -xyz water_dimer.xyz -hippo hippo_water.prm -data data.water_dimer.hippo                   # HIPPO non-periodic system
   % python tinker2lmp.py -xyz water_box.xyz -amoeba amoeba_water.prm -data data.water_box.amoeba -pbc 18.643 18.643 18.643    # AMOEBA periodic system
   % python tinker2lmp.py -xyz water_box.xyz -hippo hippo_water.prm -data data.water_box.hippo -pbc 18.643 18.643 18.643       # HIPPO periodic system
   % python tinker2lmp.py -xyz ubiquitin.xyz -amoeba amoeba_ubiquitin.prm -data data.ubiquitin.new -pbc 54.99 41.91 41.91 -bitorsion bitorsion.ubiquitin.data.new   # system with bitorsions

Switches and their arguments may be specified in any order.

The -xyz switch is required and specifies an input XYZ file as an
argument.  The format of this file is an extended XYZ format defined
and used by Tinker for its input.  Example \*.xyz files are in the
examples/amoeba directory.  The file lists the atoms in the system.
Each atom has the following information: Tinker species name (ignored
by LAMMPS), xyz coordinates, Tinker numeric type, and a list of atom
IDs the atom is bonded to.

Here is more information about the extended XYZ format defined and
used by Tinker, and links to programs that convert standard PDB files
to the extended XYZ format:

* `http://openbabel.org/docs/current/FileFormats/Tinker_XYZ_format.html <http://openbabel.org/docs/current/FileFormats/Tinker_XYZ_format.html>`_
* `https://github.com/emleddin/pdbxyz-xyzpdb <https://github.com/emleddin/pdbxyz-xyzpdb>`_
* `https://github.com/TinkerTools/tinker/blob/release/source/pdbxyz.f <https://github.com/TinkerTools/tinker/blob/release/source/pdbxyz.f>`_

The -amoeba or -hippo switch is required.  It specifies an input
AMOEBA or HIPPO PRM force field file as an argument.  This should be
the same file used by the :doc:`pair_style <pair_style>` command in
the input script.

The -data switch is required.  It specifies an output file name for
the LAMMPS data file that will be produced.

For periodic systems, the -pbc switch is required.  It specifies the
periodic box size for each dimension (x,y,z).  For a Tinker simulation
these are specified in the KEY file.

The -bitorsion switch is only needed if the system contains Tinker
bitorsion interactions.  The data for each type of bitorsion
interaction will be written to the specified file, and read by the
:doc:`fix amoeba/bitorsion <fix_amoeba_bitorsion>` command.  The data
includes 2d arrays of values to which splines are fit, and thus is not
compatible with the LAMMPS data file format.

----------

.. _howto-Ponder:

**(Ponder)** Ponder, Wu, Ren, Pande, Chodera, Schnieders, Haque,  Mobley, Lambrecht, DiStasio Jr, M. Head-Gordon, Clark,  Johnson, T. Head-Gordon, J Phys Chem B, 114, 2549-2564 (2010).

.. _howto-Rackers:

**(Rackers)** Rackers, Silva, Wang, Ponder, J Chem Theory Comput, 17, 7056-7084 (2021).

.. _howto-Ren:

**(Ren)** Ren and Ponder, J Phys Chem B, 107, 5933 (2003).

.. _howto-Shi:

**(Shi)** Shi, Xia, Zhang, Best, Wu, Ponder, Ren, J Chem Theory Comp, 9, 4046, 2013.

