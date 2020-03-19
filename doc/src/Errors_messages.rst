Error messages
==============

This is an alphabetic list of the ERROR messages LAMMPS prints out and
the reason why.  If the explanation here is not sufficient, the
documentation for the offending command may help.  Error messages also
list the source file and line number where the error was generated.
For example, a message like this:

.. parsed-literal::

   ERROR: Illegal velocity command (velocity.cpp:78)

means that line #78 in the file src/velocity.cpp generated the error.
Looking in the source code may help you figure out what went wrong.

Note that error messages from :doc:`user-contributed packages <Packages_user>` are not listed here.  If such an error
occurs and is not self-explanatory, you will need to look in the source
code or contact the author of the package.

Doc page with :doc:`WARNING messages <Errors_warnings>`

----------

*1-3 bond count is inconsistent*
   An inconsistency was detected when computing the number of 1-3
   neighbors for each atom.  This likely means something is wrong with
   the bond topologies you have defined.

*1-4 bond count is inconsistent*
   An inconsistency was detected when computing the number of 1-4
   neighbors for each atom.  This likely means something is wrong with
   the bond topologies you have defined.

*Accelerator sharing is not currently supported on system*
   Multiple MPI processes cannot share the accelerator on your
   system. For NVIDIA GPUs, see the nvidia-smi command to change this
   setting.

*All angle coeffs are not set*
   All angle coefficients must be set in the data file or by the
   angle_coeff command before running a simulation.

*All atom IDs = 0 but atom_modify id = yes*
   Self-explanatory.

*All atoms of a swapped type must have same charge.*
   Self-explanatory.

*All atoms of a swapped type must have the same charge.*
   Self-explanatory.

*All bond coeffs are not set*
   All bond coefficients must be set in the data file or by the
   bond_coeff command before running a simulation.

*All dihedral coeffs are not set*
   All dihedral coefficients must be set in the data file or by the
   dihedral_coeff command before running a simulation.

*All improper coeffs are not set*
   All improper coefficients must be set in the data file or by the
   improper_coeff command before running a simulation.

*All masses are not set*
   For atom styles that define masses for each atom type, all masses must
   be set in the data file or by the mass command before running a
   simulation.  They must also be set before using the velocity
   command.

*All mol IDs should be set for fix gcmc group atoms*
   The molecule flag is on, yet not all molecule ids in the fix group
   have been set to non-zero positive values by the user. This is an
   error since all atoms in the fix gcmc group are eligible for deletion,
   rotation, and translation and therefore must have valid molecule ids.

*All pair coeffs are not set*
   All pair coefficients must be set in the data file or by the
   pair_coeff command before running a simulation.

*All read_dump x,y,z fields must be specified for scaled, triclinic coords*
   For triclinic boxes and scaled coordinates you must specify all 3 of
   the x,y,z fields, else LAMMPS cannot reconstruct the unscaled
   coordinates.

*All universe/uloop variables must have same # of values*
   Self-explanatory.

*All variables in next command must be same style*
   Self-explanatory.

*Angle atom missing in delete_bonds*
   The delete_bonds command cannot find one or more atoms in a particular
   angle on a particular processor.  The pairwise cutoff is too short or
   the atoms are too far apart to make a valid angle.

*Angle atom missing in set command*
   The set command cannot find one or more atoms in a particular angle on
   a particular processor.  The pairwise cutoff is too short or the atoms
   are too far apart to make a valid angle.

*Angle atoms %d %d %d missing on proc %d at step %ld*
   One or more of 3 atoms needed to compute a particular angle are
   missing on this processor.  Typically this is because the pairwise
   cutoff is set too short or the angle has blown apart and an atom is
   too far away.

*Angle atoms missing on proc %d at step %ld*
   One or more of 3 atoms needed to compute a particular angle are
   missing on this processor.  Typically this is because the pairwise
   cutoff is set too short or the angle has blown apart and an atom is
   too far away.

*Angle coeff for hybrid has invalid style*
   Angle style hybrid uses another angle style as one of its
   coefficients.  The angle style used in the angle_coeff command or read
   from a restart file is not recognized.

*Angle coeffs are not set*
   No angle coefficients have been assigned in the data file or via the
   angle_coeff command.

*Angle extent > half of periodic box length*
   This error was detected by the neigh_modify check yes setting.  It is
   an error because the angle atoms are so far apart it is ambiguous how
   it should be defined.

*Angle potential must be defined for SHAKE*
   When shaking angles, an angle_style potential must be used.

*Angle style hybrid cannot have hybrid as an argument*
   Self-explanatory.

*Angle style hybrid cannot have none as an argument*
   Self-explanatory.

*Angle style hybrid cannot use same angle style twice*
   Self-explanatory.

*Angle table must range from 0 to 180 degrees*
   Self-explanatory.

*Angle table parameters did not set N*
   List of angle table parameters must include N setting.

*Angle_coeff command before angle_style is defined*
   Coefficients cannot be set in the data file or via the angle_coeff
   command until an angle_style has been assigned.

*Angle_coeff command before simulation box is defined*
   The angle_coeff command cannot be used before a read_data,
   read_restart, or create_box command.

*Angle_coeff command when no angles allowed*
   The chosen atom style does not allow for angles to be defined.

*Angle_style command when no angles allowed*
   The chosen atom style does not allow for angles to be defined.

*Angles assigned incorrectly*
   Angles read in from the data file were not assigned correctly to
   atoms.  This means there is something invalid about the topology
   definitions.

*Angles defined but no angle types*
   The data file header lists angles but no angle types.

*Append boundary must be shrink/minimum*
   The boundary style of the face where atoms are added
   must be of type m (shrink/minimum).

*Arccos of invalid value in variable formula*
   Argument of arccos() must be between -1 and 1.

*Arcsin of invalid value in variable formula*
   Argument of arcsin() must be between -1 and 1.

*Assigning body parameters to non-body atom*
   Self-explanatory.

*Assigning ellipsoid parameters to non-ellipsoid atom*
   Self-explanatory.

*Assigning line parameters to non-line atom*
   Self-explanatory.

*Assigning quat to non-body atom*
   Self-explanatory.

*Assigning tri parameters to non-tri atom*
   Self-explanatory.

*At least one atom of each swapped type must be present to define charges.*
   Self-explanatory.

*Atom IDs must be consecutive for velocity create loop all*
   Self-explanatory.

*Atom IDs must be used for molecular systems*
   Atom IDs are used to identify and find partner atoms in bonds.

*Atom count changed in fix neb*
   This is not allowed in a NEB calculation.

*Atom count is inconsistent, cannot write data file*
   The sum of atoms across processors does not equal the global number
   of atoms.  Probably some atoms have been lost.

*Atom count is inconsistent, cannot write restart file*
   Sum of atoms across processors does not equal initial total count.
   This is probably because you have lost some atoms.

*Atom in too many rigid bodies - boost MAXBODY*
   Fix poems has a parameter MAXBODY (in fix_poems.cpp) which determines
   the maximum number of rigid bodies a single atom can belong to (i.e. a
   multibody joint).  The bodies you have defined exceed this limit.

*Atom sort did not operate correctly*
   This is an internal LAMMPS error.  Please report it to the
   developers.

*Atom style hybrid cannot have hybrid as an argument*
   Self-explanatory.

*Atom style hybrid cannot use same atom style twice*
   Self-explanatory.

*Atom style template molecule must have atom types*
   The defined molecule(s) does not specify atom types.

*Atom style was redefined after using fix property/atom*
   This is not allowed.

*Atom type must be zero in fix gcmc mol command*
   Self-explanatory.

*Atom vector in equal-style variable formula*
   Atom vectors generate one value per atom which is not allowed
   in an equal-style variable.

*Atom-style variable in equal-style variable formula*
   Atom-style variables generate one value per atom which is not allowed
   in an equal-style variable.

*Atom_modify id command after simulation box is defined*
   The atom_modify id command cannot be used after a read_data,
   read_restart, or create_box command.

*Atom_modify map command after simulation box is defined*
   The atom_modify map command cannot be used after a read_data,
   read_restart, or create_box command.

*Atom_modify sort and first options cannot be used together*
   Self-explanatory.

*Atom_style command after simulation box is defined*
   The atom_style command cannot be used after a read_data,
   read_restart, or create_box command.

*Atom_style line can only be used in 2d simulations*
   Self-explanatory.

*Atom_style tri can only be used in 3d simulations*
   Self-explanatory.

*Atomfile variable could not read values*
   Check the file assigned to the variable.

*Atomfile variable in equal-style variable formula*
   Self-explanatory.

*Atomfile-style variable in equal-style variable formula*
   Self-explanatory.

*Attempt to pop empty stack in fix box/relax*
   Internal LAMMPS error.  Please report it to the developers.

*Attempt to push beyond stack limit in fix box/relax*
   Internal LAMMPS error.  Please report it to the developers.

*Attempting to rescale a 0.0 temperature*
   Cannot rescale a temperature that is already 0.0.

*Attempting to insert more particles than available lattice points*
   Self-explanatory.

*Bad FENE bond*
   Two atoms in a FENE bond have become so far apart that the bond cannot
   be computed.

*Bad TIP4P angle type for PPPM/TIP4P*
   Specified angle type is not valid.

*Bad TIP4P angle type for PPPMDisp/TIP4P*
   Specified angle type is not valid.

*Bad TIP4P bond type for PPPM/TIP4P*
   Specified bond type is not valid.

*Bad TIP4P bond type for PPPMDisp/TIP4P*
   Specified bond type is not valid.

*Bad fix ID in fix append/atoms command*
   The value of the fix_id for keyword spatial must start with "f\_".

*Bad grid of processors*
   The 3d grid of processors defined by the processors command does not
   match the number of processors LAMMPS is being run on.

*Bad kspace_modify kmax/ewald parameter*
   Kspace_modify values for the kmax/ewald keyword must be integers > 0

*Bad kspace_modify slab parameter*
   Kspace_modify value for the slab/volume keyword must be >= 2.0.

*Bad matrix inversion in mldivide3*
   This error should not occur unless the matrix is badly formed.

*Bad principal moments*
   Fix rigid did not compute the principal moments of inertia of a rigid
   group of atoms correctly.

*Bad quadratic solve for particle/line collision*
   This is an internal error.  It should normally not occur.

*Bad quadratic solve for particle/tri collision*
   This is an internal error.  It should normally not occur.

*Bad real space Coulombic cutoff in fix tune/kspace*
   Fix tune/kspace tried to find the optimal real space Coulombic cutoff using
   the Newton-Rhaphson method, but found a non-positive or NaN cutoff

*Balance command before simulation box is defined*
   The balance command cannot be used before a read_data, read_restart,
   or create_box command.

*Balance produced bad splits*
   This should not occur.  It means two or more cutting plane locations
   are on top of each other or out of order.  Report the problem to the
   developers.

*Balance rcb cannot be used with comm_style brick*
   Comm_style tiled must be used instead.

*Balance shift string is invalid*
   The string can only contain the characters "x", "y", or "z".

*Bias compute does not calculate a velocity bias*
   The specified compute must compute a bias for temperature.

*Bias compute does not calculate temperature*
   The specified compute must compute temperature.

*Bias compute group does not match compute group*
   The specified compute must operate on the same group as the parent
   compute.

*Big particle in fix srd cannot be point particle*
   Big particles must be extended spheroids or ellipsoids.

*Bigint setting in lmptype.h is invalid*
   Size of bigint is less than size of tagint.

*Bigint setting in lmptype.h is not compatible*
   Format of bigint stored in restart file is not consistent with LAMMPS
   version you are running.  See the settings in src/lmptype.h

*Bitmapped lookup tables require int/float be same size*
   Cannot use pair tables on this machine, because of word sizes.  Use
   the pair_modify command with table 0 instead.

*Bitmapped table in file does not match requested table*
   Setting for bitmapped table in pair_coeff command must match table
   in file exactly.

*Bitmapped table is incorrect length in table file*
   Number of table entries is not a correct power of 2.

*Bond and angle potentials must be defined for TIP4P*
   Cannot use TIP4P pair potential unless bond and angle potentials
   are defined.

*Bond atom missing in box size check*
   The 2nd atoms needed to compute a particular bond is missing on this
   processor.  Typically this is because the pairwise cutoff is set too
   short or the bond has blown apart and an atom is too far away.

*Bond atom missing in delete_bonds*
   The delete_bonds command cannot find one or more atoms in a particular
   bond on a particular processor.  The pairwise cutoff is too short or
   the atoms are too far apart to make a valid bond.

*Bond atom missing in image check*
   The 2nd atom in a particular bond is missing on this processor.
   Typically this is because the pairwise cutoff is set too short or the
   bond has blown apart and an atom is too far away.

*Bond atom missing in set command*
   The set command cannot find one or more atoms in a particular bond on
   a particular processor.  The pairwise cutoff is too short or the atoms
   are too far apart to make a valid bond.

*Bond atoms %d %d missing on proc %d at step %ld*
   The 2nd atom needed to compute a particular bond is missing on this
   processor.  Typically this is because the pairwise cutoff is set too
   short or the bond has blown apart and an atom is too far away.

*Bond atoms missing on proc %d at step %ld*
   The 2nd atom needed to compute a particular bond is missing on this
   processor.  Typically this is because the pairwise cutoff is set too
   short or the bond has blown apart and an atom is too far away.

*Bond coeff for hybrid has invalid style*
   Bond style hybrid uses another bond style as one of its coefficients.
   The bond style used in the bond_coeff command or read from a restart
   file is not recognized.

*Bond coeffs are not set*
   No bond coefficients have been assigned in the data file or via the
   bond_coeff command.

*Bond extent > half of periodic box length*
   This error was detected by the neigh_modify check yes setting.  It is
   an error because the bond atoms are so far apart it is ambiguous how
   it should be defined.

*Bond potential must be defined for SHAKE*
   Cannot use fix shake unless bond potential is defined.

*Bond style hybrid cannot have hybrid as an argument*
   Self-explanatory.

*Bond style hybrid cannot have none as an argument*
   Self-explanatory.

*Bond style hybrid cannot use same bond style twice*
   Self-explanatory.

*Bond style quartic cannot be used with 3,4-body interactions*
   No angle, dihedral, or improper styles can be defined when using
   bond style quartic.

*Bond style quartic cannot be used with atom style template*
   This bond style can change the bond topology which is not
   allowed with this atom style.

*Bond style quartic requires special_bonds = 1,1,1*
   This is a restriction of the current bond quartic implementation.

*Bond table parameters did not set N*
   List of bond table parameters must include N setting.

*Bond table values are not increasing*
   The values in the tabulated file must be monotonically increasing.

*BondAngle coeff for hybrid angle has invalid format*
   No "ba" field should appear in data file entry.

*BondBond coeff for hybrid angle has invalid format*
   No "bb" field should appear in data file entry.

*Bond_coeff command before bond_style is defined*
   Coefficients cannot be set in the data file or via the bond_coeff
   command until an bond_style has been assigned.

*Bond_coeff command before simulation box is defined*
   The bond_coeff command cannot be used before a read_data,
   read_restart, or create_box command.

*Bond_coeff command when no bonds allowed*
   The chosen atom style does not allow for bonds to be defined.

*Bond_style command when no bonds allowed*
   The chosen atom style does not allow for bonds to be defined.

*Bonds assigned incorrectly*
   Bonds read in from the data file were not assigned correctly to atoms.
   This means there is something invalid about the topology definitions.

*Bonds defined but no bond types*
   The data file header lists bonds but no bond types.

*Bond/react: Cannot use fix bond/react with non-molecular systems*
   Only systems with bonds that can be changed can be used. Atom_style
   template does not qualify.

*Bond/react: Invalid template atom ID in map file*
   Atom IDs in molecule templates range from 1 to the number of atoms in the template.

*Bond/react: Rmax cutoff is longer than pairwise cutoff*
   This is not allowed because bond creation is done using the pairwise
   neighbor list.

*Bond/react: Molecule template ID for fix bond/react does not exist*
   A valid molecule template must have been created with the molecule
   command.

*Bond/react: Reaction templates must contain the same number of atoms*
   There should be a one-to-one correspondence between atoms in the
   pre-reacted and post-reacted templates, as specified by the map file.

*Bond/react: Unknown section in map file*
   Please ensure reaction map files are properly formatted.

*Bond/react: Atom affected by reaction too close to template edge*
   This means an atom which changes type or connectivity during the
   reaction is too close to an 'edge' atom defined in the map
   file. This could cause incorrect assignment of bonds, angle, etc.
   Generally, this means you must include more atoms in your templates,
   such that there are at least two atoms between each atom involved in
   the reaction and an edge atom.

*Bond/react: Fix bond/react needs ghost atoms from farther away*
   This is because a processor needs to map the entire unreacted
   molecule template onto simulation atoms it knows about. The
   comm_modify cutoff command can be used to extend the communication
   range.

*Bond/react: A deleted atom cannot remain bonded to an atom that is not deleted*
   Self-explanatory.

*Bond/react: First neighbors of chiral atoms must be of mutually different types*
   Self-explanatory.

*Bond/react: Chiral atoms must have exactly four first neighbors*
   Self-explanatory.

*Bond/react: Molecule template 'Coords' section required for chiralIDs keyword*
   The coordinates of atoms in the pre-reacted template are used to determine
   chirality.

*Bond/react special bond generation overflow*
   The number of special bonds per-atom created by a reaction exceeds the
   system setting. See the read_data or create_box command for how to
   specify this value.

*Bond/react topology/atom exceed system topology/atom*
   The number of bonds, angles etc per-atom created by a reaction exceeds
   the system setting. See the read_data or create_box command for how to
   specify this value.

*Both restart files must use % or neither*
   Self-explanatory.

*Both restart files must use MPI-IO or neither*
   Self-explanatory.

*Both sides of boundary must be periodic*
   Cannot specify a boundary as periodic only on the lo or hi side.  Must
   be periodic on both sides.

*Boundary command after simulation box is defined*
   The boundary command cannot be used after a read_data, read_restart,
   or create_box command.

*Box bounds are invalid*
   The box boundaries specified in the read_data file are invalid.  The
   lo value must be less than the hi value for all 3 dimensions.

*Box command after simulation box is defined*
   The box command cannot be used after a read_data, read_restart, or
   create_box command.

*CPU neighbor lists must be used for ellipsoid/sphere mix.*
   When using Gay-Berne or RE-squared pair styles with both ellipsoidal and
   spherical particles, the neighbor list must be built on the CPU

*Can not specify Pxy/Pxz/Pyz in fix box/relax with non-triclinic box*
   Only triclinic boxes can be used with off-diagonal pressure components.
   See the region prism command for details.

*Can not specify Pxy/Pxz/Pyz in fix nvt/npt/nph with non-triclinic box*
   Only triclinic boxes can be used with off-diagonal pressure components.
   See the region prism command for details.

*Can only use -plog with multiple partitions*
   Self-explanatory.  See doc page discussion of command-line switches.

*Can only use -pscreen with multiple partitions*
   Self-explanatory.  See doc page discussion of command-line switches.

*Can only use Kokkos supported regions with Kokkos package*
   Self-explanatory.

*Can only use NEB with 1-processor replicas*
   This is current restriction for NEB as implemented in LAMMPS.

*Can only use TAD with 1-processor replicas for NEB*
   This is current restriction for NEB as implemented in LAMMPS.

*Cannot (yet) do analytic differentiation with pppm/gpu*
   This is a current restriction of this command.

*Cannot (yet) request ghost atoms with Kokkos half neighbor list*
   This feature is not yet supported.

*Cannot (yet) use 'electron' units with dipoles*
   This feature is not yet supported.

*Cannot (yet) use Ewald with triclinic box and slab correction*
   This feature is not yet supported.

*Cannot (yet) use K-space slab correction with compute group/group for triclinic systems*
   This option is not yet supported.

*Cannot (yet) use MSM with 2d simulation*
   This feature is not yet supported.

*Cannot (yet) use PPPM with triclinic box and TIP4P*
   This feature is not yet supported.

*Cannot (yet) use PPPM with triclinic box and kspace_modify diff ad*
   This feature is not yet supported.

*Cannot (yet) use PPPM with triclinic box and slab correction*
   This feature is not yet supported.

*Cannot (yet) use kspace slab correction with long-range dipoles and non-neutral systems or per-atom energy*
   This feature is not yet supported.

*Cannot (yet) use kspace_modify diff ad with compute group/group*
   This option is not yet supported.

*Cannot (yet) use kspace_style pppm/stagger with triclinic systems*
   This feature is not yet supported.

*Cannot (yet) use molecular templates with Kokkos*
   Self-explanatory.

*Cannot (yet) use respa with Kokkos*
   Self-explanatory.

*Cannot (yet) use rigid bodies with fix deform and Kokkos*
   Self-explanatory.

*Cannot (yet) use rigid bodies with fix nh and Kokkos*
   Self-explanatory.

*Cannot (yet) use single precision with MSM (remove -DFFT_SINGLE from Makefile and re-compile)*
   Single precision cannot be used with MSM.

*Cannot add atoms to fix move variable*
   Atoms can not be added afterwards to this fix option.

*Cannot append atoms to a triclinic box*
   The simulation box must be defined with edges aligned with the
   Cartesian axes.

*Cannot balance in z dimension for 2d simulation*
   Self-explanatory.

*Cannot change box ortho/triclinic with certain fixes defined*
   This is because those fixes store the shape of the box.  You need to
   use unfix to discard the fix, change the box, then redefine a new
   fix.

*Cannot change box ortho/triclinic with dumps defined*
   This is because some dumps store the shape of the box.  You need to
   use undump to discard the dump, change the box, then redefine a new
   dump.

*Cannot change box tilt factors for orthogonal box*
   Cannot use tilt factors unless the simulation box is non-orthogonal.

*Cannot change box to orthogonal when tilt is non-zero*
   Self-explanatory.

*Cannot change box z boundary to non-periodic for a 2d simulation*
   Self-explanatory.

*Cannot change dump_modify every for dump dcd*
   The frequency of writing dump dcd snapshots cannot be changed.

*Cannot change dump_modify every for dump xtc*
   The frequency of writing dump xtc snapshots cannot be changed.

*Cannot change timestep once fix srd is setup*
   This is because various SRD properties depend on the timestep
   size.

*Cannot change timestep with fix pour*
   This is because fix pour pre-computes the time delay for particles to
   fall out of the insertion volume due to gravity.

*Cannot change to comm_style brick from tiled layout*
   Self-explanatory.

*Cannot change_box after reading restart file with per-atom info*
   This is because the restart file info cannot be migrated with the
   atoms.  You can get around this by performing a 0-timestep run which
   will assign the restart file info to actual atoms.

*Cannot change_box in xz or yz for 2d simulation*
   Self-explanatory.

*Cannot change_box in z dimension for 2d simulation*
   Self-explanatory.

*Cannot clear group all*
   This operation is not allowed.

*Cannot close restart file - MPI error: %s*
   This error was generated by MPI when reading/writing an MPI-IO restart
   file.

*Cannot compute initial g_ewald_disp*
   LAMMPS failed to compute an initial guess for the PPPM_disp g_ewald_6
   factor that partitions the computation between real space and k-space
   for Dispersion interactions.

*Cannot create an atom map unless atoms have IDs*
   The simulation requires a mapping from global atom IDs to local atoms,
   but the atoms that have been defined have no IDs.

*Cannot create atoms with undefined lattice*
   Must use the lattice command before using the create_atoms
   command.

*Cannot create/grow a vector/array of pointers for %s*
   LAMMPS code is making an illegal call to the templated memory
   allocaters, to create a vector or array of pointers.

*Cannot create_atoms after reading restart file with per-atom info*
   The per-atom info was stored to be used when by a fix that you may
   re-define.  If you add atoms before re-defining the fix, then there
   will not be a correct amount of per-atom info.

*Cannot create_box after simulation box is defined*
   A simulation box can only be defined once.

*Cannot currently use pair reax with pair hybrid*
   This is not yet supported.

*Cannot currently use pppm/gpu with fix balance.*
   Self-explanatory.

*Cannot delete group all*
   Self-explanatory.

*Cannot delete group currently used by a compute*
   Self-explanatory.

*Cannot delete group currently used by a dump*
   Self-explanatory.

*Cannot delete group currently used by a fix*
   Self-explanatory.

*Cannot delete group currently used by atom_modify first*
   Self-explanatory.

*Cannot delete_atoms bond yes for non-molecular systems*
   Self-explanatory.

*Cannot displace_atoms after reading restart file with per-atom info*
   This is because the restart file info cannot be migrated with the
   atoms.  You can get around this by performing a 0-timestep run which
   will assign the restart file info to actual atoms.

*Cannot do GCMC on atoms in atom_modify first group*
   This is a restriction due to the way atoms are organized in a list to
   enable the atom_modify first command.

*Cannot do atom/swap on atoms in atom_modify first group*
   This is a restriction due to the way atoms are organized in a list to
   enable the atom_modify first command.

*Cannot dump sort on atom IDs with no atom IDs defined*
   Self-explanatory.

*Cannot dump sort when multiple dump files are written*
   In this mode, each processor dumps its atoms to a file, so
   no sorting is allowed.

*Cannot embed Python when also extending Python with LAMMPS*
   When running LAMMPS via Python through the LAMMPS library interface
   you cannot also user the input script python command.

*Cannot evaporate atoms in atom_modify first group*
   This is a restriction due to the way atoms are organized in
   a list to enable the atom_modify first command.

*Cannot find create_bonds group ID*
   Self-explanatory.

*Cannot find delete_bonds group ID*
   Group ID used in the delete_bonds command does not exist.

*Cannot find specified group ID for core particles*
   Self-explanatory.

*Cannot find specified group ID for shell particles*
   Self-explanatory.

*Cannot have both pair_modify shift and tail set to yes*
   These 2 options are contradictory.

*Cannot intersect groups using a dynamic group*
   This operation is not allowed.

*Cannot mix molecular and molecule template atom styles*
   Self-explanatory.

*Cannot open -reorder file*
   Self-explanatory.

*Cannot open ADP potential file %s*
   The specified ADP potential file cannot be opened.  Check that the
   path and name are correct.

*Cannot open AIREBO potential file %s*
   The specified AIREBO potential file cannot be opened.  Check that the
   path and name are correct.

*Cannot open BOP potential file %s*
   The specified BOP potential file cannot be opened.  Check that the
   path and name are correct.

*Cannot open COMB potential file %s*
   The specified COMB potential file cannot be opened.  Check that the
   path and name are correct.

*Cannot open COMB3 lib.comb3 file*
   The COMB3 library file cannot be opened.  Check that the path and name
   are correct.

*Cannot open COMB3 potential file %s*
   The specified COMB3 potential file cannot be opened.  Check that the
   path and name are correct.

*Cannot open EAM potential file %s*
   The specified EAM potential file cannot be opened.  Check that the
   path and name are correct.

*Cannot open EIM potential file %s*
   The specified EIM potential file cannot be opened.  Check that the
   path and name are correct.

*Cannot open LCBOP potential file %s*
   The specified LCBOP potential file cannot be opened.  Check that the
   path and name are correct.

*Cannot open MEAM potential file %s*
   The specified MEAM potential file cannot be opened.  Check that the
   path and name are correct.

*Cannot open SNAP coefficient file %s*
   The specified SNAP coefficient file cannot be opened.  Check that the
   path and name are correct.

*Cannot open SNAP parameter file %s*
   The specified SNAP parameter file cannot be opened.  Check that the
   path and name are correct.

*Cannot open Stillinger-Weber potential file %s*
   The specified SW potential file cannot be opened.  Check that the path
   and name are correct.

*Cannot open Tersoff potential file %s*
   The specified potential file cannot be opened.  Check that the path
   and name are correct.

*Cannot open Vashishta potential file %s*
   The specified Vashishta potential file cannot be opened.  Check that the path
   and name are correct.

*Cannot open balance output file*
   Self-explanatory.

*Cannot open coul/streitz potential file %s*
   The specified coul/streitz potential file cannot be opened.  Check
   that the path and name are correct.

*Cannot open custom file*
   Self-explanatory.

*Cannot open data file %s*
   The specified file cannot be opened.  Check that the path and name are
   correct.

*Cannot open dir to search for restart file*
   Using a "\*" in the name of the restart file will open the current
   directory to search for matching file names.

*Cannot open dump file*
   Self-explanatory.

*Cannot open dump file %s*
   The output file for the dump command cannot be opened.  Check that the
   path and name are correct.

*Cannot open file %s*
   The specified file cannot be opened.  Check that the path and name are
   correct. If the file is a compressed file, also check that the gzip
   executable can be found and run.

*Cannot open file variable file %s*
   The specified file cannot be opened.  Check that the path and name are
   correct.

*Cannot open fix ave/chunk file %s*
   The specified file cannot be opened.  Check that the path and name are
   correct.

*Cannot open fix ave/correlate file %s*
   The specified file cannot be opened.  Check that the path and name are
   correct.

*Cannot open fix ave/histo file %s*
   The specified file cannot be opened.  Check that the path and name are
   correct.

*Cannot open fix ave/time file %s*
   The specified file cannot be opened.  Check that the path and name are
   correct.

*Cannot open fix balance output file*
   Self-explanatory.

*Cannot open fix poems file %s*
   The specified file cannot be opened.  Check that the path and name are
   correct.

*Cannot open fix print file %s*
   The output file generated by the fix print command cannot be opened

*Cannot open fix qeq parameter file %s*
   The specified file cannot be opened.  Check that the path and name are
   correct.

*Cannot open fix qeq/comb file %s*
   The output file for the fix qeq/combs command cannot be opened.
   Check that the path and name are correct.

*Cannot open fix reax/bonds file %s*
   The output file for the fix reax/bonds command cannot be opened.
   Check that the path and name are correct.

*Cannot open fix rigid infile %s*
   The specified file cannot be opened.  Check that the path and name are
   correct.

*Cannot open fix rigid restart file %s*
   The specified file cannot be opened.  Check that the path and name are
   correct.

*Cannot open fix rigid/small infile %s*
   The specified file cannot be opened.  Check that the path and name are
   correct.

*Cannot open fix tmd file %s*
   The output file for the fix tmd command cannot be opened.  Check that
   the path and name are correct.

*Cannot open fix ttm file %s*
   The output file for the fix ttm command cannot be opened.  Check that
   the path and name are correct.

*Cannot open gzipped file*
   LAMMPS was compiled without support for reading and writing gzipped
   files through a pipeline to the gzip program with -DLAMMPS_GZIP.

*Cannot open input script %s*
   Self-explanatory.

*Cannot open log.cite file*
   This file is created when you use some LAMMPS features, to indicate
   what paper you should cite on behalf of those who implemented
   the feature.  Check that you have write privileges into the directory
   you are running in.

*Cannot open log.lammps for writing*
   The default LAMMPS log file cannot be opened.  Check that the
   directory you are running in allows for files to be created.

*Cannot open logfile*
   The LAMMPS log file named in a command-line argument cannot be opened.
   Check that the path and name are correct.

*Cannot open logfile %s*
   The LAMMPS log file specified in the input script cannot be opened.
   Check that the path and name are correct.

*Cannot open molecule file %s*
   The specified file cannot be opened.  Check that the path and name are
   correct.

*Cannot open nb3b/harmonic potential file %s*
   The specified potential file cannot be opened.  Check that the path
   and name are correct.

*Cannot open pair_write file*
   The specified output file for pair energies and forces cannot be
   opened.  Check that the path and name are correct.

*Cannot open polymorphic potential file %s*
   The specified polymorphic potential file cannot be opened.  Check that
   the path and name are correct.

*Cannot open print file %s*
   Self-explanatory.

*Cannot open processors output file*
   Self-explanatory.

*Cannot open restart file %s*
   Self-explanatory.

*Cannot open restart file for reading - MPI error: %s*
   This error was generated by MPI when reading/writing an MPI-IO restart
   file.

*Cannot open restart file for writing - MPI error: %s*
   This error was generated by MPI when reading/writing an MPI-IO restart
   file.

*Cannot open screen file*
   The screen file specified as a command-line argument cannot be
   opened.  Check that the directory you are running in allows for files
   to be created.

*Cannot open temporary file for world counter.*
   Self-explanatory.

*Cannot open universe log file*
   For a multi-partition run, the master log file cannot be opened.
   Check that the directory you are running in allows for files to be
   created.

*Cannot open universe screen file*
   For a multi-partition run, the master screen file cannot be opened.
   Check that the directory you are running in allows for files to be
   created.

*Cannot read from restart file - MPI error: %s*
   This error was generated by MPI when reading/writing an MPI-IO restart
   file.

*Cannot read_data without add keyword after simulation box is defined*
   Self-explanatory.

*Cannot read_restart after simulation box is defined*
   The read_restart command cannot be used after a read_data,
   read_restart, or create_box command.

*Cannot redefine variable as a different style*
   An equal-style variable can be re-defined but only if it was
   originally an equal-style variable.

*Cannot replicate 2d simulation in z dimension*
   The replicate command cannot replicate a 2d simulation in the z
   dimension.

*Cannot replicate with fixes that store atom quantities*
   Either fixes are defined that create and store atom-based vectors or a
   restart file was read which included atom-based vectors for fixes.
   The replicate command cannot duplicate that information for new atoms.
   You should use the replicate command before fixes are applied to the
   system.

*Cannot reset timestep with a dynamic region defined*
   Dynamic regions (see the region command) have a time dependence.
   Thus you cannot change the timestep when one or more of these
   are defined.

*Cannot reset timestep with a time-dependent fix defined*
   You cannot reset the timestep when a fix that keeps track of elapsed
   time is in place.

*Cannot run 2d simulation with non-periodic Z dimension*
   Use the boundary command to make the z dimension periodic in order to
   run a 2d simulation.

*Cannot set bond topology types for atom style template*
   The bond, angle, etc types cannot be changed for this atom style since
   they are static settings in the molecule template files.

*Cannot set both respa pair and inner/middle/outer*
   In the rRESPA integrator, you must compute pairwise potentials either
   all together (pair), or in pieces (inner/middle/outer).  You can't do
   both.

*Cannot set cutoff/multi before simulation box is defined*
   Self-explanatory.

*Cannot set dpd/theta for this atom style*
   Self-explanatory.

*Cannot set dump_modify flush for dump xtc*
   Self-explanatory.

*Cannot set mass for this atom style*
   This atom style does not support mass settings for each atom type.
   Instead they are defined on a per-atom basis in the data file.

*Cannot set meso/cv for this atom style*
   Self-explanatory.

*Cannot set meso/e for this atom style*
   Self-explanatory.

*Cannot set meso/rho for this atom style*
   Self-explanatory.

*Cannot set non-zero image flag for non-periodic dimension*
   Self-explanatory.

*Cannot set non-zero z velocity for 2d simulation*
   Self-explanatory.

*Cannot set quaternion for atom that has none*
   Self-explanatory.

*Cannot set quaternion with xy components for 2d system*
   Self-explanatory.

*Cannot set respa hybrid and any of pair/inner/middle/outer*
   In the rRESPA integrator, you must compute pairwise potentials either
   all together (pair), with different cutoff regions (inner/middle/outer),
   or per hybrid sub-style (hybrid).  You cannot mix those.

*Cannot set respa middle without inner/outer*
   In the rRESPA integrator, you must define both a inner and outer
   setting in order to use a middle setting.

*Cannot set restart file size - MPI error: %s*
   This error was generated by MPI when reading/writing an MPI-IO restart
   file.

*Cannot set smd/contact/radius for this atom style*
   Self-explanatory.

*Cannot set smd/mass/density for this atom style*
   Self-explanatory.

*Cannot set temperature for fix rigid/nph*
   The temp keyword cannot be specified.

*Cannot set theta for atom that is not a line*
   Self-explanatory.

*Cannot set this attribute for this atom style*
   The attribute being set does not exist for the defined atom style.

*Cannot set variable z velocity for 2d simulation*
   Self-explanatory.

*Cannot skew triclinic box in z for 2d simulation*
   Self-explanatory.

*Cannot subtract groups using a dynamic group*
   This operation is not allowed.

*Cannot union groups using a dynamic group*
   This operation is not allowed.

*Cannot use -kokkos on without KOKKOS installed*
   Self-explanatory.

*Cannot use -reorder after -partition*
   Self-explanatory.  See doc page discussion of command-line switches.

*Cannot use Ewald with 2d simulation*
   The kspace style ewald cannot be used in 2d simulations.  You can use
   2d Ewald in a 3d simulation; see the kspace_modify command.

*Cannot use Ewald/disp solver on system with no charge, dipole, or LJ particles*
   No atoms in system have a non-zero charge or dipole, or are LJ
   particles.  Change charges/dipoles or change options of the kspace
   solver/pair style.

*Cannot use EwaldDisp with 2d simulation*
   This is a current restriction of this command.

*Cannot use Kokkos pair style with rRESPA inner/middle*
   Self-explanatory.

*Cannot use NEB unless atom map exists*
   Use the atom_modify command to create an atom map.

*Cannot use NEB with a single replica*
   Self-explanatory.

*Cannot use NEB with atom_modify sort enabled*
   This is current restriction for NEB implemented in LAMMPS.

*Cannot use PPPM with 2d simulation*
   The kspace style pppm cannot be used in 2d simulations.  You can use
   2d PPPM in a 3d simulation; see the kspace_modify command.

*Cannot use PPPMDisp with 2d simulation*
   The kspace style pppm/disp cannot be used in 2d simulations.  You can
   use 2d pppm/disp in a 3d simulation; see the kspace_modify command.

*Cannot use PRD with a changing box*
   The current box dimensions are not copied between replicas

*Cannot use PRD with a time-dependent fix defined*
   PRD alters the timestep in ways that will mess up these fixes.

*Cannot use PRD with a time-dependent region defined*
   PRD alters the timestep in ways that will mess up these regions.

*Cannot use PRD with atom_modify sort enabled*
   This is a current restriction of PRD.  You must turn off sorting,
   which is enabled by default, via the atom_modify command.

*Cannot use PRD with multi-processor replicas unless atom map exists*
   Use the atom_modify command to create an atom map.

*Cannot use TAD unless atom map exists for NEB*
   See atom_modify map command to set this.

*Cannot use TAD with a single replica for NEB*
   NEB requires multiple replicas.

*Cannot use TAD with atom_modify sort enabled for NEB*
   This is a current restriction of NEB.

*Cannot use a damped dynamics min style with fix box/relax*
   This is a current restriction in LAMMPS.  Use another minimizer
   style.

*Cannot use a damped dynamics min style with per-atom DOF*
   This is a current restriction in LAMMPS.  Use another minimizer
   style.

*Cannot use append/atoms in periodic dimension*
   The boundary style of the face where atoms are added can not be of
   type p (periodic).

*Cannot use atomfile-style variable unless atom map exists*
   Self-explanatory.  See the atom_modify command to create a map.

*Cannot use both com and bias with compute temp/chunk*
   Self-explanatory.

*Cannot use chosen neighbor list style with buck/coul/cut/kk*
   Self-explanatory.

*Cannot use chosen neighbor list style with buck/coul/long/kk*
   Self-explanatory.

*Cannot use chosen neighbor list style with buck/kk*
   That style is not supported by Kokkos.

*Cannot use chosen neighbor list style with coul/cut/kk*
   That style is not supported by Kokkos.

*Cannot use chosen neighbor list style with coul/debye/kk*
   Self-explanatory.

*Cannot use chosen neighbor list style with coul/dsf/kk*
   That style is not supported by Kokkos.

*Cannot use chosen neighbor list style with coul/wolf/kk*
   That style is not supported by Kokkos.

*Cannot use chosen neighbor list style with lj/charmm/coul/charmm/implicit/kk*
   Self-explanatory.

*Cannot use chosen neighbor list style with lj/charmm/coul/charmm/kk*
   Self-explanatory.

*Cannot use chosen neighbor list style with lj/charmm/coul/long/kk*
   Self-explanatory.

*Cannot use chosen neighbor list style with lj/class2/coul/cut/kk*
   Self-explanatory.

*Cannot use chosen neighbor list style with lj/class2/coul/long/kk*
   Self-explanatory.

*Cannot use chosen neighbor list style with lj/class2/kk*
   Self-explanatory.

*Cannot use chosen neighbor list style with lj/cut/coul/cut/kk*
   That style is not supported by Kokkos.

*Cannot use chosen neighbor list style with lj/cut/coul/debye/kk*
   Self-explanatory.

*Cannot use chosen neighbor list style with lj/cut/coul/long/kk*
   That style is not supported by Kokkos.

*Cannot use chosen neighbor list style with lj/cut/kk*
   That style is not supported by Kokkos.

*Cannot use chosen neighbor list style with lj/expand/kk*
   Self-explanatory.

*Cannot use chosen neighbor list style with lj/gromacs/coul/gromacs/kk*
   Self-explanatory.

*Cannot use chosen neighbor list style with lj/gromacs/kk*
   Self-explanatory.

*Cannot use chosen neighbor list style with lj/sdk/kk*
   That style is not supported by Kokkos.

*Cannot use chosen neighbor list style with pair eam/kk*
   That style is not supported by Kokkos.

*Cannot use chosen neighbor list style with pair eam/kk/alloy*
   Self-explanatory.

*Cannot use chosen neighbor list style with pair eam/kk/fs*
   Self-explanatory.

*Cannot use chosen neighbor list style with pair sw/kk*
   Self-explanatory.

*Cannot use chosen neighbor list style with tersoff/kk*
   Self-explanatory.

*Cannot use chosen neighbor list style with tersoff/zbl/kk*
   Self-explanatory.

*Cannot use compute chunk/atom bin z for 2d model*
   Self-explanatory.

*Cannot use compute cluster/atom unless atoms have IDs*
   Atom IDs are used to identify clusters.

*Cannot use create_atoms rotate unless single style*
   Self-explanatory.

*Cannot use create_bonds unless atoms have IDs*
   This command requires a mapping from global atom IDs to local atoms,
   but the atoms that have been defined have no IDs.

*Cannot use create_bonds with non-molecular system*
   Self-explanatory.

*Cannot use cwiggle in variable formula between runs*
   This is a function of elapsed time.

*Cannot use delete_atoms bond yes with atom_style template*
   This is because the bonds for that atom style are hardwired in the
   molecule template.

*Cannot use delete_atoms unless atoms have IDs*
   Your atoms do not have IDs, so the delete_atoms command cannot be
   used.

*Cannot use delete_bonds with non-molecular system*
   Your choice of atom style does not have bonds.

*Cannot use dump_modify fileper without % in dump file name*
   Self-explanatory.

*Cannot use dump_modify nfile without % in dump file name*
   Self-explanatory.

*Cannot use dynamic group with fix adapt atom*
   This is not yet supported.

*Cannot use fix TMD unless atom map exists*
   Using this fix requires the ability to lookup an atom index, which is
   provided by an atom map.  An atom map does not exist (by default) for
   non-molecular problems.  Using the atom_modify map command will force
   an atom map to be created.

*Cannot use fix bond/break with non-molecular systems*
   Only systems with bonds that can be changed can be used.  Atom_style
   template does not qualify.

*Cannot use fix bond/create with non-molecular systems*
   Only systems with bonds that can be changed can be used.  Atom_style
   template does not qualify.

*Cannot use fix bond/swap with non-molecular systems*
   Only systems with bonds that can be changed can be used.  Atom_style
   template does not qualify.

*Cannot use fix box/relax on a 2nd non-periodic dimension*
   When specifying an off-diagonal pressure component, the 2nd of the two
   dimensions must be periodic.  E.g. if the xy component is specified,
   then the y dimension must be periodic.

*Cannot use fix box/relax on a non-periodic dimension*
   When specifying a diagonal pressure component, the dimension must be
   periodic.

*Cannot use fix box/relax with both relaxation and scaling on a tilt factor*
   When specifying scaling on a tilt factor component, that component can not
   also be controlled by the barostat. E.g. if scalexy yes is specified and
   also keyword tri or xy, this is wrong.

*Cannot use fix box/relax with tilt factor scaling on a 2nd non-periodic dimension*
   When specifying scaling on a tilt factor component, the 2nd of the two
   dimensions must be periodic.  E.g. if the xy component is specified,
   then the y dimension must be periodic.

*Cannot use fix deform on a shrink-wrapped boundary*
   The x, y, z options cannot be applied to shrink-wrapped
   dimensions.

*Cannot use fix deform tilt on a shrink-wrapped 2nd dim*
   This is because the shrink-wrapping will change the value
   of the strain implied by the tilt factor.

*Cannot use fix deform trate on a box with zero tilt*
   The trate style alters the current strain.

*Cannot use fix deposit rigid and not molecule*
   Self-explanatory.

*Cannot use fix deposit rigid and shake*
   These two attributes are conflicting.

*Cannot use fix deposit shake and not molecule*
   Self-explanatory.

*Cannot use fix enforce2d with 3d simulation*
   Self-explanatory.

*Cannot use fix gcmc in a 2d simulation*
   Fix gcmc is set up to run in 3d only. No 2d simulations with fix gcmc
   are allowed.

*Cannot use fix gcmc shake and not molecule*
   Self-explanatory.

*Cannot use fix msst without per-type mass defined*
   Self-explanatory.

*Cannot use fix npt and fix deform on same component of stress tensor*
   This would be changing the same box dimension twice.

*Cannot use fix nvt/npt/nph on a 2nd non-periodic dimension*
   When specifying an off-diagonal pressure component, the 2nd of the two
   dimensions must be periodic.  E.g. if the xy component is specified,
   then the y dimension must be periodic.

*Cannot use fix nvt/npt/nph on a non-periodic dimension*
   When specifying a diagonal pressure component, the dimension must be
   periodic.

*Cannot use fix nvt/npt/nph with both xy dynamics and xy scaling*
   Self-explanatory.

*Cannot use fix nvt/npt/nph with both xz dynamics and xz scaling*
   Self-explanatory.

*Cannot use fix nvt/npt/nph with both yz dynamics and yz scaling*
   Self-explanatory.

*Cannot use fix nvt/npt/nph with xy scaling when y is non-periodic dimension*
   The 2nd dimension in the barostatted tilt factor must be periodic.

*Cannot use fix nvt/npt/nph with xz scaling when z is non-periodic dimension*
   The 2nd dimension in the barostatted tilt factor must be periodic.

*Cannot use fix nvt/npt/nph with yz scaling when z is non-periodic dimension*
   The 2nd dimension in the barostatted tilt factor must be periodic.

*Cannot use fix pour rigid and not molecule*
   Self-explanatory.

*Cannot use fix pour rigid and shake*
   These two attributes are conflicting.

*Cannot use fix pour shake and not molecule*
   Self-explanatory.

*Cannot use fix pour with triclinic box*
   This option is not yet supported.

*Cannot use fix press/berendsen and fix deform on same component of stress tensor*
   These commands both change the box size/shape, so you cannot use both
   together.

*Cannot use fix press/berendsen on a non-periodic dimension*
   Self-explanatory.

*Cannot use fix press/berendsen with triclinic box*
   Self-explanatory.

*Cannot use fix reax/bonds without pair_style reax*
   Self-explanatory.

*Cannot use fix rigid npt/nph and fix deform on same component of stress tensor*
   This would be changing the same box dimension twice.

*Cannot use fix rigid npt/nph on a non-periodic dimension*
   When specifying a diagonal pressure component, the dimension must be
   periodic.

*Cannot use fix rigid/small npt/nph on a non-periodic dimension*
   When specifying a diagonal pressure component, the dimension must be
   periodic.

*Cannot use fix shake with non-molecular system*
   Your choice of atom style does not have bonds.

*Cannot use fix ttm with 2d simulation*
   This is a current restriction of this fix due to the grid it creates.

*Cannot use fix ttm with triclinic box*
   This is a current restriction of this fix due to the grid it creates.

*Cannot use fix tune/kspace without a kspace style*
   Self-explanatory.

*Cannot use fix tune/kspace without a pair style*
   This fix (tune/kspace) can only be used when a pair style has been specified.

*Cannot use fix wall in periodic dimension*
   Self-explanatory.

*Cannot use fix wall zlo/zhi for a 2d simulation*
   Self-explanatory.

*Cannot use fix wall/reflect in periodic dimension*
   Self-explanatory.

*Cannot use fix wall/reflect zlo/zhi for a 2d simulation*
   Self-explanatory.

*Cannot use fix wall/srd in periodic dimension*
   Self-explanatory.

*Cannot use fix wall/srd more than once*
   Nor is their a need to since multiple walls can be specified
   in one command.

*Cannot use fix wall/srd without fix srd*
   Self-explanatory.

*Cannot use fix wall/srd zlo/zhi for a 2d simulation*
   Self-explanatory.

*Cannot use fix_deposit unless atoms have IDs*
   Self-explanatory.

*Cannot use fix_pour unless atoms have IDs*
   Self-explanatory.

*Cannot use include command within an if command*
   Self-explanatory.

*Cannot use lines with fix srd unless overlap is set*
   This is because line segments are connected to each other.

*Cannot use multiple fix wall commands with pair brownian*
   Self-explanatory.

*Cannot use multiple fix wall commands with pair lubricate*
   Self-explanatory.

*Cannot use multiple fix wall commands with pair lubricate/poly*
   Self-explanatory.

*Cannot use multiple fix wall commands with pair lubricateU*
   Self-explanatory.

*Cannot use neigh_modify exclude with GPU neighbor builds*
   This is a current limitation of the GPU implementation
   in LAMMPS.

*Cannot use neighbor bins - box size << cutoff*
   Too many neighbor bins will be created.  This typically happens when
   the simulation box is very small in some dimension, compared to the
   neighbor cutoff.  Use the "nsq" style instead of "bin" style.

*Cannot use newton pair with beck/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with born/coul/long/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with born/coul/wolf/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with born/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with buck/coul/cut/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with buck/coul/long/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with buck/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with colloid/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with coul/cut/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with coul/debye/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with coul/dsf/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with coul/long/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with dipole/cut/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with dipole/sf/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with dpd/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with dpd/tstat/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with eam/alloy/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with eam/fs/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with eam/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with gauss/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with gayberne/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with lj/charmm/coul/long/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with lj/class2/coul/long/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with lj/class2/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with lj/cubic/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with lj/cut/coul/cut/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with lj/cut/coul/debye/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with lj/cut/coul/dsf/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with lj/cut/coul/long/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with lj/cut/coul/msm/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with lj/cut/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with lj/expand/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with lj/gromacs/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with lj/sdk/coul/long/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with lj/sdk/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with lj96/cut/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with mie/cut/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with morse/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with resquared/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with soft/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with table/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with yukawa/colloid/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with yukawa/gpu pair style*
   Self-explanatory.

*Cannot use newton pair with zbl/gpu pair style*
   Self-explanatory.

*Cannot use non-zero forces in an energy minimization*
   Fix setforce cannot be used in this manner.  Use fix addforce
   instead.

*Cannot use non-periodic boundares with fix ttm*
   This fix requires a fully periodic simulation box.

*Cannot use non-periodic boundaries with Ewald*
   For kspace style ewald, all 3 dimensions must have periodic boundaries
   unless you use the kspace_modify command to define a 2d slab with a
   non-periodic z dimension.

*Cannot use non-periodic boundaries with EwaldDisp*
   For kspace style ewald/disp, all 3 dimensions must have periodic
   boundaries unless you use the kspace_modify command to define a 2d
   slab with a non-periodic z dimension.

*Cannot use non-periodic boundaries with PPPM*
   For kspace style pppm, all 3 dimensions must have periodic boundaries
   unless you use the kspace_modify command to define a 2d slab with a
   non-periodic z dimension.

*Cannot use non-periodic boundaries with PPPMDisp*
   For kspace style pppm/disp, all 3 dimensions must have periodic
   boundaries unless you use the kspace_modify command to define a 2d
   slab with a non-periodic z dimension.

*Cannot use order greater than 8 with pppm/gpu.*
   Self-explanatory.

*Cannot use package gpu neigh yes with triclinic box*
   This is a current restriction in LAMMPS.

*Cannot use pair tail corrections with 2d simulations*
   The correction factors are only currently defined for 3d systems.

*Cannot use processors part command without using partitions*
   See the command-line -partition switch.

*Cannot use ramp in variable formula between runs*
   This is because the ramp() function is time dependent.

*Cannot use read_data add before simulation box is defined*
   Self-explanatory.

*Cannot use read_data extra with add flag*
   Self-explanatory.

*Cannot use read_data offset without add flag*
   Self-explanatory.

*Cannot use read_data shift without add flag*
   Self-explanatory.

*Cannot use region INF or EDGE when box does not exist*
   Regions that extend to the box boundaries can only be used after the
   create_box command has been used.

*Cannot use set atom with no atom IDs defined*
   Atom IDs are not defined, so they cannot be used to identify an atom.

*Cannot use set mol with no molecule IDs defined*
   Self-explanatory.

*Cannot use swiggle in variable formula between runs*
   This is a function of elapsed time.

*Cannot use tris with fix srd unless overlap is set*
   This is because triangles are connected to each other.

*Cannot use variable energy with constant efield in fix efield*
   LAMMPS computes the energy itself when the E-field is constant.

*Cannot use variable energy with constant force in fix addforce*
   This is because for constant force, LAMMPS can compute the change
   in energy directly.

*Cannot use variable every setting for dump dcd*
   The format of DCD dump files requires snapshots be output
   at a constant frequency.

*Cannot use variable every setting for dump xtc*
   The format of this file requires snapshots at regular intervals.

*Cannot use vdisplace in variable formula between runs*
   This is a function of elapsed time.

*Cannot use velocity bias command without temp keyword*
   Self-explanatory.

*Cannot use velocity create loop all unless atoms have IDs*
   Atoms in the simulation to do not have IDs, so this style
   of velocity creation cannot be performed.

*Cannot use wall in periodic dimension*
   Self-explanatory.

*Cannot use write_restart fileper without % in restart file name*
   Self-explanatory.

*Cannot use write_restart nfile without % in restart file name*
   Self-explanatory.

*Cannot wiggle and shear fix wall/gran*
   Cannot specify both options at the same time.

*Cannot write to restart file - MPI error: %s*
   This error was generated by MPI when reading/writing an MPI-IO restart
   file.

*Cannot yet use KSpace solver with grid with comm style tiled*
   This is current restriction in LAMMPS.

*Cannot yet use comm_style tiled with multi-mode comm*
   Self-explanatory.

*Cannot yet use comm_style tiled with triclinic box*
   Self-explanatory.

*Cannot yet use compute tally with Kokkos*
   This feature is not yet supported.

*Cannot yet use fix bond/break with this improper style*
   This is a current restriction in LAMMPS.

*Cannot yet use fix bond/create with this improper style*
   This is a current restriction in LAMMPS.

*Cannot yet use minimize with Kokkos*
   This feature is not yet supported.

*Cannot yet use pair hybrid with Kokkos*
   This feature is not yet supported.

*Cannot zero Langevin force of 0 atoms*
   The group has zero atoms, so you cannot request its force
   be zeroed.

*Cannot zero gld force for zero atoms*
   There are no atoms currently in the group.

*Cannot zero momentum of no atoms*
   Self-explanatory.

*Change_box command before simulation box is defined*
   Self-explanatory.

*Change_box volume used incorrectly*
   The "dim volume" option must be used immediately following one or two
   settings for "dim1 ..." (and optionally "dim2 ...") and must be for a
   different dimension, i.e. dim != dim1 and dim != dim2.

*Chunk/atom compute does not exist for compute angmom/chunk*
   Self-explanatory.

*Chunk/atom compute does not exist for compute com/chunk*
   Self-explanatory.

*Chunk/atom compute does not exist for compute gyration/chunk*
   Self-explanatory.

*Chunk/atom compute does not exist for compute inertia/chunk*
   Self-explanatory.

*Chunk/atom compute does not exist for compute msd/chunk*
   Self-explanatory.

*Chunk/atom compute does not exist for compute omega/chunk*
   Self-explanatory.

*Chunk/atom compute does not exist for compute property/chunk*
   Self-explanatory.

*Chunk/atom compute does not exist for compute temp/chunk*
   Self-explanatory.

*Chunk/atom compute does not exist for compute torque/chunk*
   Self-explanatory.

*Chunk/atom compute does not exist for compute vcm/chunk*
   Self-explanatory.

*Chunk/atom compute does not exist for fix ave/chunk*
   Self-explanatory.

*Comm tiled invalid index in box drop brick*
   Internal error check in comm_style tiled which should not occur.
   Contact the developers.

*Comm tiled mis-match in box drop brick*
   Internal error check in comm_style tiled which should not occur.
   Contact the developers.

*Comm_modify group != atom_modify first group*
   Self-explanatory.

*Communication cutoff for comm_style tiled cannot exceed periodic box length*
   Self-explanatory.

*Communication cutoff too small for SNAP micro load balancing*
   This can happen if you change the neighbor skin after your pair_style
   command or if your box dimensions grow during a run. You can set the
   cutoff explicitly via the comm_modify cutoff command.

*Compute %s does not allow use of dynamic group*
   Dynamic groups have not yet been enabled for this compute.

*Compute ID for compute chunk /atom does not exist*
   Self-explanatory.

*Compute ID for compute chunk/atom does not exist*
   Self-explanatory.

*Compute gyration ID does not exist for compute gyration/shape*
   Self-explanatory. Provide a valid compute ID.

*Compute gyration/shape compute ID does not point to a gyration compute*
   Self-explanatory. Provide and ID of a compute gyration command.

*Compute ID for compute reduce does not exist*
   Self-explanatory.

*Compute ID for compute slice does not exist*
   Self-explanatory.

*Compute ID for fix ave/atom does not exist*
   Self-explanatory.

*Compute ID for fix ave/chunk does not exist*
   Self-explanatory.

*Compute ID for fix ave/correlate does not exist*
   Self-explanatory.

*Compute ID for fix ave/histo does not exist*
   Self-explanatory.

*Compute ID for fix ave/time does not exist*
   Self-explanatory.

*Compute ID for fix numdiff does not exist*
   Self-explanatory.

*Compute ID for fix store/state does not exist*
   Self-explanatory.

*Compute ID for fix vector does not exist*
   Self-explanatory.

*Compute ID must be alphanumeric or underscore characters*
   Self-explanatory.

*Compute angle/local used when angles are not allowed*
   The atom style does not support angles.

*Compute angmom/chunk does not use chunk/atom compute*
   The style of the specified compute is not chunk/atom.

*Compute body/local requires atom style body*
   Self-explanatory.

*Compute bond/local used when bonds are not allowed*
   The atom style does not support bonds.

*Compute centro/atom requires a pair style be defined*
   This is because the computation of the centro-symmetry values
   uses a pairwise neighbor list.

*Compute chunk/atom bin/cylinder radius is too large for periodic box*
   Radius cannot be bigger than 1/2 of a non-axis  periodic dimension.

*Compute chunk/atom bin/sphere radius is too large for periodic box*
   Radius cannot be bigger than 1/2 of any periodic dimension.

*Compute chunk/atom compute array is accessed out-of-range*
   The index for the array is out of bounds.

*Compute chunk/atom compute does not calculate a per-atom array*
   Self-explanatory.

*Compute chunk/atom compute does not calculate a per-atom vector*
   Self-explanatory.

*Compute chunk/atom compute does not calculate per-atom values*
   Self-explanatory.

*Compute chunk/atom cylinder axis must be z for 2d*
   Self-explanatory.

*Compute chunk/atom fix array is accessed out-of-range*
   the index for the array is out of bounds.

*Compute chunk/atom fix does not calculate a per-atom array*
   Self-explanatory.

*Compute chunk/atom fix does not calculate a per-atom vector*
   Self-explanatory.

*Compute chunk/atom fix does not calculate per-atom values*
   Self-explanatory.

*Compute chunk/atom for triclinic boxes requires units reduced*
   Self-explanatory.

*Compute chunk/atom ids once but nchunk is not once*
   You cannot assign chunks IDs to atom permanently if the number of
   chunks may change.

*Compute chunk/atom molecule for non-molecular system*
   Self-explanatory.

*Compute chunk/atom sphere z origin must be 0.0 for 2d*
   Self-explanatory.

*Compute chunk/atom stores no IDs for compute property/chunk*
   It will only store IDs if its compress option is enabled.

*Compute chunk/atom stores no coord1 for compute property/chunk*
   Only certain binning options for compute chunk/atom store coordinates.

*Compute chunk/atom stores no coord2 for compute property/chunk*
   Only certain binning options for compute chunk/atom store coordinates.

*Compute chunk/atom stores no coord3 for compute property/chunk*
   Only certain binning options for compute chunk/atom store coordinates.

*Compute chunk/atom variable is not atom-style variable*
   Self-explanatory.

*Compute chunk/atom without bins cannot use discard mixed*
   That discard option only applies to the binning styles.

*Compute cluster/atom cutoff is longer than pairwise cutoff*
   Cannot identify clusters beyond cutoff.

*Compute cluster/atom requires a pair style be defined*
   This is so that the pair style defines a cutoff distance which
   is used to find clusters.

*Compute cna/atom cutoff is longer than pairwise cutoff*
   Self-explanatory.

*Compute cna/atom requires a pair style be defined*
   Self-explanatory.

*Compute com/chunk does not use chunk/atom compute*
   The style of the specified compute is not chunk/atom.

*Compute contact/atom requires a pair style be defined*
   Self-explanatory.

*Compute contact/atom requires atom style sphere*
   Self-explanatory.

*Compute coord/atom cutoff is longer than pairwise cutoff*
   Cannot compute coordination at distances longer than the pair cutoff,
   since those atoms are not in the neighbor list.

*Compute coord/atom requires a pair style be defined*
   Self-explanatory.

*Compute damage/atom requires peridynamic potential*
   Damage is a Peridynamic-specific metric.  It requires you
   to be running a Peridynamics simulation.

*Compute dihedral/local used when dihedrals are not allowed*
   The atom style does not support dihedrals.

*Compute dilatation/atom cannot be used with this pair style*
   Self-explanatory.

*Compute dilatation/atom requires Peridynamic pair style*
   Self-explanatory.

*Compute does not allow an extra compute or fix to be reset*
   This is an internal LAMMPS error.  Please report it to the
   developers.

*Compute erotate/asphere requires atom style ellipsoid or line or tri*
   Self-explanatory.

*Compute erotate/asphere requires extended particles*
   This compute cannot be used with point particles.

*Compute erotate/rigid with non-rigid fix-ID*
   Self-explanatory.

*Compute erotate/sphere requires atom style sphere*
   Self-explanatory.

*Compute erotate/sphere/atom requires atom style sphere*
   Self-explanatory.

*Compute event/displace has invalid fix event assigned*
   This is an internal LAMMPS error.  Please report it to the
   developers.

*Compute group/group group ID does not exist*
   Self-explanatory.

*Compute gyration/chunk does not use chunk/atom compute*
   The style of the specified compute is not chunk/atom.

*Compute heat/flux compute ID does not compute ke/atom*
   Self-explanatory.

*Compute heat/flux compute ID does not compute pe/atom*
   Self-explanatory.

*Compute heat/flux compute ID does not compute stress/atom*
   Self-explanatory.

*Compute hexorder/atom cutoff is longer than pairwise cutoff*
   Cannot compute order parameter beyond cutoff.

*Compute hexorder/atom requires a pair style be defined*
   Self-explanatory.

*Compute improper/local used when impropers are not allowed*
   The atom style does not support impropers.

*Compute inertia/chunk does not use chunk/atom compute*
   The style of the specified compute is not chunk/atom.

*Compute ke/rigid with non-rigid fix-ID*
   Self-explanatory.

*Compute msd/chunk does not use chunk/atom compute*
   The style of the specified compute is not chunk/atom.

*Compute msd/chunk nchunk is not static*
   This is required because the MSD cannot be computed consistently if
   the number of chunks is changing.  Compute chunk/atom allows setting
   nchunk to be static.

*Compute nve/asphere requires atom style ellipsoid*
   Self-explanatory.

*Compute nvt/nph/npt asphere requires atom style ellipsoid*
   Self-explanatory.

*Compute nvt/nph/npt body requires atom style body*
   Self-explanatory.

*Compute omega/chunk does not use chunk/atom compute*
   The style of the specified compute is not chunk/atom.

*Compute orientorder/atom cutoff is longer than pairwise cutoff*
   Cannot compute order parameter beyond cutoff.

*Compute orientorder/atom requires a pair style be defined*
   Self-explanatory.

*Compute pair must use group all*
   Pair styles accumulate energy on all atoms.

*Compute pe must use group all*
   Energies computed by potentials (pair, bond, etc) are computed on all
   atoms.

*Compute plasticity/atom cannot be used with this pair style*
   Self-explanatory.

*Compute plasticity/atom requires Peridynamic pair style*
   Self-explanatory.

*Compute pressure must use group all*
   Virial contributions computed by potentials (pair, bond, etc) are
   computed on all atoms.

*Compute pressure requires temperature ID to include kinetic energy*
   The keflag cannot be used unless a temperature compute is provided.

*Compute pressure temperature ID does not compute temperature*
   The compute ID assigned to a pressure computation must compute
   temperature.

*Compute property/atom floating point vector does not exist*
   The command is accessing a vector added by the fix property/atom
   command, that does not exist.

*Compute property/atom for atom property that isn't allocated*
   Self-explanatory.

*Compute property/atom integer vector does not exist*
   The command is accessing a vector added by the fix property/atom
   command, that does not exist.

*Compute property/chunk does not use chunk/atom compute*
   The style of the specified compute is not chunk/atom.

*Compute property/local cannot use these inputs together*
   Only inputs that generate the same number of datums can be used
   together.  E.g. bond and angle quantities cannot be mixed.

*Compute property/local does not (yet) work with atom_style template*
   Self-explanatory.

*Compute property/local for property that isn't allocated*
   Self-explanatory.

*Compute rdf requires a pair style be defined*
   Self-explanatory.

*Compute reduce compute array is accessed out-of-range*
   An index for the array is out of bounds.

*Compute reduce compute calculates global values*
   A compute that calculates peratom or local values is required.

*Compute reduce compute does not calculate a local array*
   Self-explanatory.

*Compute reduce compute does not calculate a local vector*
   Self-explanatory.

*Compute reduce compute does not calculate a per-atom array*
   Self-explanatory.

*Compute reduce compute does not calculate a per-atom vector*
   Self-explanatory.

*Compute reduce fix array is accessed out-of-range*
   An index for the array is out of bounds.

*Compute reduce fix calculates global values*
   A fix that calculates peratom or local values is required.

*Compute reduce fix does not calculate a local array*
   Self-explanatory.

*Compute reduce fix does not calculate a local vector*
   Self-explanatory.

*Compute reduce fix does not calculate a per-atom array*
   Self-explanatory.

*Compute reduce fix does not calculate a per-atom vector*
   Self-explanatory.

*Compute reduce replace requires min or max mode*
   Self-explanatory.

*Compute reduce variable is not atom-style variable*
   Self-explanatory.

*Compute slice compute array is accessed out-of-range*
   An index for the array is out of bounds.

*Compute slice compute does not calculate a global array*
   Self-explanatory.

*Compute slice compute does not calculate a global vector*
   Self-explanatory.

*Compute slice compute does not calculate global vector or array*
   Self-explanatory.

*Compute slice compute vector is accessed out-of-range*
   The index for the vector is out of bounds.

*Compute slice fix array is accessed out-of-range*
   An index for the array is out of bounds.

*Compute slice fix does not calculate a global array*
   Self-explanatory.

*Compute slice fix does not calculate a global vector*
   Self-explanatory.

*Compute slice fix does not calculate global vector or array*
   Self-explanatory.

*Compute slice fix vector is accessed out-of-range*
   The index for the vector is out of bounds.

*Compute sna/atom cutoff is longer than pairwise cutoff*
   Self-explanatory.

*Compute sna/atom requires a pair style be defined*
   Self-explanatory.

*Compute snad/atom cutoff is longer than pairwise cutoff*
   Self-explanatory.

*Compute snad/atom requires a pair style be defined*
   Self-explanatory.

*Compute snav/atom cutoff is longer than pairwise cutoff*
   Self-explanatory.

*Compute snav/atom requires a pair style be defined*
   Self-explanatory.

*Compute stress/atom temperature ID does not compute temperature*
   The specified compute must compute temperature.

*Compute temp/asphere requires atom style ellipsoid*
   Self-explanatory.

*Compute temp/asphere requires extended particles*
   This compute cannot be used with point particles.

*Compute temp/body requires atom style body*
   Self-explanatory.

*Compute temp/body requires bodies*
   This compute can only be applied to body particles.

*Compute temp/chunk does not use chunk/atom compute*
   The style of the specified compute is not chunk/atom.

*Compute temp/cs requires ghost atoms store velocity*
   Use the comm_modify vel yes command to enable this.

*Compute temp/cs used when bonds are not allowed*
   This compute only works on pairs of bonded particles.

*Compute temp/partial cannot use vz for 2d systemx*
   Self-explanatory.

*Compute temp/profile cannot bin z for 2d systems*
   Self-explanatory.

*Compute temp/profile cannot use vz for 2d systemx*
   Self-explanatory.

*Compute temp/sphere requires atom style sphere*
   Self-explanatory.

*Compute ti kspace style does not exist*
   Self-explanatory.

*Compute ti pair style does not exist*
   Self-explanatory.

*Compute ti tail when pair style does not compute tail corrections*
   Self-explanatory.

*Compute torque/chunk does not use chunk/atom compute*
   The style of the specified compute is not chunk/atom.

*Compute used in dump between runs is not current*
   The compute was not invoked on the current timestep, therefore it
   cannot be used in a dump between runs.

*Compute used in variable between runs is not current*
   Computes cannot be invoked by a variable in between runs.  Thus they
   must have been evaluated on the last timestep of the previous run in
   order for their value(s) to be accessed.  See the doc page for the
   variable command for more info.

*Compute used in variable thermo keyword between runs is not current*
   Some thermo keywords rely on a compute to calculate their value(s).
   Computes cannot be invoked by a variable in between runs.  Thus they
   must have been evaluated on the last timestep of the previous run in
   order for their value(s) to be accessed.  See the doc page for the
   variable command for more info.

*Compute vcm/chunk does not use chunk/atom compute*
   The style of the specified compute is not chunk/atom.

*Computed temperature for fix temp/berendsen cannot be 0.0*
   Self-explanatory.

*Computed temperature for fix temp/rescale cannot be 0.0*
   Cannot rescale the temperature to a new value if the current
   temperature is 0.0.

*Core/shell partner atom not found*
   Could not find one of the atoms in the bond pair.

*Core/shell partners were not all found*
   Could not find or more atoms in the bond pairs.

*Could not adjust g_ewald_6*
   The Newton-Raphson solver failed to converge to a good value for
   g_ewald.  This error should not occur for typical problems.  Please
   send an email to the developers.

*Could not compute g_ewald*
   The Newton-Raphson solver failed to converge to a good value for
   g_ewald.  This error should not occur for typical problems.  Please
   send an email to the developers.

*Could not compute grid size*
   The code is unable to compute a grid size consistent with the desired
   accuracy.  This error should not occur for typical problems.  Please
   send an email to the developers.

*Could not compute grid size for Coulomb interaction*
   The code is unable to compute a grid size consistent with the desired
   accuracy.  This error should not occur for typical problems.  Please
   send an email to the developers.

*Could not compute grid size for Dispersion*
   The code is unable to compute a grid size consistent with the desired
   accuracy.  This error should not occur for typical problems.  Please
   send an email to the developers.

*Could not create 3d FFT plan*
   The FFT setup for the PPPM solver failed, typically due
   to lack of memory.  This is an unusual error.  Check the
   size of the FFT grid you are requesting.

*Could not create 3d grid of processors*
   The specified constraints did not allow a Px by Py by Pz grid to be
   created where Px \* Py \* Pz = P = total number of processors.

*Could not create 3d remap plan*
   The FFT setup in pppm failed.

*Could not create Python function arguments*
   This is an internal Python error, possibly because the number
   of inputs to the function is too large.

*Could not create numa grid of processors*
   The specified constraints did not allow this style of grid to be
   created.  Usually this is because the total processor count is not a
   multiple of the cores/node or the user specified processor count is >
   1 in one of the dimensions.

*Could not create twolevel 3d grid of processors*
   The specified constraints did not allow this style of grid to be
   created.

*Could not evaluate Python function input variable*
   Self-explanatory.

*Could not find Python function*
   The provided Python code was run successfully, but it not
   define a callable function with the required name.

*Could not find atom_modify first group ID*
   Self-explanatory.

*Could not find change_box group ID*
   Group ID used in the change_box command does not exist.

*Could not find compute ID for PRD*
   Self-explanatory.

*Could not find compute ID for TAD*
   Self-explanatory.

*Could not find compute ID for temperature bias*
   Self-explanatory.

*Could not find compute ID to delete*
   Self-explanatory.

*Could not find compute displace/atom fix ID*
   Self-explanatory.

*Could not find compute event/displace fix ID*
   Self-explanatory.

*Could not find compute group ID*
   Self-explanatory.

*Could not find compute heat/flux compute ID*
   Self-explanatory.

*Could not find compute msd fix ID*
   Self-explanatory.

*Could not find compute msd/chunk fix ID*
   The compute creates an internal fix, which has been deleted.

*Could not find compute pressure temperature ID*
   The compute ID for calculating temperature does not exist.

*Could not find compute stress/atom temperature ID*
   Self-explanatory.

*Could not find compute vacf fix ID*
   Self-explanatory.

*Could not find compute/voronoi surface group ID*
   Self-explanatory.

*Could not find compute_modify ID*
   Self-explanatory.

*Could not find custom per-atom property ID*
   Self-explanatory.

*Could not find delete_atoms group ID*
   Group ID used in the delete_atoms command does not exist.

*Could not find delete_atoms region ID*
   Region ID used in the delete_atoms command does not exist.

*Could not find displace_atoms group ID*
   Group ID used in the displace_atoms command does not exist.

*Could not find dump custom compute ID*
   Self-explanatory.

*Could not find dump custom fix ID*
   Self-explanatory.

*Could not find dump custom variable name*
   Self-explanatory.

*Could not find dump group ID*
   A group ID used in the dump command does not exist.

*Could not find dump local compute ID*
   Self-explanatory.

*Could not find dump local fix ID*
   Self-explanatory.

*Could not find dump modify compute ID*
   Self-explanatory.

*Could not find dump modify custom atom floating point property ID*
   Self-explanatory.

*Could not find dump modify custom atom integer property ID*
   Self-explanatory.

*Could not find dump modify fix ID*
   Self-explanatory.

*Could not find dump modify variable name*
   Self-explanatory.

*Could not find fix ID to delete*
   Self-explanatory.

*Could not find fix adapt storage fix ID*
   This should not happen unless you explicitly deleted
   a secondary fix that fix adapt created internally.

*Could not find fix gcmc exclusion group ID*
   Self-explanatory.

*Could not find fix gcmc rotation group ID*
   Self-explanatory.

*Could not find fix group ID*
   A group ID used in the fix command does not exist.

*Could not find fix msst compute ID*
   Self-explanatory.

*Could not find fix poems group ID*
   A group ID used in the fix poems command does not exist.

*Could not find fix recenter group ID*
   A group ID used in the fix recenter command does not exist.

*Could not find fix rigid group ID*
   A group ID used in the fix rigid command does not exist.

*Could not find fix srd group ID*
   Self-explanatory.

*Could not find fix_modify ID*
   A fix ID used in the fix_modify command does not exist.

*Could not find fix_modify pressure ID*
   The compute ID for computing pressure does not exist.

*Could not find fix_modify temperature ID*
   The compute ID for computing temperature does not exist.

*Could not find group clear group ID*
   Self-explanatory.

*Could not find group delete group ID*
   Self-explanatory.

*Could not find pair fix ID*
   A fix is created internally by the pair style to store shear
   history information.  You cannot delete it.

*Could not find set group ID*
   Group ID specified in set command does not exist.

*Could not find specified fix gcmc group ID*
   Self-explanatory.

*Could not find thermo compute ID*
   Compute ID specified in thermo_style command does not exist.

*Could not find thermo custom compute ID*
   The compute ID needed by thermo style custom to compute a requested
   quantity does not exist.

*Could not find thermo custom fix ID*
   The fix ID needed by thermo style custom to compute a requested
   quantity does not exist.

*Could not find thermo custom variable name*
   Self-explanatory.

*Could not find thermo fix ID*
   Fix ID specified in thermo_style command does not exist.

*Could not find thermo variable name*
   Self-explanatory.

*Could not find thermo_modify pressure ID*
   The compute ID needed by thermo style custom to compute pressure does
   not exist.

*Could not find thermo_modify temperature ID*
   The compute ID needed by thermo style custom to compute temperature does
   not exist.

*Could not find undump ID*
   A dump ID used in the undump command does not exist.

*Could not find velocity group ID*
   A group ID used in the velocity command does not exist.

*Could not find velocity temperature ID*
   The compute ID needed by the velocity command to compute temperature
   does not exist.

*Could not find/initialize a specified accelerator device*
   Could not initialize at least one of the devices specified for the gpu
   package

*Could not grab element entry from EIM potential file*
   Self-explanatory

*Could not grab global entry from EIM potential file*
   Self-explanatory.

*Could not grab pair entry from EIM potential file*
   Self-explanatory.

*Could not initialize embedded Python*
   The main module in Python was not accessible.

*Could not open Python file*
   The specified file of Python code cannot be opened.  Check that the
   path and name are correct.

*Could not process Python file*
   The Python code in the specified file was not run successfully by
   Python, probably due to errors in the Python code.

*Could not process Python string*
   The Python code in the here string was not run successfully by Python,
   probably due to errors in the Python code.

*Coulomb PPPMDisp order has been reduced below minorder*
   The default minimum order is 2.  This can be reset by the
   kspace_modify minorder command.

*Coulombic cutoff not supported in pair_style buck/long/coul/coul*
   Must use long-range Coulombic interactions.

*Coulombic cutoff not supported in pair_style lj/long/coul/long*
   Must use long-range Coulombic interactions.

*Coulombic cutoff not supported in pair_style lj/long/tip4p/long*
   Must use long-range Coulombic interactions.

*Coulombic cutoffs of pair hybrid sub-styles do not match*
   If using a Kspace solver, all Coulombic cutoffs of long pair styles must
   be the same.

*Coulombic cut not supported in pair_style lj/long/dipole/long*
   Must use long-range Coulombic interactions.

*Cound not find dump_modify ID*
   Self-explanatory.

*Create_atoms command before simulation box is defined*
   The create_atoms command cannot be used before a read_data,
   read_restart, or create_box command.

*Create_atoms molecule has atom IDs, but system does not*
   The atom_style id command can be used to force atom IDs to be stored.

*Create_atoms molecule must have atom types*
   The defined molecule does not specify atom types.

*Create_atoms molecule must have coordinates*
   The defined molecule does not specify coordinates.

*Create_atoms region ID does not exist*
   A region ID used in the create_atoms command does not exist.

*Create_bonds command before simulation box is defined*
   Self-explanatory.

*Create_bonds command requires no kspace_style be defined*
   This is so that atom pairs that are already bonded to not appear
   in the neighbor list.

*Create_bonds command requires special_bonds 1-2 weights be 0.0*
   This is so that atom pairs that are already bonded to not appear in
   the neighbor list.

*Create_bonds max distance > neighbor cutoff*
   Can only create bonds for atom pairs that will be in neighbor list.

*Create_bonds requires a pair style be defined*
   Self-explanatory.

*Create_box region ID does not exist*
   Self-explanatory.

*Create_box region does not support a bounding box*
   Not all regions represent bounded volumes.  You cannot use
   such a region with the create_box command.

*Custom floating point vector for fix store/state does not exist*
   The command is accessing a vector added by the fix property/atom
   command, that does not exist.

*Custom integer vector for fix store/state does not exist*
   The command is accessing a vector added by the fix property/atom
   command, that does not exist.

*Custom per-atom property ID is not floating point*
   Self-explanatory.

*Custom per-atom property ID is not integer*
   Self-explanatory.

*Cut-offs missing in pair_style lj/long/dipole/long*
   Self-explanatory.

*Cutoffs missing in pair_style buck/long/coul/long*
   Self-explanatory.

*Cutoffs missing in pair_style lj/long/coul/long*
   Self-explanatory.

*Cyclic loop in joint connections*
   Fix poems cannot (yet) work with coupled bodies whose joints connect
   the bodies in a ring (or cycle).

*Degenerate lattice primitive vectors*
   Invalid set of 3 lattice vectors for lattice command.

*Delete region ID does not exist*
   Self-explanatory.

*Delete_atoms command before simulation box is defined*
   The delete_atoms command cannot be used before a read_data,
   read_restart, or create_box command.

*Delete_atoms cutoff > max neighbor cutoff*
   Can only delete atoms in atom pairs that will be in neighbor list.

*Delete_atoms mol yes requires atom attribute molecule*
   Cannot use this option with a non-molecular system.

*Delete_atoms requires a pair style be defined*
   This is because atom deletion within a cutoff uses a pairwise
   neighbor list.

*Delete_bonds command before simulation box is defined*
   The delete_bonds command cannot be used before a read_data,
   read_restart, or create_box command.

*Delete_bonds command with no atoms existing*
   No atoms are yet defined so the delete_bonds command cannot be used.

*Deposition region extends outside simulation box*
   Self-explanatory.

*Did not assign all atoms correctly*
   Atoms read in from a data file were not assigned correctly to
   processors.  This is likely due to some atom coordinates being
   outside a non-periodic simulation box.

*Did not assign all restart atoms correctly*
   Atoms read in from the restart file were not assigned correctly to
   processors.  This is likely due to some atom coordinates being outside
   a non-periodic simulation box.  Normally this should not happen.  You
   may wish to use the "remap" option on the read_restart command to see
   if this helps.

*Did not find all elements in MEAM library file*
   The requested elements were not found in the MEAM file.

*Did not find fix shake partner info*
   Could not find bond partners implied by fix shake command.  This error
   can be triggered if the delete_bonds command was used before fix
   shake, and it removed bonds without resetting the 1-2, 1-3, 1-4
   weighting list via the special keyword.

*Did not find keyword in table file*
   Keyword used in pair_coeff command was not found in table file.

*Did not set pressure for fix rigid/nph*
   The press keyword must be specified.

*Did not set temp for fix rigid/nvt/small*
   Self-explanatory.

*Did not set temp or press for fix rigid/npt/small*
   Self-explanatory.

*Did not set temperature for fix rigid/nvt*
   The temp keyword must be specified.

*Did not set temperature or pressure for fix rigid/npt*
   The temp and press keywords must be specified.

*Dihedral atom missing in delete_bonds*
   The delete_bonds command cannot find one or more atoms in a particular
   dihedral on a particular processor.  The pairwise cutoff is too short
   or the atoms are too far apart to make a valid dihedral.

*Dihedral atom missing in set command*
   The set command cannot find one or more atoms in a particular dihedral
   on a particular processor.  The pairwise cutoff is too short or the
   atoms are too far apart to make a valid dihedral.

*Dihedral atoms %d %d %d %d missing on proc %d at step %ld*
   One or more of 4 atoms needed to compute a particular dihedral are
   missing on this processor.  Typically this is because the pairwise
   cutoff is set too short or the dihedral has blown apart and an atom is
   too far away.

*Dihedral atoms missing on proc %d at step %ld*
   One or more of 4 atoms needed to compute a particular dihedral are
   missing on this processor.  Typically this is because the pairwise
   cutoff is set too short or the dihedral has blown apart and an atom is
   too far away.

*Dihedral charmm is incompatible with Pair style*
   Dihedral style charmm must be used with a pair style charmm
   in order for the 1-4 epsilon/sigma parameters to be defined.

*Dihedral coeff for hybrid has invalid style*
   Dihedral style hybrid uses another dihedral style as one of its
   coefficients.  The dihedral style used in the dihedral_coeff command
   or read from a restart file is not recognized.

*Dihedral coeffs are not set*
   No dihedral coefficients have been assigned in the data file or via
   the dihedral_coeff command.

*Dihedral style hybrid cannot have hybrid as an argument*
   Self-explanatory.

*Dihedral style hybrid cannot have none as an argument*
   Self-explanatory.

*Dihedral style hybrid cannot use same dihedral style twice*
   Self-explanatory.

*Dihedral/improper extent > half of periodic box length*
   This error was detected by the neigh_modify check yes setting.  It is
   an error because the dihedral atoms are so far apart it is ambiguous
   how it should be defined.

*Dihedral_coeff command before dihedral_style is defined*
   Coefficients cannot be set in the data file or via the dihedral_coeff
   command until an dihedral_style has been assigned.

*Dihedral_coeff command before simulation box is defined*
   The dihedral_coeff command cannot be used before a read_data,
   read_restart, or create_box command.

*Dihedral_coeff command when no dihedrals allowed*
   The chosen atom style does not allow for dihedrals to be defined.

*Dihedral_style command when no dihedrals allowed*
   The chosen atom style does not allow for dihedrals to be defined.

*Dihedrals assigned incorrectly*
   Dihedrals read in from the data file were not assigned correctly to
   atoms.  This means there is something invalid about the topology
   definitions.

*Dihedrals defined but no dihedral types*
   The data file header lists dihedrals but no dihedral types.

*Dimension command after simulation box is defined*
   The dimension command cannot be used after a read_data,
   read_restart, or create_box command.

*Dispersion PPPMDisp order has been reduced below minorder*
   The default minimum order is 2.  This can be reset by the
   kspace_modify minorder command.

*Displace_atoms command before simulation box is defined*
   The displace_atoms command cannot be used before a read_data,
   read_restart, or create_box command.

*Distance must be > 0 for compute event/displace*
   Self-explanatory.

*Divide by 0 in influence function*
   This should not normally occur.  It is likely a problem with your
   model.

*Divide by 0 in influence function of pair peri/lps*
   This should not normally occur.  It is likely a problem with your
   model.

*Divide by 0 in variable formula*
   Self-explanatory.

*Domain too large for neighbor bins*
   The domain has become extremely large so that neighbor bins cannot be
   used.  Most likely, one or more atoms have been blown out of the
   simulation box to a great distance.

*Double precision is not supported on this accelerator*
   Self-explanatory

*Dump atom/gz only writes compressed files*
   The dump atom/gz output file name must have a .gz suffix.

*Dump cfg arguments can not mix xs\|ys\|zs with xsu\|ysu\|zsu*
   Self-explanatory.

*Dump cfg arguments must start with 'mass type xs ys zs' or 'mass type xsu ysu zsu'*
   This is a requirement of the CFG output format.  See the dump cfg doc
   page for more details.

*Dump cfg requires one snapshot per file*
   Use the wildcard "\*" character in the filename.

*Dump cfg/gz only writes compressed files*
   The dump cfg/gz output file name must have a .gz suffix.

*Dump custom and fix not computed at compatible times*
   The fix must produce per-atom quantities on timesteps that dump custom
   needs them.

*Dump custom compute does not calculate per-atom array*
   Self-explanatory.

*Dump custom compute does not calculate per-atom vector*
   Self-explanatory.

*Dump custom compute does not compute per-atom info*
   Self-explanatory.

*Dump custom compute vector is accessed out-of-range*
   Self-explanatory.

*Dump custom fix does not compute per-atom array*
   Self-explanatory.

*Dump custom fix does not compute per-atom info*
   Self-explanatory.

*Dump custom fix does not compute per-atom vector*
   Self-explanatory.

*Dump custom fix vector is accessed out-of-range*
   Self-explanatory.

*Dump custom variable is not atom-style variable*
   Only atom-style variables generate per-atom quantities, needed for
   dump output.

*Dump custom/gz only writes compressed files*
   The dump custom/gz output file name must have a .gz suffix.

*Dump dcd of non-matching # of atoms*
   Every snapshot written by dump dcd must contain the same # of atoms.

*Dump dcd requires sorting by atom ID*
   Use the dump_modify sort command to enable this.

*Dump every variable returned a bad timestep*
   The variable must return a timestep greater than the current timestep.

*Dump file MPI-IO output not allowed with % in filename*
   This is because a % signifies one file per processor and MPI-IO
   creates one large file for all processors.

*Dump file does not contain requested snapshot*
   Self-explanatory.

*Dump file is incorrectly formatted*
   Self-explanatory.

*Dump image body yes requires atom style body*
   Self-explanatory.

*Dump image bond not allowed with no bond types*
   Self-explanatory.

*Dump image cannot perform sorting*
   Self-explanatory.

*Dump image line requires atom style line*
   Self-explanatory.

*Dump image persp option is not yet supported*
   Self-explanatory.

*Dump image requires one snapshot per file*
   Use a "\*" in the filename.

*Dump image tri requires atom style tri*
   Self-explanatory.

*Dump local and fix not computed at compatible times*
   The fix must produce per-atom quantities on timesteps that dump local
   needs them.

*Dump local attributes contain no compute or fix*
   Self-explanatory.

*Dump local cannot sort by atom ID*
   This is because dump local does not really dump per-atom info.

*Dump local compute does not calculate local array*
   Self-explanatory.

*Dump local compute does not calculate local vector*
   Self-explanatory.

*Dump local compute does not compute local info*
   Self-explanatory.

*Dump local compute vector is accessed out-of-range*
   Self-explanatory.

*Dump local count is not consistent across input fields*
   Every column of output must be the same length.

*Dump local fix does not compute local array*
   Self-explanatory.

*Dump local fix does not compute local info*
   Self-explanatory.

*Dump local fix does not compute local vector*
   Self-explanatory.

*Dump local fix vector is accessed out-of-range*
   Self-explanatory.

*Dump modify bcolor not allowed with no bond types*
   Self-explanatory.

*Dump modify bdiam not allowed with no bond types*
   Self-explanatory.

*Dump modify compute ID does not compute per-atom array*
   Self-explanatory.

*Dump modify compute ID does not compute per-atom info*
   Self-explanatory.

*Dump modify compute ID does not compute per-atom vector*
   Self-explanatory.

*Dump modify compute ID vector is not large enough*
   Self-explanatory.

*Dump modify element names do not match atom types*
   Number of element names must equal number of atom types.

*Dump modify fix ID does not compute per-atom array*
   Self-explanatory.

*Dump modify fix ID does not compute per-atom info*
   Self-explanatory.

*Dump modify fix ID does not compute per-atom vector*
   Self-explanatory.

*Dump modify fix ID vector is not large enough*
   Self-explanatory.

*Dump modify variable is not atom-style variable*
   Self-explanatory.

*Dump sort column is invalid*
   Self-explanatory.

*Dump xtc requires sorting by atom ID*
   Use the dump_modify sort command to enable this.

*Dump xyz/gz only writes compressed files*
   The dump xyz/gz output file name must have a .gz suffix.

*Dump_modify buffer yes not allowed for this style*
   Self-explanatory.

*Dump_modify format string is too short*
   There are more fields to be dumped in a line of output than your
   format string specifies.

*Dump_modify region ID does not exist*
   Self-explanatory.

*Dumping an atom property that isn't allocated*
   The chosen atom style does not define the per-atom quantity being
   dumped.

*Duplicate atom IDs exist*
   Self-explanatory.

*Duplicate fields in read_dump command*
   Self-explanatory.

*Duplicate particle in PeriDynamic bond - simulation box is too small*
   This is likely because your box length is shorter than 2 times
   the bond length.

*Electronic temperature dropped below zero*
   Something has gone wrong with the fix ttm electron temperature model.

*Element not defined in potential file*
   The specified element is not in the potential file.

*Empty brackets in variable*
   There is no variable syntax that uses empty brackets.  Check
   the variable doc page.

*Energy was not tallied on needed timestep*
   You are using a thermo keyword that requires potentials to
   have tallied energy, but they didn't on this timestep.  See the
   variable doc page for ideas on how to make this work.

*Epsilon or sigma reference not set by pair style in PPPMDisp*
   Self-explanatory.

*Epsilon or sigma reference not set by pair style in ewald/n*
   The pair style is not providing the needed epsilon or sigma values.

*Error in vdw spline: inner radius > outer radius*
   A pre-tabulated spline is invalid.  Likely a problem with the
   potential parameters.

*Error writing averaged chunk data*
   Something in the output to the file triggered an error.

*Error writing file header*
   Something in the output to the file triggered an error.

*Error writing out correlation data*
   Something in the output to the file triggered an error.

*Error writing out histogram data*
   Something in the output to the file triggered an error.

*Error writing out time averaged data*
   Something in the output to the file triggered an error.

*Failed to allocate %ld bytes for array %s*
   Your LAMMPS simulation has run out of memory.  You need to run a
   smaller simulation or on more processors.

*Failed to open FFmpeg pipeline to file %s*
   The specified file cannot be opened.  Check that the path and name are
   correct and writable and that the FFmpeg executable can be found and run.

*Failed to reallocate %ld bytes for array %s*
   Your LAMMPS simulation has run out of memory.  You need to run a
   smaller simulation or on more processors.

*Fewer SRD bins than processors in some dimension*
   This is not allowed.  Make your SRD bin size smaller.

*File variable could not read value*
   Check the file assigned to the variable.

*Final box dimension due to fix deform is < 0.0*
   Self-explanatory.

*Fix %s does not allow use of dynamic group*
   Dynamic groups have not yet been enabled for this fix.

*Fix ID for compute chunk/atom does not exist*
   Self-explanatory.

*Fix ID for compute erotate/rigid does not exist*
   Self-explanatory.

*Fix ID for compute ke/rigid does not exist*
   Self-explanatory.

*Fix ID for compute reduce does not exist*
   Self-explanatory.

*Fix ID for compute slice does not exist*
   Self-explanatory.

*Fix ID for fix ave/atom does not exist*
   Self-explanatory.

*Fix ID for fix ave/chunk does not exist*
   Self-explanatory.

*Fix ID for fix ave/correlate does not exist*
   Self-explanatory.

*Fix ID for fix ave/histo does not exist*
   Self-explanatory.

*Fix ID for fix ave/time does not exist*
   Self-explanatory.

*Fix ID for fix store/state does not exist*
   Self-explanatory

*Fix ID for fix vector does not exist*
   Self-explanatory.

*Fix ID for read_data does not exist*
   Self-explanatory.

*Fix ID for velocity does not exist*
   Self-explanatory.

*Fix ID must be alphanumeric or underscore characters*
   Self-explanatory.

*Fix SRD: bad bin assignment for SRD advection*
   Something has gone wrong in your SRD model; try using more
   conservative settings.

*Fix SRD: bad search bin assignment*
   Something has gone wrong in your SRD model; try using more
   conservative settings.

*Fix SRD: bad stencil bin for big particle*
   Something has gone wrong in your SRD model; try using more
   conservative settings.

*Fix SRD: too many big particles in bin*
   Reset the ATOMPERBIN parameter at the top of fix_srd.cpp
   to a larger value, and re-compile the code.

*Fix SRD: too many walls in bin*
   This should not happen unless your system has been setup incorrectly.

*Fix adapt interface to this pair style not supported*
   New coding for the pair style would need to be done.

*Fix adapt kspace style does not exist*
   Self-explanatory.

*Fix adapt pair style does not exist*
   Self-explanatory

*Fix adapt pair style param not supported*
   The pair style does not know about the parameter you specified.

*Fix adapt requires atom attribute charge*
   The atom style being used does not specify an atom charge.

*Fix adapt requires atom attribute diameter*
   The atom style being used does not specify an atom diameter.

*Fix adapt type pair range is not valid for pair hybrid sub-style*
   Self-explanatory.

*Fix append/atoms requires a lattice be defined*
   Use the lattice command for this purpose.

*Fix ave/atom compute array is accessed out-of-range*
   Self-explanatory.

*Fix ave/atom compute does not calculate a per-atom array*
   Self-explanatory.

*Fix ave/atom compute does not calculate a per-atom vector*
   A compute used by fix ave/atom must generate per-atom values.

*Fix ave/atom compute does not calculate per-atom values*
   A compute used by fix ave/atom must generate per-atom values.

*Fix ave/atom fix array is accessed out-of-range*
   Self-explanatory.

*Fix ave/atom fix does not calculate a per-atom array*
   Self-explanatory.

*Fix ave/atom fix does not calculate a per-atom vector*
   A fix used by fix ave/atom must generate per-atom values.

*Fix ave/atom fix does not calculate per-atom values*
   A fix used by fix ave/atom must generate per-atom values.

*Fix ave/atom variable is not atom-style variable*
   A variable used by fix ave/atom must generate per-atom values.

*Fix ave/chunk compute does not calculate a per-atom array*
   Self-explanatory.

*Fix ave/chunk compute does not calculate a per-atom vector*
   Self-explanatory.

*Fix ave/chunk compute does not calculate per-atom values*
   Self-explanatory.

*Fix ave/chunk compute vector is accessed out-of-range*
   Self-explanatory.

*Fix ave/chunk does not use chunk/atom compute*
   The specified compute is not for a compute chunk/atom command.

*Fix ave/chunk fix does not calculate a per-atom array*
   Self-explanatory.

*Fix ave/chunk fix does not calculate a per-atom vector*
   Self-explanatory.

*Fix ave/chunk fix does not calculate per-atom values*
   Self-explanatory.

*Fix ave/chunk fix vector is accessed out-of-range*
   Self-explanatory.

*Fix ave/chunk variable is not atom-style variable*
   Self-explanatory.

*Fix ave/correlate compute does not calculate a scalar*
   Self-explanatory.

*Fix ave/correlate compute does not calculate a vector*
   Self-explanatory.

*Fix ave/correlate compute vector is accessed out-of-range*
   The index for the vector is out of bounds.

*Fix ave/correlate fix does not calculate a scalar*
   Self-explanatory.

*Fix ave/correlate fix does not calculate a vector*
   Self-explanatory.

*Fix ave/correlate fix vector is accessed out-of-range*
   The index for the vector is out of bounds.

*Fix ave/correlate variable is not equal-style variable*
   Self-explanatory.

*Fix ave/histo cannot input local values in scalar mode*
   Self-explanatory.

*Fix ave/histo cannot input per-atom values in scalar mode*
   Self-explanatory.

*Fix ave/histo compute array is accessed out-of-range*
   Self-explanatory.

*Fix ave/histo compute does not calculate a global array*
   Self-explanatory.

*Fix ave/histo compute does not calculate a global scalar*
   Self-explanatory.

*Fix ave/histo compute does not calculate a global vector*
   Self-explanatory.

*Fix ave/histo compute does not calculate a local array*
   Self-explanatory.

*Fix ave/histo compute does not calculate a local vector*
   Self-explanatory.

*Fix ave/histo compute does not calculate a per-atom array*
   Self-explanatory.

*Fix ave/histo compute does not calculate a per-atom vector*
   Self-explanatory.

*Fix ave/histo compute does not calculate local values*
   Self-explanatory.

*Fix ave/histo compute does not calculate per-atom values*
   Self-explanatory.

*Fix ave/histo compute vector is accessed out-of-range*
   Self-explanatory.

*Fix ave/histo fix array is accessed out-of-range*
   Self-explanatory.

*Fix ave/histo fix does not calculate a global array*
   Self-explanatory.

*Fix ave/histo fix does not calculate a global scalar*
   Self-explanatory.

*Fix ave/histo fix does not calculate a global vector*
   Self-explanatory.

*Fix ave/histo fix does not calculate a local array*
   Self-explanatory.

*Fix ave/histo fix does not calculate a local vector*
   Self-explanatory.

*Fix ave/histo fix does not calculate a per-atom array*
   Self-explanatory.

*Fix ave/histo fix does not calculate a per-atom vector*
   Self-explanatory.

*Fix ave/histo fix does not calculate local values*
   Self-explanatory.

*Fix ave/histo fix does not calculate per-atom values*
   Self-explanatory.

*Fix ave/histo fix vector is accessed out-of-range*
   Self-explanatory.

*Fix ave/histo input is invalid compute*
   Self-explanatory.

*Fix ave/histo input is invalid fix*
   Self-explanatory.

*Fix ave/histo input is invalid variable*
   Self-explanatory.

*Fix ave/histo inputs are not all global, peratom, or local*
   All inputs in a single fix ave/histo command must be of the
   same style.

*Fix ave/histo/weight value and weight vector lengths do not match*
   Self-explanatory.

*Fix ave/time cannot set output array intensive/extensive from these inputs*
   One of more of the vector inputs has individual elements which are
   flagged as intensive or extensive.  Such an input cannot be flagged as
   all intensive/extensive when turned into an array by fix ave/time.

*Fix ave/time cannot use variable with vector mode*
   Variables produce scalar values.

*Fix ave/time columns are inconsistent lengths*
   Self-explanatory.

*Fix ave/time compute array is accessed out-of-range*
   An index for the array is out of bounds.

*Fix ave/time compute does not calculate a scalar*
   Self-explanatory.

*Fix ave/time compute does not calculate a vector*
   Self-explanatory.

*Fix ave/time compute does not calculate an array*
   Self-explanatory.

*Fix ave/time compute vector is accessed out-of-range*
   The index for the vector is out of bounds.

*Fix ave/time fix array cannot be variable length*
   Self-explanatory.

*Fix ave/time fix array is accessed out-of-range*
   An index for the array is out of bounds.

*Fix ave/time fix does not calculate a scalar*
   Self-explanatory.

*Fix ave/time fix does not calculate a vector*
   Self-explanatory.

*Fix ave/time fix does not calculate an array*
   Self-explanatory.

*Fix ave/time fix vector cannot be variable length*
   Self-explanatory.

*Fix ave/time fix vector is accessed out-of-range*
   The index for the vector is out of bounds.

*Fix ave/time variable is not equal-style variable*
   Self-explanatory.

*Fix balance rcb cannot be used with comm_style brick*
   Comm_style tiled must be used instead.

*Fix balance shift string is invalid*
   The string can only contain the characters "x", "y", or "z".

*Fix bond/break needs ghost atoms from further away*
   This is because the fix needs to walk bonds to a certain distance to
   acquire needed info, The comm_modify cutoff command can be used to
   extend the communication range.

*Fix bond/create angle type is invalid*
   Self-explanatory.

*Fix bond/create cutoff is longer than pairwise cutoff*
   This is not allowed because bond creation is done using the
   pairwise neighbor list.

*Fix bond/create dihedral type is invalid*
   Self-explanatory.

*Fix bond/create improper type is invalid*
   Self-explanatory.

*Fix bond/create induced too many angles/dihedrals/impropers per atom*
   See the read_data command for info on using the "extra/angle/per/atom",
   (or dihedral, improper) keywords to allow for additional
   angles, dihedrals, and impropers to be formed.

*Fix bond/create needs ghost atoms from further away*
   This is because the fix needs to walk bonds to a certain distance to
   acquire needed info, The comm_modify cutoff command can be used to
   extend the communication range.

*Fix bond/swap cannot use dihedral or improper styles*
   These styles cannot be defined when using this fix.

*Fix bond/swap requires pair and bond styles*
   Self-explanatory.

*Fix bond/swap requires special_bonds = 0,1,1*
   Self-explanatory.

*Fix box/relax generated negative box length*
   The pressure being applied is likely too large.  Try applying
   it incrementally, to build to the high pressure.

*Fix command before simulation box is defined*
   The fix command cannot be used before a read_data, read_restart, or
   create_box command.

*Fix deform cannot use yz variable with xy*
   The yz setting cannot be a variable if xy deformation is also
   specified.  This is because LAMMPS cannot determine if the yz setting
   will induce a box flip which would be invalid if xy is also changing.

*Fix deform is changing yz too much with xy*
   When both yz and xy are changing, it induces changes in xz if the
   box must flip from one tilt extreme to another.  Thus it is not
   allowed for yz to grow so much that a flip is induced.

*Fix deform tilt factors require triclinic box*
   Cannot deform the tilt factors of a simulation box unless it
   is a triclinic (non-orthogonal) box.

*Fix deform volume setting is invalid*
   Cannot use volume style unless other dimensions are being controlled.

*Fix deposit and fix rigid/small not using same molecule template ID*
   Self-explanatory.

*Fix deposit and fix shake not using same molecule template ID*
   Self-explanatory.

*Fix deposit molecule must have atom types*
   The defined molecule does not specify atom types.

*Fix deposit molecule must have coordinates*
   The defined molecule does not specify coordinates.

*Fix deposit molecule template ID must be same as atom_style template ID*
   When using atom_style template, you cannot deposit molecules that are
   not in that template.

*Fix deposit region cannot be dynamic*
   Only static regions can be used with fix deposit.

*Fix deposit region does not support a bounding box*
   Not all regions represent bounded volumes.  You cannot use
   such a region with the fix deposit command.

*Fix deposit shake fix does not exist*
   Self-explanatory.

*Fix efield requires atom attribute q or mu*
   The atom style defined does not have this attribute.

*Fix efield with dipoles cannot use atom-style variables*
   This option is not supported.

*Fix evaporate molecule requires atom attribute molecule*
   The atom style being used does not define a molecule ID.

*Fix external callback function not set*
   This must be done by an external program in order to use this fix.

*Fix for fix ave/atom not computed at compatible time*
   Fixes generate their values on specific timesteps.  Fix ave/atom is
   requesting a value on a non-allowed timestep.

*Fix for fix ave/chunk not computed at compatible time*
   Fixes generate their values on specific timesteps.  Fix ave/chunk is
   requesting a value on a non-allowed timestep.

*Fix for fix ave/correlate not computed at compatible time*
   Fixes generate their values on specific timesteps.  Fix ave/correlate
   is requesting a value on a non-allowed timestep.

*Fix for fix ave/histo not computed at compatible time*
   Fixes generate their values on specific timesteps.  Fix ave/histo is
   requesting a value on a non-allowed timestep.

*Fix for fix ave/spatial not computed at compatible time*
   Fixes generate their values on specific timesteps.  Fix ave/spatial is
   requesting a value on a non-allowed timestep.

*Fix for fix ave/time not computed at compatible time*
   Fixes generate their values on specific timesteps.  Fix ave/time
   is requesting a value on a non-allowed timestep.

*Fix for fix store/state not computed at compatible time*
   Fixes generate their values on specific timesteps.  Fix store/state
   is requesting a value on a non-allowed timestep.

*Fix for fix vector not computed at compatible time*
   Fixes generate their values on specific timesteps.  Fix vector is
   requesting a value on a non-allowed timestep.

*Fix freeze requires atom attribute torque*
   The atom style defined does not have this attribute.

*Fix gcmc and fix shake not using same molecule template ID*
   Self-explanatory.

*Fix gcmc atom has charge, but atom style does not*
   Self-explanatory.

*Fix gcmc cannot exchange individual atoms belonging to a molecule*
   This is an error since you should not delete only one atom of a
   molecule.  The user has specified atomic (non-molecular) gas
   exchanges, but an atom belonging to a molecule could be deleted.

*Fix gcmc does not (yet) work with atom_style template*
   Self-explanatory.

*Fix gcmc molecule command requires that atoms have molecule attributes*
   Should not choose the gcmc molecule feature if no molecules are being
   simulated. The general molecule flag is off, but gcmc's molecule flag
   is on.

*Fix gcmc molecule has charges, but atom style does not*
   Self-explanatory.

*Fix gcmc molecule must have atom types*
   The defined molecule does not specify atom types.

*Fix gcmc molecule must have coordinates*
   The defined molecule does not specify coordinates.

*Fix gcmc molecule template ID must be same as atom_style template ID*
   When using atom_style template, you cannot insert molecules that are
   not in that template.

*Fix gcmc put atom outside box*
   This should not normally happen.  Contact the developers.

*Fix gcmc ran out of available atom IDs*
   See the setting for tagint in the src/lmptype.h file.

*Fix gcmc ran out of available molecule IDs*
   See the setting for tagint in the src/lmptype.h file.

*Fix gcmc region cannot be dynamic*
   Only static regions can be used with fix gcmc.

*Fix gcmc region does not support a bounding box*
   Not all regions represent bounded volumes.  You cannot use
   such a region with the fix gcmc command.

*Fix gcmc region extends outside simulation box*
   Self-explanatory.

*Fix gcmc shake fix does not exist*
   Self-explanatory.

*Fix gld c coefficients must be >= 0*
   Self-explanatory.

*Fix gld needs more prony series coefficients*
   Self-explanatory.

*Fix gld prony terms must be > 0*
   Self-explanatory.

*Fix gld series type must be pprony for now*
   Self-explanatory.

*Fix gld start temperature must be >= 0*
   Self-explanatory.

*Fix gld stop temperature must be >= 0*
   Self-explanatory.

*Fix gld tau coefficients must be > 0*
   Self-explanatory.

*Fix heat group has no atoms*
   Self-explanatory.

*Fix heat kinetic energy of an atom went negative*
   This will cause the velocity rescaling about to be performed by fix
   heat to be invalid.

*Fix heat kinetic energy went negative*
   This will cause the velocity rescaling about to be performed by fix
   heat to be invalid.

*Fix in variable not computed at compatible time*
   Fixes generate their values on specific timesteps.  The variable is
   requesting the values on a non-allowed timestep.

*Fix langevin angmom is not yet implemented with kokkos*
   This option is not yet available.

*Fix langevin angmom requires atom style ellipsoid*
   Self-explanatory.

*Fix langevin angmom requires extended particles*
   This fix option cannot be used with point particles.

*Fix langevin gjf and respa are not compatible*
   Self-explanatory.

*Fix langevin gjf cannot have period equal to dt/2*
   If the period is equal to dt/2 then division by zero will happen.

*Fix langevin gjf should come before fix nve*
   Self-explanatory.

*Fix langevin gjf with tbias is not yet implemented with kokkos*
   This option is not yet available.

*Fix langevin omega is not yet implemented with kokkos*
   This option is not yet available.

*Fix langevin omega requires atom style sphere*
   Self-explanatory.

*Fix langevin omega requires extended particles*
   One of the particles has radius 0.0.

*Fix langevin period must be > 0.0*
   The time window for temperature relaxation must be > 0

*Fix langevin variable returned negative temperature*
   Self-explanatory.

*Fix momentum group has no atoms*
   Self-explanatory.

*Fix move cannot define z or vz variable for 2d problem*
   Self-explanatory.

*Fix move cannot rotate aroung non z-axis for 2d problem*
   Self-explanatory.

*Fix move cannot set linear z motion for 2d problem*
   Self-explanatory.

*Fix move cannot set wiggle z motion for 2d problem*
   Self-explanatory.

*Fix msst compute ID does not compute potential energy*
   Self-explanatory.

*Fix msst compute ID does not compute pressure*
   Self-explanatory.

*Fix msst compute ID does not compute temperature*
   Self-explanatory.

*Fix msst requires a periodic box*
   Self-explanatory.

*Fix msst tscale must satisfy 0 <= tscale < 1*
   Self-explanatory.

*Fix npt/nph has tilted box too far in one step - periodic cell is too far from equilibrium state*
   Self-explanatory.  The change in the box tilt is too extreme
   on a short timescale.

*Fix numdiff requires an atom map, see atom_modify*
   Self-explanatory. Efficient loop over all atoms for numerical
   difference requires an atom map.

*Fix numdiff requires consecutive atom IDs*
   Self-explanatory. Efficient loop over all atoms for numerical
   difference requires consecutive atom IDs.

*Fix nve/asphere requires extended particles*
   This fix can only be used for particles with a shape setting.

*Fix nve/asphere/noforce requires atom style ellipsoid*
   Self-explanatory.

*Fix nve/asphere/noforce requires extended particles*
   One of the particles is not an ellipsoid.

*Fix nve/body requires atom style body*
   Self-explanatory.

*Fix nve/body requires bodies*
   This fix can only be used for particles that are bodies.

*Fix nve/line can only be used for 2d simulations*
   Self-explanatory.

*Fix nve/line requires atom style line*
   Self-explanatory.

*Fix nve/line requires line particles*
   Self-explanatory.

*Fix nve/sphere dipole requires atom attribute mu*
   An atom style with this attribute is needed.

*Fix nve/sphere requires atom style sphere*
   Self-explanatory.

*Fix nve/sphere requires extended particles*
   This fix can only be used for particles of a finite size.

*Fix nve/tri can only be used for 3d simulations*
   Self-explanatory.

*Fix nve/tri requires atom style tri*
   Self-explanatory.

*Fix nve/tri requires tri particles*
   Self-explanatory.

*Fix nvt/nph/npt asphere requires extended particles*
   The shape setting for a particle in the fix group has shape = 0.0,
   which means it is a point particle.

*Fix nvt/nph/npt body requires bodies*
   Self-explanatory.

*Fix nvt/nph/npt sphere requires atom style sphere*
   Self-explanatory.

*Fix nvt/npt/nph damping parameters must be > 0.0*
   Self-explanatory.

*Fix nvt/npt/nph dilate group ID does not exist*
   Self-explanatory.

*Fix nvt/sphere requires extended particles*
   This fix can only be used for particles of a finite size.

*Fix orient/fcc file open failed*
   The fix orient/fcc command could not open a specified file.

*Fix orient/fcc file read failed*
   The fix orient/fcc command could not read the needed parameters from a
   specified file.

*Fix orient/fcc found self twice*
   The neighbor lists used by fix orient/fcc are messed up.  If this
   error occurs, it is likely a bug, so send an email to the
   `developers <https://lammps.sandia.gov/authors.html>`_.

*Fix peri neigh does not exist*
   Somehow a fix that the pair style defines has been deleted.

*Fix pour and fix rigid/small not using same molecule template ID*
   Self-explanatory.

*Fix pour and fix shake not using same molecule template ID*
   Self-explanatory.

*Fix pour insertion count per timestep is 0*
   Self-explanatory.

*Fix pour molecule must have atom types*
   The defined molecule does not specify atom types.

*Fix pour molecule must have coordinates*
   The defined molecule does not specify coordinates.

*Fix pour molecule template ID must be same as atom style template ID*
   When using atom_style template, you cannot pour molecules that are
   not in that template.

*Fix pour polydisperse fractions do not sum to 1.0*
   Self-explanatory.

*Fix pour region ID does not exist*
   Self-explanatory.

*Fix pour region cannot be dynamic*
   Only static regions can be used with fix pour.

*Fix pour region does not support a bounding box*
   Not all regions represent bounded volumes.  You cannot use
   such a region with the fix pour command.

*Fix pour requires atom attributes radius, rmass*
   The atom style defined does not have these attributes.

*Fix pour rigid fix does not exist*
   Self-explanatory.

*Fix pour shake fix does not exist*
   Self-explanatory.

*Fix press/berendsen damping parameters must be > 0.0*
   Self-explanatory.

*Fix property/atom cannot specify mol twice*
   Self-explanatory.

*Fix property/atom cannot specify q twice*
   Self-explanatory.

*Fix property/atom mol when atom_style already has molecule attribute*
   Self-explanatory.

*Fix property/atom q when atom_style already has charge attribute*
   Self-explanatory.

*Fix property/atom vector name already exists*
   The name for an integer or floating-point vector must be unique.

*Fix qeq has negative upper Taper radius cutoff*
   Self-explanatory.

*Fix qeq/comb group has no atoms*
   Self-explanatory.

*Fix qeq/comb requires atom attribute q*
   An atom style with charge must be used to perform charge equilibration.

*Fix qeq/dynamic group has no atoms*
   Self-explanatory.

*Fix qeq/dynamic requires atom attribute q*
   Self-explanatory.

*Fix qeq/fire group has no atoms*
   Self-explanatory.

*Fix qeq/fire requires atom attribute q*
   Self-explanatory.

*Fix qeq/point group has no atoms*
   Self-explanatory.

*Fix qeq/point has insufficient QEq matrix size*
   Occurs when number of neighbor atoms for an atom increased too much
   during a run.  Increase SAFE_ZONE and MIN_CAP in fix_qeq.h and
   re-compile.

*Fix qeq/point requires atom attribute q*
   Self-explanatory.

*Fix qeq/shielded group has no atoms*
   Self-explanatory.

*Fix qeq/shielded has insufficient QEq matrix size*
   Occurs when number of neighbor atoms for an atom increased too much
   during a run.  Increase SAFE_ZONE and MIN_CAP in fix_qeq.h and
   re-compile.

*Fix qeq/shielded requires atom attribute q*
   Self-explanatory.

*Fix qeq/slater could not extract params from pair coul/streitz*
   This should not happen unless pair coul/streitz has been altered.

*Fix qeq/slater group has no atoms*
   Self-explanatory.

*Fix qeq/slater has insufficient QEq matrix size*
   Occurs when number of neighbor atoms for an atom increased too much
   during a run.  Increase SAFE_ZONE and MIN_CAP in fix_qeq.h and
   re-compile.

*Fix qeq/slater requires atom attribute q*
   Self-explanatory.

*Fix reax/bonds numbonds > nsbmax_most*
   The limit of the number of bonds expected by the ReaxFF force field
   was exceeded.

*Fix recenter group has no atoms*
   Self-explanatory.

*Fix restrain requires an atom map, see atom_modify*
   Self-explanatory.

*Fix rigid atom has non-zero image flag in a non-periodic dimension*
   Image flags for non-periodic dimensions should not be set.

*Fix rigid file has no lines*
   Self-explanatory.

*Fix rigid langevin period must be > 0.0*
   Self-explanatory.

*Fix rigid molecule requires atom attribute molecule*
   Self-explanatory.

*Fix rigid npt/nph dilate group ID does not exist*
   Self-explanatory.

*Fix rigid npt/nph does not yet allow triclinic box*
   This is a current restriction in LAMMPS.

*Fix rigid npt/nph period must be > 0.0*
   Self-explanatory.

*Fix rigid npt/small t_chain should not be less than 1*
   Self-explanatory.

*Fix rigid npt/small t_order must be 3 or 5*
   Self-explanatory.

*Fix rigid nvt/npt/nph damping parameters must be > 0.0*
   Self-explanatory.

*Fix rigid nvt/small t_chain should not be less than 1*
   Self-explanatory.

*Fix rigid nvt/small t_iter should not be less than 1*
   Self-explanatory.

*Fix rigid nvt/small t_order must be 3 or 5*
   Self-explanatory.

*Fix rigid xy torque cannot be on for 2d simulation*
   Self-explanatory.

*Fix rigid z force cannot be on for 2d simulation*
   Self-explanatory.

*Fix rigid/npt period must be > 0.0*
   Self-explanatory.

*Fix rigid/npt temperature order must be 3 or 5*
   Self-explanatory.

*Fix rigid/npt/small period must be > 0.0*
   Self-explanatory.

*Fix rigid/nvt period must be > 0.0*
   Self-explanatory.

*Fix rigid/nvt temperature order must be 3 or 5*
   Self-explanatory.

*Fix rigid/nvt/small period must be > 0.0*
   Self-explanatory.

*Fix rigid/small atom has non-zero image flag in a non-periodic dimension*
   Image flags for non-periodic dimensions should not be set.

*Fix rigid/small langevin period must be > 0.0*
   Self-explanatory.

*Fix rigid/small molecule must have atom types*
   The defined molecule does not specify atom types.

*Fix rigid/small molecule must have coordinates*
   The defined molecule does not specify coordinates.

*Fix rigid/small npt/nph period must be > 0.0*
   Self-explanatory.

*Fix rigid/small nvt/npt/nph damping parameters must be > 0.0*
   Self-explanatory.

*Fix rigid/small nvt/npt/nph dilate group ID does not exist*
   Self-explanatory.

*Fix rigid/small requires an atom map, see atom_modify*
   Self-explanatory.

*Fix rigid/small requires atom attribute molecule*
   Self-explanatory.

*Fix rigid: Bad principal moments*
   The principal moments of inertia computed for a rigid body
   are not within the required tolerances.

*Fix shake cannot be used with minimization*
   Cannot use fix shake while doing an energy minimization since
   it turns off bonds that should contribute to the energy.

*Fix shake molecule template must have shake info*
   The defined molecule does not specify SHAKE information.

*Fix spring couple group ID does not exist*
   Self-explanatory.

*Fix srd can only currently be used with comm_style brick*
   This is a current restriction in LAMMPS.

*Fix srd lamda must be >= 0.6 of SRD grid size*
   This is a requirement for accuracy reasons.

*Fix srd no-slip requires atom attribute torque*
   This is because the SRD collisions will impart torque to the solute
   particles.

*Fix srd requires SRD particles all have same mass*
   Self-explanatory.

*Fix srd requires ghost atoms store velocity*
   Use the comm_modify vel yes command to enable this.

*Fix srd requires newton pair on*
   Self-explanatory.

*Fix store/state compute array is accessed out-of-range*
   Self-explanatory.

*Fix store/state compute does not calculate a per-atom array*
   The compute calculates a per-atom vector.

*Fix store/state compute does not calculate a per-atom vector*
   The compute calculates a per-atom vector.

*Fix store/state compute does not calculate per-atom values*
   Computes that calculate global or local quantities cannot be used
   with fix store/state.

*Fix store/state fix array is accessed out-of-range*
   Self-explanatory.

*Fix store/state fix does not calculate a per-atom array*
   The fix calculates a per-atom vector.

*Fix store/state fix does not calculate a per-atom vector*
   The fix calculates a per-atom array.

*Fix store/state fix does not calculate per-atom values*
   Fixes that calculate global or local quantities cannot be used with
   fix store/state.

*Fix store/state for atom property that isn't allocated*
   Self-explanatory.

*Fix store/state variable is not atom-style variable*
   Only atom-style variables calculate per-atom quantities.

*Fix temp/berendsen period must be > 0.0*
   Self-explanatory.

*Fix temp/berendsen variable returned negative temperature*
   Self-explanatory.

*Fix temp/csld is not compatible with fix rattle or fix shake*
   These two commands cannot currently be used together with fix temp/csld.

*Fix temp/csld variable returned negative temperature*
   Self-explanatory.

*Fix temp/csvr variable returned negative temperature*
   Self-explanatory.

*Fix temp/rescale variable returned negative temperature*
   Self-explanatory.

*Fix tfmc displacement length must be > 0*
   Self-explanatory.

*Fix tfmc is not compatible with fix shake*
   These two commands cannot currently be used together.

*Fix tfmc temperature must be > 0*
   Self-explanatory.

*Fix thermal/conductivity swap value must be positive*
   Self-explanatory.

*Fix tmd must come after integration fixes*
   Any fix tmd command must appear in the input script after all time
   integration fixes (nve, nvt, npt).  See the fix tmd documentation for
   details.

*Fix ttm electron temperatures must be > 0.0*
   Self-explanatory.

*Fix ttm electronic_density must be > 0.0*
   Self-explanatory.

*Fix ttm electronic_specific_heat must be > 0.0*
   Self-explanatory.

*Fix ttm electronic_thermal_conductivity must be >= 0.0*
   Self-explanatory.

*Fix ttm gamma_p must be > 0.0*
   Self-explanatory.

*Fix ttm gamma_s must be >= 0.0*
   Self-explanatory.

*Fix ttm number of nodes must be > 0*
   Self-explanatory.

*Fix ttm v_0 must be >= 0.0*
   Self-explanatory.

*Fix used in compute chunk/atom not computed at compatible time*
   The chunk/atom compute cannot query the output of the fix on a timestep
   it is needed.

*Fix used in compute reduce not computed at compatible time*
   Fixes generate their values on specific timesteps.  Compute reduce is
   requesting a value on a non-allowed timestep.

*Fix used in compute slice not computed at compatible time*
   Fixes generate their values on specific timesteps.  Compute slice is
   requesting a value on a non-allowed timestep.

*Fix vector cannot set output array intensive/extensive from these inputs*
   The inputs to the command have conflicting intensive/extensive attributes.
   You need to use more than one fix vector command.

*Fix vector compute does not calculate a scalar*
   Self-explanatory.

*Fix vector compute does not calculate a vector*
   Self-explanatory.

*Fix vector compute vector is accessed out-of-range*
   Self-explanatory.

*Fix vector fix does not calculate a scalar*
   Self-explanatory.

*Fix vector fix does not calculate a vector*
   Self-explanatory.

*Fix vector fix vector is accessed out-of-range*
   Self-explanatory.

*Fix vector variable is not equal-style variable*
   Self-explanatory.

*Fix viscosity swap value must be positive*
   Self-explanatory.

*Fix viscosity vtarget value must be positive*
   Self-explanatory.

*Fix wall cutoff <= 0.0*
   Self-explanatory.

*Fix wall/colloid requires atom style sphere*
   Self-explanatory.

*Fix wall/colloid requires extended particles*
   One of the particles has radius 0.0.

*Fix wall/gran is incompatible with Pair style*
   Must use a granular pair style to define the parameters needed for
   this fix.

*Fix wall/gran requires atom style sphere*
   Self-explanatory.

*Fix wall/piston command only available at zlo*
   The face keyword must be zlo.

*Fix wall/region colloid requires atom style sphere*
   Self-explanatory.

*Fix wall/region colloid requires extended particles*
   One of the particles has radius 0.0.

*Fix wall/region cutoff <= 0.0*
   Self-explanatory.

*Fix_modify pressure ID does not compute pressure*
   The compute ID assigned to the fix must compute pressure.

*Fix_modify temperature ID does not compute temperature*
   The compute ID assigned to the fix must compute temperature.

*For triclinic deformation, specified target stress must be hydrostatic*
   Triclinic pressure control is allowed using the tri keyword, but
   non-hydrostatic pressure control can not be used in this case.

*Found no restart file matching pattern*
   When using a "\*" in the restart file name, no matching file was found.

*GPU library not compiled for this accelerator*
   Self-explanatory.

*GPU package does not (yet) work with atom_style template*
   Self-explanatory.

*GPU particle split must be set to 1 for this pair style.*
   For this pair style, you cannot run part of the force calculation on
   the host.  See the package command.

*GPUs are requested but Kokkos has not been compiled for CUDA*
   Re-compile Kokkos with CUDA support to use GPUs.

*Ghost velocity forward comm not yet implemented with Kokkos*
   This is a current restriction.

*Gmask function in equal-style variable formula*
   Gmask is per-atom operation.

*Gravity changed since fix pour was created*
   The gravity vector defined by fix gravity must be static.

*Gravity must point in -y to use with fix pour in 2d*
   Self-explanatory.

*Gravity must point in -z to use with fix pour in 3d*
   Self-explanatory.

*Grmask function in equal-style variable formula*
   Grmask is per-atom operation.

*Group ID does not exist*
   A group ID used in the group command does not exist.

*Group ID in variable formula does not exist*
   Self-explanatory.

*Group all cannot be made dynamic*
   This operation is not allowed.

*Group command before simulation box is defined*
   The group command cannot be used before a read_data, read_restart, or
   create_box command.

*Group dynamic cannot reference itself*
   Self-explanatory.

*Group dynamic parent group cannot be dynamic*
   Self-explanatory.

*Group dynamic parent group does not exist*
   Self-explanatory.

*Group region ID does not exist*
   A region ID used in the group command does not exist.

*If read_dump purges it cannot replace or trim*
   These operations are not compatible.  See the read_dump doc
   page for details.

*Illegal ... command*
   Self-explanatory.  Check the input script syntax and compare to the
   documentation for the command.  You can use -echo screen as a
   command-line option when running LAMMPS to see the offending line.

*Illegal COMB parameter*
   One or more of the coefficients defined in the potential file is
   invalid.

*Illegal COMB3 parameter*
   One or more of the coefficients defined in the potential file is
   invalid.

*Illegal Stillinger-Weber parameter*
   One or more of the coefficients defined in the potential file is
   invalid.

*Illegal Tersoff parameter*
   One or more of the coefficients defined in the potential file is
   invalid.

*Illegal Vashishta parameter*
   One or more of the coefficients defined in the potential file is
   invalid.

*Illegal compute voronoi/atom command (occupation and (surface or edges))*
   Self-explanatory.

*Illegal coul/streitz parameter*
   One or more of the coefficients defined in the potential file is
   invalid.

*Illegal dump_modify sfactor value (must be > 0.0)*
   Self-explanatory.

*Illegal dump_modify tfactor value (must be > 0.0)*
   Self-explanatory.

*Illegal fix gcmc gas mass <= 0*
   The computed mass of the designated gas molecule or atom type was less
   than or equal to zero.

*Illegal fix tfmc random seed*
   Seeds can only be nonzero positive integers.

*Illegal fix wall/piston velocity*
   The piston velocity must be positive.

*Illegal integrate style*
   Self-explanatory.

*Illegal nb3b/harmonic parameter*
   One or more of the coefficients defined in the potential file is
   invalid.

*Illegal number of angle table entries*
   There must be at least 2 table entries.

*Illegal number of bond table entries*
   There must be at least 2 table entries.

*Illegal number of pair table entries*
   There must be at least 2 table entries.

*Illegal or unset periodicity in restart*
   This error should not normally occur unless the restart file is invalid.

*Illegal range increment value*
   The increment must be >= 1.

*Illegal simulation box*
   The lower bound of the simulation box is greater than the upper bound.

*Illegal size double vector read requested*
   This error should not normally occur unless the restart file is invalid.

*Illegal size integer vector read requested*
   This error should not normally occur unless the restart file is invalid.

*Illegal size string or corrupt restart*
   This error should not normally occur unless the restart file is invalid.

*Imageint setting in lmptype.h is invalid*
   Imageint must be as large or larger than smallint.

*Imageint setting in lmptype.h is not compatible*
   Format of imageint stored in restart file is not consistent with
   LAMMPS version you are running.  See the settings in src/lmptype.h

*Improper atom missing in delete_bonds*
   The delete_bonds command cannot find one or more atoms in a particular
   improper on a particular processor.  The pairwise cutoff is too short
   or the atoms are too far apart to make a valid improper.

*Improper atom missing in set command*
   The set command cannot find one or more atoms in a particular improper
   on a particular processor.  The pairwise cutoff is too short or the
   atoms are too far apart to make a valid improper.

*Improper atoms %d %d %d %d missing on proc %d at step %ld*
   One or more of 4 atoms needed to compute a particular improper are
   missing on this processor.  Typically this is because the pairwise
   cutoff is set too short or the improper has blown apart and an atom is
   too far away.

*Improper atoms missing on proc %d at step %ld*
   One or more of 4 atoms needed to compute a particular improper are
   missing on this processor.  Typically this is because the pairwise
   cutoff is set too short or the improper has blown apart and an atom is
   too far away.

*Improper coeff for hybrid has invalid style*
   Improper style hybrid uses another improper style as one of its
   coefficients.  The improper style used in the improper_coeff command
   or read from a restart file is not recognized.

*Improper coeffs are not set*
   No improper coefficients have been assigned in the data file or via
   the improper_coeff command.

*Improper style hybrid cannot have hybrid as an argument*
   Self-explanatory.

*Improper style hybrid cannot have none as an argument*
   Self-explanatory.

*Improper style hybrid cannot use same improper style twice*
   Self-explanatory.

*Improper_coeff command before improper_style is defined*
   Coefficients cannot be set in the data file or via the improper_coeff
   command until an improper_style has been assigned.

*Improper_coeff command before simulation box is defined*
   The improper_coeff command cannot be used before a read_data,
   read_restart, or create_box command.

*Improper_coeff command when no impropers allowed*
   The chosen atom style does not allow for impropers to be defined.

*Improper_style command when no impropers allowed*
   The chosen atom style does not allow for impropers to be defined.

*Impropers assigned incorrectly*
   Impropers read in from the data file were not assigned correctly to
   atoms.  This means there is something invalid about the topology
   definitions.

*Impropers defined but no improper types*
   The data file header lists improper but no improper types.

*Incompatible KIM Simulator Model*
   The requested KIM Simulator Model was defined for a different MD code
   and thus is not compatible with LAMMPS.

*Incompatible units for KIM Simulator Model*
   The selected unit style is not compatible with the requested KIM
   Simulator Model.

*Incomplete use of variables in create_atoms command*
   The var and set options must be used together.

*Inconsistent iparam/jparam values in fix bond/create command*
   If itype and jtype are the same, then their maxbond and newtype
   settings must also be the same.

*Inconsistent line segment in data file*
   The end points of the line segment are not equal distances from the
   center point which is the atom coordinate.

*Inconsistent triangle in data file*
   The centroid of the triangle as defined by the corner points is not
   the atom coordinate.

*Inconsistent use of finite-size particles by molecule template molecules*
   Not all of the molecules define a radius for their constituent
   particles.

*Incorrect # of floating-point values in Bodies section of data file*
   See doc page for body style.

*Incorrect # of integer values in Bodies section of data file*
   See doc page for body style.

*Incorrect %s format in data file*
   A section of the data file being read by fix property/atom does
   not have the correct number of values per line.

*Incorrect SNAP parameter file*
   The file cannot be parsed correctly, check its internal syntax.

*Incorrect args for angle coefficients*
   Self-explanatory.  Check the input script or data file.

*Incorrect args for bond coefficients*
   Self-explanatory.  Check the input script or data file.

*Incorrect args for dihedral coefficients*
   Self-explanatory.  Check the input script or data file.

*Incorrect args for improper coefficients*
   Self-explanatory.  Check the input script or data file.

*Incorrect args for pair coefficients*
   Self-explanatory.  Check the input script or data file.

*Incorrect args in pair_style command*
   Self-explanatory.

*Incorrect atom format in data file*
   Number of values per atom line in the data file is not consistent with
   the atom style.

*Incorrect atom format in neb file*
   The number of fields per line is not what expected.

*Incorrect bonus data format in data file*
   See the read_data doc page for a description of how various kinds of
   bonus data must be formatted for certain atom styles.

*Incorrect boundaries with slab Ewald*
   Must have periodic x,y dimensions and non-periodic z dimension to use
   2d slab option with Ewald.

*Incorrect boundaries with slab EwaldDisp*
   Must have periodic x,y dimensions and non-periodic z dimension to use
   2d slab option with Ewald.

*Incorrect boundaries with slab PPPM*
   Must have periodic x,y dimensions and non-periodic z dimension to use
   2d slab option with PPPM.

*Incorrect boundaries with slab PPPMDisp*
   Must have periodic x,y dimensions and non-periodic z dimension to use
   2d slab option with pppm/disp.

*Incorrect conversion in format string*
   A format style variable was not using either a %f, a %g, or a %e conversion.
   Or an immediate variable with format suffix was not using either
   a %f, a %g or a %e conversion in the format suffix.

*Incorrect element names in ADP potential file*
   The element names in the ADP file do not match those requested.

*Incorrect element names in EAM potential file*
   The element names in the EAM file do not match those requested.

*Incorrect format of ... section in data file*
   Number or type of values per line in the given section of the data file
   is not consistent with the requirements for this section.

*Incorrect format in COMB potential file*
   Incorrect number of words per line in the potential file.

*Incorrect format in COMB3 potential file*
   Incorrect number of words per line in the potential file.

*Incorrect format in MEAM potential file*
   Incorrect number of words per line in the potential file.

*Incorrect format in SNAP coefficient file*
   Incorrect number of words per line in the coefficient file.

*Incorrect format in SNAP parameter file*
   Incorrect number of words per line in the parameter file.

*Incorrect format in Stillinger-Weber potential file*
   Incorrect number of words per line in the potential file.

*Incorrect format in TMD target file*
   Format of file read by fix tmd command is incorrect.

*Incorrect format in Tersoff potential file*
   Incorrect number of words per line in the potential file.

*Incorrect format in Vashishta potential file*
   Incorrect number of words per line in the potential file.

*Incorrect format in coul/streitz potential file*
   Incorrect number of words per line in the potential file.

*Incorrect format in nb3b/harmonic potential file*
   Incorrect number of words per line in the potential file.

*Incorrect integer value in Bodies section of data file*
   See doc page for body style.

*Incorrect multiplicity arg for dihedral coefficients*
   Self-explanatory.  Check the input script or data file.

*Incorrect number of elements in potential file*
   Self-explanatory.

*Incorrect rigid body format in fix rigid file*
   The number of fields per line is not what expected.

*Incorrect rigid body format in fix rigid/small file*
   The number of fields per line is not what expected.

*Incorrect sign arg for dihedral coefficients*
   Self-explanatory.  Check the input script or data file.

*Incorrect table format check for element types*
   Self-explanatory.

*Incorrect velocity format in data file*
   Each atom style defines a format for the Velocity section
   of the data file.  The read-in lines do not match.

*Incorrect weight arg for dihedral coefficients*
   Self-explanatory.  Check the input script or data file.

*Index between variable brackets must be positive*
   Self-explanatory.

*Indexed per-atom vector in variable formula without atom map*
   Accessing a value from an atom vector requires the ability to lookup
   an atom index, which is provided by an atom map.  An atom map does not
   exist (by default) for non-molecular problems.  Using the atom_modify
   map command will force an atom map to be created.

*Initial temperatures not all set in fix ttm*
   Self-explanatory.

*Input line quote not followed by white-space*
   An end quote must be followed by white-space.

*Insertion region extends outside simulation box*
   Self-explanatory.

*Insufficient Jacobi rotations for POEMS body*
   Eigensolve for rigid body was not sufficiently accurate.

*Insufficient Jacobi rotations for body nparticle*
   Eigensolve for rigid body was not sufficiently accurate.

*Insufficient Jacobi rotations for rigid body*
   Eigensolve for rigid body was not sufficiently accurate.

*Insufficient Jacobi rotations for rigid molecule*
   Eigensolve for rigid body was not sufficiently accurate.

*Insufficient Jacobi rotations for triangle*
   The calculation of the inertia tensor of the triangle failed.  This
   should not happen if it is a reasonably shaped triangle.

*Insufficient memory on accelerator*
   There is insufficient memory on one of the devices specified for the gpu
   package

*Internal error in atom_style body*
   This error should not occur.  Contact the developers.

*Invalid -reorder N value*
   Self-explanatory.

*Invalid Angles section in molecule file*
   Self-explanatory.

*Invalid Bonds section in molecule file*
   Self-explanatory.

*Invalid Boolean syntax in if command*
   Self-explanatory.

*Invalid Charges section in molecule file*
   Self-explanatory.

*Invalid Coords section in molecule file*
   Self-explanatory.

*Invalid Diameters section in molecule file*
   Self-explanatory.

*Invalid Dihedrals section in molecule file*
   Self-explanatory.

*Invalid Impropers section in molecule file*
   Self-explanatory.

*Invalid Kokkos command-line args*
   Self-explanatory.  See Section 2.7 of the manual for details.

*Invalid LAMMPS restart file*
   The file does not appear to be a LAMMPS restart file since
   it doesn't contain the correct magic string at the beginning.

*Invalid Masses section in molecule file*
   Self-explanatory.

*Invalid molecule ID in molecule file*
   Molecule ID must be a non-zero positive integer.

*Invalid Molecules section in molecule file*
   Self-explanatory.

*Invalid REAX atom type*
   There is a mis-match between LAMMPS atom types and the elements
   listed in the ReaxFF force field file.

*Invalid Special Bond Counts section in molecule file*
   Self-explanatory.

*Invalid Types section in molecule file*
   Self-explanatory.

*Invalid angle count in molecule file*
   Self-explanatory.

*Invalid angle table length*
   Length must be 2 or greater.

*Invalid angle type in Angles section of data file*
   Angle type must be positive integer and within range of specified angle
   types.

*Invalid angle type in Angles section of molecule file*
   Self-explanatory.

*Invalid angle type index for fix shake*
   Self-explanatory.

*Invalid args for non-hybrid pair coefficients*
   "NULL" is only supported in pair_coeff calls when using pair hybrid

*Invalid argument to factorial %d*
   N must be >= 0 and <= 167, otherwise the factorial result is too
   large.

*Invalid atom ID in %s section of data file*
   An atom in a section of the data file being read by fix property/atom
   has an invalid atom ID that is <= 0 or > the maximum existing atom ID.

*Invalid atom ID in Angles section of data file*
   Atom IDs must be positive integers and within range of defined
   atoms.

*Invalid atom ID in Angles section of molecule file*
   Self-explanatory.

*Invalid atom ID in Atoms section of data file*
   Atom IDs must be positive integers.

*Invalid atom ID in Bodies section of data file*
   Atom IDs must be positive integers and within range of defined
   atoms.

*Invalid atom ID in Bonds section of data file*
   Atom IDs must be positive integers and within range of defined
   atoms.

*Invalid atom ID in Bonds section of molecule file*
   Self-explanatory.

*Invalid atom ID in Bonus section of data file*
   Atom IDs must be positive integers and within range of defined
   atoms.

*Invalid atom ID in Dihedrals section of data file*
   Atom IDs must be positive integers and within range of defined
   atoms.

*Invalid atom ID in Fragments section of molecule file*
   Self-explanatory.

*Invalid atom ID in Impropers section of data file*
   Atom IDs must be positive integers and within range of defined
   atoms.

*Invalid atom ID in Velocities section of data file*
   Atom IDs must be positive integers and within range of defined
   atoms.

*Invalid atom ID in dihedrals section of molecule file*
   Self-explanatory.

*Invalid atom ID in impropers section of molecule file*
   Self-explanatory.

*Invalid atom ID in variable file*
   Self-explanatory.

*Invalid atom IDs in neb file*
   An ID in the file was not found in the system.

*Invalid atom diameter in molecule file*
   Diameters must be >= 0.0.

*Invalid atom mass for fix shake*
   Mass specified in fix shake command must be > 0.0.

*Invalid atom mass in molecule file*
   Masses must be > 0.0.

*Invalid atom type in Atoms section of data file*
   Atom types must range from 1 to specified # of types.

*Invalid atom type in create_atoms command*
   The create_box command specified the range of valid atom types.
   An invalid type is being requested.

*Invalid atom type in create_atoms mol command*
   The atom types in the defined molecule are added to the value
   specified in the create_atoms command, as an offset.  The final value
   for each atom must be between 1 to N, where N is the number of atom
   types.

*Invalid atom type in fix atom/swap command*
   The atom type specified in the atom/swap command does not exist.

*Invalid atom type in fix bond/create command*
   Self-explanatory.

*Invalid atom type in fix deposit command*
   Self-explanatory.

*Invalid atom type in fix deposit mol command*
   The atom types in the defined molecule are added to the value
   specified in the create_atoms command, as an offset.  The final value
   for each atom must be between 1 to N, where N is the number of atom
   types.

*Invalid atom type in fix gcmc command*
   The atom type specified in the gcmc command does not exist.

*Invalid atom type in fix pour command*
   Self-explanatory.

*Invalid atom type in fix pour mol command*
   The atom types in the defined molecule are added to the value
   specified in the create_atoms command, as an offset.  The final value
   for each atom must be between 1 to N, where N is the number of atom
   types.

*Invalid atom type in molecule file*
   Atom types must range from 1 to specified # of types.

*Invalid atom type in neighbor exclusion list*
   Atom types must range from 1 to Ntypes inclusive.

*Invalid atom type index for fix shake*
   Atom types must range from 1 to Ntypes inclusive.

*Invalid atom types in pair_write command*
   Atom types must range from 1 to Ntypes inclusive.

*Invalid atom vector in variable formula*
   The atom vector is not recognized.

*Invalid atom_style body command*
   No body style argument was provided.

*Invalid atom_style command*
   Self-explanatory.

*Invalid attribute in dump custom command*
   Self-explanatory.

*Invalid attribute in dump local command*
   Self-explanatory.

*Invalid attribute in dump modify command*
   Self-explanatory.

*Invalid basis setting in create_atoms command*
   The basis index must be between 1 to N where N is the number of basis
   atoms in the lattice.  The type index must be between 1 to N where N
   is the number of atom types.

*Invalid basis setting in fix append/atoms command*
   The basis index must be between 1 to N where N is the number of basis
   atoms in the lattice.  The type index must be between 1 to N where N
   is the number of atom types.

*Invalid bin bounds in compute chunk/atom*
   The lo/hi values are inconsistent.

*Invalid bin bounds in fix ave/spatial*
   The lo/hi values are inconsistent.

*Invalid body nparticle command*
   Arguments in atom-style command are not correct.

*Invalid bond count in molecule file*
   Self-explanatory.

*Invalid bond table length*
   Length must be 2 or greater.

*Invalid bond type in Bonds section of data file*
   Bond type must be positive integer and within range of specified bond
   types.

*Invalid bond type in Bonds section of molecule file*
   Self-explanatory.

*Invalid bond type in create_bonds command*
   Self-explanatory.

*Invalid bond type in fix bond/break command*
   Self-explanatory.

*Invalid bond type in fix bond/create command*
   Self-explanatory.

*Invalid bond type index for fix shake*
   Self-explanatory.  Check the fix shake command in the input script.

*Invalid coeffs for this dihedral style*
   Cannot set class 2 coeffs in data file for this dihedral style.

*Invalid color in dump_modify command*
   The specified color name was not in the list of recognized colors.
   See the dump_modify doc page.

*Invalid color map min/max values*
   The min/max values are not consistent with either each other or
   with values in the color map.

*Invalid command-line argument*
   One or more command-line arguments is invalid.  Check the syntax of
   the command you are using to launch LAMMPS.

*Invalid compute ID in variable formula*
   The compute is not recognized.

*Invalid create_atoms rotation vector for 2d model*
   The rotation vector can only have a z component.

*Invalid custom OpenCL parameter string.*
   There are not enough or too many parameters in the custom string for package
   GPU.

*Invalid cutoff in comm_modify command*
   Specified cutoff must be >= 0.0.

*Invalid cutoffs in pair_write command*
   Inner cutoff must be larger than 0.0 and less than outer cutoff.

*Invalid d1 or d2 value for pair colloid coeff*
   Neither d1 or d2 can be < 0.

*Invalid data file section: Angle Coeffs*
   Atom style does not allow angles.

*Invalid data file section: AngleAngle Coeffs*
   Atom style does not allow impropers.

*Invalid data file section: AngleAngleTorsion Coeffs*
   Atom style does not allow dihedrals.

*Invalid data file section: AngleTorsion Coeffs*
   Atom style does not allow dihedrals.

*Invalid data file section: Angles*
   Atom style does not allow angles.

*Invalid data file section: Bodies*
   Atom style does not allow bodies.

*Invalid data file section: Bond Coeffs*
   Atom style does not allow bonds.

*Invalid data file section: BondAngle Coeffs*
   Atom style does not allow angles.

*Invalid data file section: BondBond Coeffs*
   Atom style does not allow angles.

*Invalid data file section: BondBond13 Coeffs*
   Atom style does not allow dihedrals.

*Invalid data file section: Bonds*
   Atom style does not allow bonds.

*Invalid data file section: Dihedral Coeffs*
   Atom style does not allow dihedrals.

*Invalid data file section: Dihedrals*
   Atom style does not allow dihedrals.

*Invalid data file section: Ellipsoids*
   Atom style does not allow ellipsoids.

*Invalid data file section: EndBondTorsion Coeffs*
   Atom style does not allow dihedrals.

*Invalid data file section: Improper Coeffs*
   Atom style does not allow impropers.

*Invalid data file section: Impropers*
   Atom style does not allow impropers.

*Invalid data file section: Lines*
   Atom style does not allow lines.

*Invalid data file section: MiddleBondTorsion Coeffs*
   Atom style does not allow dihedrals.

*Invalid data file section: Triangles*
   Atom style does not allow triangles.

*Invalid delta_conf in tad command*
   The value must be between 0 and 1 inclusive.

*Invalid density in Atoms section of data file*
   Density value cannot be <= 0.0.

*Invalid density in set command*
   Density must be > 0.0.

*Invalid diameter in set command*
   Self-explanatory.

*Invalid dihedral count in molecule file*
   Self-explanatory.

*Invalid dihedral type in Dihedrals section of data file*
   Dihedral type must be positive integer and within range of specified
   dihedral types.

*Invalid dihedral type in dihedrals section of molecule file*
   Self-explanatory.

*Invalid dipole length in set command*
   Self-explanatory.

*Invalid displace_atoms rotate axis for 2d*
   Axis must be in z direction.

*Invalid dump dcd filename*
   Filenames used with the dump dcd style cannot be binary or compressed
   or cause multiple files to be written.

*Invalid dump frequency*
   Dump frequency must be 1 or greater.

*Invalid dump image element name*
   The specified element name was not in the standard list of elements.
   See the dump_modify doc page.

*Invalid dump image filename*
   The file produced by dump image cannot be binary and must
   be for a single processor.

*Invalid dump image persp value*
   Persp value must be >= 0.0.

*Invalid dump image theta value*
   Theta must be between 0.0 and 180.0 inclusive.

*Invalid dump image zoom value*
   Zoom value must be > 0.0.

*Invalid dump movie filename*
   The file produced by dump movie cannot be binary or compressed
   and must be a single file for a single processor.

*Invalid dump xtc filename*
   Filenames used with the dump xtc style cannot be binary or compressed
   or cause multiple files to be written.

*Invalid dump xyz filename*
   Filenames used with the dump xyz style cannot be binary or cause files
   to be written by each processor.

*Invalid dump_modify threshold operator*
   Operator keyword used for threshold specification in not recognized.

*Invalid entry in -reorder file*
   Self-explanatory.

*Invalid fix ID in variable formula*
   The fix is not recognized.

*Invalid fix ave/time off column*
   Self-explanatory.

*Invalid fix box/relax command for a 2d simulation*
   Fix box/relax styles involving the z dimension cannot be used in
   a 2d simulation.

*Invalid fix box/relax command pressure settings*
   If multiple dimensions are coupled, those dimensions must be specified.

*Invalid fix box/relax pressure settings*
   Settings for coupled dimensions must be the same.

*Invalid fix nvt/npt/nph command for a 2d simulation*
   Cannot control z dimension in a 2d model.

*Invalid fix nvt/npt/nph command pressure settings*
   If multiple dimensions are coupled, those dimensions must be
   specified.

*Invalid fix nvt/npt/nph pressure settings*
   Settings for coupled dimensions must be the same.

*Invalid fix press/berendsen for a 2d simulation*
   The z component of pressure cannot be controlled for a 2d model.

*Invalid fix press/berendsen pressure settings*
   Settings for coupled dimensions must be the same.

*Invalid fix qeq parameter file*
   Element index > number of atom types.

*Invalid fix rigid npt/nph command for a 2d simulation*
   Cannot control z dimension in a 2d model.

*Invalid fix rigid npt/nph command pressure settings*
   If multiple dimensions are coupled, those dimensions must be
   specified.

*Invalid fix rigid/small npt/nph command for a 2d simulation*
   Cannot control z dimension in a 2d model.

*Invalid fix rigid/small npt/nph command pressure settings*
   If multiple dimensions are coupled, those dimensions must be
   specified.

*Invalid flag in force field section of restart file*
   Unrecognized entry in restart file.

*Invalid flag in header section of restart file*
   Unrecognized entry in restart file.

*Invalid flag in peratom section of restart file*
   The format of this section of the file is not correct.

*Invalid flag in type arrays section of restart file*
   Unrecognized entry in restart file.

*Invalid frequency in temper command*
   Nevery must be > 0.

*Invalid group ID in neigh_modify command*
   A group ID used in the neigh_modify command does not exist.

*Invalid group function in variable formula*
   Group function is not recognized.

*Invalid group in comm_modify command*
   Self-explanatory.

*Invalid image up vector*
   Up vector cannot be (0,0,0).

*Invalid immediate variable*
   Syntax of immediate value is incorrect.

*Invalid improper count in molecule file*
   Self-explanatory.

*Invalid improper type in Impropers section of data file*
   Improper type must be positive integer and within range of specified
   improper types.

*Invalid improper type in impropers section of molecule file*
   Self-explanatory.

*Invalid index for non-body particles in compute body/local command*
   Only indices 1,2,3 can be used for non-body particles.

*Invalid index in compute body/local command*
   Self-explanatory.

*Invalid is_active() function in variable formula*
   Self-explanatory.

*Invalid is_available() function in variable formula*
   Self-explanatory.

*Invalid is_defined() function in variable formula*
   Self-explanatory.

*Invalid keyword in angle table parameters*
   Self-explanatory.

*Invalid keyword in bond table parameters*
   Self-explanatory.

*Invalid keyword in compute angle/local command*
   Self-explanatory.

*Invalid keyword in compute bond/local command*
   Self-explanatory.

*Invalid keyword in compute dihedral/local command*
   Self-explanatory.

*Invalid keyword in compute improper/local command*
   Self-explanatory.

*Invalid keyword in compute pair/local command*
   Self-explanatory.

*Invalid keyword in compute property/atom command*
   Self-explanatory.

*Invalid keyword in compute property/chunk command*
   Self-explanatory.

*Invalid keyword in compute property/local command*
   Self-explanatory.

*Invalid keyword in dump cfg command*
   Self-explanatory.

*Invalid keyword in pair table parameters*
   Keyword used in list of table parameters is not recognized.

*Invalid length in set command*
   Self-explanatory.

*Invalid mass in set command*
   Self-explanatory.

*Invalid mass line in data file*
   Self-explanatory.

*Invalid mass value*
   Self-explanatory.

*Invalid math function in variable formula*
   Self-explanatory.

*Invalid math/group/special function in variable formula*
   Self-explanatory.

*Invalid option in lattice command for non-custom style*
   Certain lattice keywords are not supported unless the
   lattice style is "custom".

*Invalid order of forces within respa levels*
   For respa, ordering of force computations within respa levels must
   obey certain rules.  E.g. bonds cannot be compute less frequently than
   angles, pairwise forces cannot be computed less frequently than
   kspace, etc.

*Invalid pair table cutoff*
   Cutoffs in pair_coeff command are not valid with read-in pair table.

*Invalid pair table length*
   Length of read-in pair table is invalid

*Invalid param file for fix qeq/shielded*
   Invalid value of gamma.

*Invalid param file for fix qeq/slater*
   Zeta value is 0.0.

*Invalid partitions in processors part command*
   Valid partitions are numbered 1 to N and the sender and receiver
   cannot be the same partition.

*Invalid python command*
   Self-explanatory.  Check the input script syntax and compare to the
   documentation for the command.  You can use -echo screen as a
   command-line option when running LAMMPS to see the offending line.

*Invalid radius in Atoms section of data file*
   Radius must be >= 0.0.

*Invalid random number seed in fix ttm command*
   Random number seed must be > 0.

*Invalid random number seed in set command*
   Random number seed must be > 0.

*Invalid replace values in compute reduce*
   Self-explanatory.

*Invalid rigid body ID in fix rigid file*
   The ID does not match the number of an existing ID of rigid bodies
   that are defined by the fix rigid command.

*Invalid rigid body ID in fix rigid/small file*
   The ID does not match the number of an existing ID of rigid bodies
   that are defined by the fix rigid/small command.

*Invalid run command N value*
   The number of timesteps must fit in a 32-bit integer.  If you want to
   run for more steps than this, perform multiple shorter runs.

*Invalid run command start/stop value*
   Self-explanatory.

*Invalid run command upto value*
   Self-explanatory.

*Invalid seed for Marsaglia random # generator*
   The initial seed for this random number generator must be a positive
   integer less than or equal to 900 million.

*Invalid seed for Park random # generator*
   The initial seed for this random number generator must be a positive
   integer.

*Invalid shake angle type in molecule file*
   Self-explanatory.

*Invalid shake atom in molecule file*
   Self-explanatory.

*Invalid shake bond type in molecule file*
   Self-explanatory.

*Invalid shake flag in molecule file*
   Self-explanatory.

*Invalid shape in Ellipsoids section of data file*
   Self-explanatory.

*Invalid shape in Triangles section of data file*
   Two or more of the triangle corners are duplicate points.

*Invalid shape in set command*
   Self-explanatory.

*Invalid shear direction for fix wall/gran*
   Self-explanatory.

*Invalid special atom index in molecule file*
   Self-explanatory.

*Invalid special function in variable formula*
   Self-explanatory.

*Invalid style in pair_write command*
   Self-explanatory.  Check the input script.

*Invalid syntax in variable formula*
   Self-explanatory.

*Invalid t_event in prd command*
   Self-explanatory.

*Invalid t_event in tad command*
   The value must be greater than 0.

*Invalid template atom in Atoms section of data file*
   The atom indices must be between 1 to N, where N is the number of
   atoms in the template molecule the atom belongs to.

*Invalid template index in Atoms section of data file*
   The template indices must be between 1 to N, where N is the number of
   molecules in the template.

*Invalid thermo keyword in variable formula*
   The keyword is not recognized.

*Invalid threads_per_atom specified.*
   For 3-body potentials on the GPU, the threads_per_atom setting cannot be
   greater than 4 for NVIDIA GPUs.

*Invalid timestep reset for fix ave/atom*
   Resetting the timestep has invalidated the sequence of timesteps this
   fix needs to process.

*Invalid timestep reset for fix ave/chunk*
   Resetting the timestep has invalidated the sequence of timesteps this
   fix needs to process.

*Invalid timestep reset for fix ave/correlate*
   Resetting the timestep has invalidated the sequence of timesteps this
   fix needs to process.

*Invalid timestep reset for fix ave/histo*
   Resetting the timestep has invalidated the sequence of timesteps this
   fix needs to process.

*Invalid timestep reset for fix ave/spatial*
   Resetting the timestep has invalidated the sequence of timesteps this
   fix needs to process.

*Invalid timestep reset for fix ave/time*
   Resetting the timestep has invalidated the sequence of timesteps this
   fix needs to process.

*Invalid tmax in tad command*
   The value must be greater than 0.0.

*Invalid type for mass set*
   Mass command must set a type from 1-N where N is the number of atom
   types.

*Invalid use of library file() function*
   This function is called through the library interface.  This
   error should not occur.  Contact the developers if it does.

*Invalid value in set command*
   The value specified for the setting is invalid, likely because it is
   too small or too large.

*Invalid variable evaluation in variable formula*
   A variable used in a formula could not be evaluated.

*Invalid variable in next command*
   Self-explanatory.

*Invalid variable name*
   Variable name used in an input script line is invalid.

*Invalid variable name in variable formula*
   Variable name is not recognized.

*Invalid variable style in special function next*
   Only file-style or atomfile-style variables can be used with next().

*Invalid variable style with next command*
   Variable styles *equal* and *world* cannot be used in a next
   command.

*Invalid volume in set command*
   Volume must be > 0.0.

*Invalid wiggle direction for fix wall/gran*
   Self-explanatory.

*Invoked angle equil angle on angle style none*
   Self-explanatory.

*Invoked angle single on angle style none*
   Self-explanatory.

*Invoked bond equil distance on bond style none*
   Self-explanatory.

*Invoked bond single on bond style none*
   Self-explanatory.

*Invoked pair single on pair style none*
   A command (e.g. a dump) attempted to invoke the single() function on a
   pair style none, which is illegal.  You are probably attempting to
   compute per-atom quantities with an undefined pair style.

*Invoking coulombic in pair style lj/coul requires atom attribute q*
   The atom style defined does not have this attribute.

*Invoking coulombic in pair style lj/long/dipole/long requires atom attribute q*
   The atom style defined does not have these attributes.

*KIM Simulator Model has no Model definition*
   There is no model definition (key: model-defn) in the KIM Simulator
   Model.  Please contact the OpenKIM database maintainers to verify
   and potentially correct this.

*KOKKOS package does not yet support comm_style tiled*
   Self-explanatory.

*KOKKOS package requires a kokkos enabled atom_style*
   Self-explanatory.

*KSpace accuracy must be > 0*
   The kspace accuracy designated in the input must be greater than zero.

*KSpace accuracy too large to estimate G vector*
   Reduce the accuracy request or specify gewald explicitly
   via the kspace_modify command.

*KSpace accuracy too low*
   Requested accuracy must be less than 1.0.

*KSpace solver requires a pair style*
   No pair style is defined.

*KSpace style does not yet support triclinic geometries*
   The specified kspace style does not allow for non-orthogonal
   simulation boxes.

*KSpace style has not yet been set*
   Cannot use kspace_modify command until a kspace style is set.

*KSpace style is incompatible with Pair style*
   Setting a kspace style requires that a pair style with matching
   long-range Coulombic or dispersion components be used.

*Keyword %s in MEAM parameter file not recognized*
   Self-explanatory.

*Kokkos has been compiled for CUDA but no GPUs are requested*
   One or more GPUs must be used when Kokkos is compiled for CUDA.

*Kspace_modify mesh parameter must be all zero or all positive*
   Valid kspace mesh parameters are >0. The code will try to auto-detect
   suitable values when all three mesh sizes are set to zero (the default).

*Kspace_modify mesh/disp parameter must be all zero or all positive*
   Valid kspace mesh/disp parameters are >0. The code will try to auto-detect
   suitable values when all three mesh sizes are set to zero **and**
   the required accuracy via *force/disp/real* as well as
   *force/disp/kspace* is set.

*Kspace style does not support compute group/group*
   Self-explanatory.

*Kspace style pppm/disp/tip4p requires newton on*
   Self-explanatory.

*Kspace style pppm/tip4p requires newton on*
   Self-explanatory.

*Kspace style requires atom attribute q*
   The atom style defined does not have these attributes.

*Kspace_modify eigtol must be smaller than one*
   Self-explanatory.

*LAMMPS is not built with Python embedded*
   This is done by including the PYTHON package before LAMMPS is built.
   This is required to use python-style variables.

*LAMMPS unit_style lj not supported by KIM models*
   Self-explanatory. Check the input script or data file.

*LJ6 off not supported in pair_style buck/long/coul/long*
   Self-explanatory.

*Label wasn't found in input script*
   Self-explanatory.

*Lattice orient vectors are not orthogonal*
   The three specified lattice orientation vectors must be mutually
   orthogonal.

*Lattice orient vectors are not right-handed*
   The three specified lattice orientation vectors must create a
   right-handed coordinate system such that a1 cross a2 = a3.

*Lattice primitive vectors are collinear*
   The specified lattice primitive vectors do not for a unit cell with
   non-zero volume.

*Lattice settings are not compatible with 2d simulation*
   One or more of the specified lattice vectors has a non-zero z
   component.

*Lattice spacings are invalid*
   Each x,y,z spacing must be > 0.

*Lattice style incompatible with simulation dimension*
   2d simulation can use sq, sq2, or hex lattice.  3d simulation can use
   sc, bcc, or fcc lattice.

*Log of zero/negative value in variable formula*
   Self-explanatory.

*Lost atoms via balance: original %ld current %ld*
   This should not occur.  Report the problem to the developers.

*Lost atoms: original %ld current %ld*
   Lost atoms are checked for each time thermo output is done.  See the
   thermo_modify lost command for options.  Lost atoms usually indicate
   bad dynamics, e.g. atoms have been blown far out of the simulation
   box, or moved further than one processor's sub-domain away before
   reneighboring.

*MEAM library error %d*
   A call to the MEAM Fortran library returned an error.

*MPI_LMP_BIGINT and bigint in lmptype.h are not compatible*
   The size of the MPI datatype does not match the size of a bigint.

*MPI_LMP_TAGINT and tagint in lmptype.h are not compatible*
   The size of the MPI datatype does not match the size of a tagint.

*MSM can only currently be used with comm_style brick*
   This is a current restriction in LAMMPS.

*MSM grid is too large*
   The global MSM grid is larger than OFFSET in one or more dimensions.
   OFFSET is currently set to 16384.  You likely need to decrease the
   requested accuracy.

*MSM order must be 4, 6, 8, or 10*
   This is a limitation of the MSM implementation in LAMMPS:
   the MSM order can only be 4, 6, 8, or 10.

*Mass command before simulation box is defined*
   The mass command cannot be used before a read_data, read_restart, or
   create_box command.

*Matrix factorization to split dispersion coefficients failed*
   This should not normally happen.  Contact the developers.

*Min_style command before simulation box is defined*
   The min_style command cannot be used before a read_data, read_restart,
   or create_box command.

*Minimization could not find thermo_pe compute*
   This compute is created by the thermo command.  It must have been
   explicitly deleted by a uncompute command.

*Minimize command before simulation box is defined*
   The minimize command cannot be used before a read_data, read_restart,
   or create_box command.

*Mismatched brackets in variable*
   Self-explanatory.

*Mismatched compute in variable formula*
   A compute is referenced incorrectly or a compute that produces per-atom
   values is used in an equal-style variable formula.

*Mismatched fix in variable formula*
   A fix is referenced incorrectly or a fix that produces per-atom
   values is used in an equal-style variable formula.

*Mismatched variable in variable formula*
   A variable is referenced incorrectly or an atom-style variable that
   produces per-atom values is used in an equal-style variable
   formula.

*Modulo 0 in variable formula*
   Self-explanatory.

*Molecule IDs too large for compute chunk/atom*
   The IDs must not be larger than can be stored in a 32-bit integer
   since chunk IDs are 32-bit integers.

*Molecule auto special bond generation overflow*
   Counts exceed maxspecial setting for other atoms in system.

*Molecule file has angles but no nangles setting*
   Self-explanatory.

*Molecule file has body params but no setting for them*
   Self-explanatory.

*Molecule file has bonds but no nbonds setting*
   Self-explanatory.

*Molecule file has dihedrals but no ndihedrals setting*
   Self-explanatory.

*Molecule file has impropers but no nimpropers setting*
   Self-explanatory.

*Molecule file has no Body Doubles section*
   Self-explanatory.

*Molecule file has no Body Integers section*
   Self-explanatory.

*Molecule file has special flags but no bonds*
   Self-explanatory.

*Molecule file needs both Special Bond sections*
   Self-explanatory.

*Molecule file requires atom style body*
   Self-explanatory.

*Molecule file shake flags not before shake atoms*
   The order of the two sections is important.

*Molecule file shake flags not before shake bonds*
   The order of the two sections is important.

*Molecule file shake info is incomplete*
   All 3 SHAKE sections are needed.

*Molecule file special list does not match special count*
   The number of values in an atom's special list does not match count.

*Molecule file z center-of-mass must be 0.0 for 2d*
   Self-explanatory.

*Molecule file z coord must be 0.0 for 2d*
   Self-explanatory.

*Molecule natoms must be 1 for body particle*
   Self-explanatory.

*Molecule sizescale must be 1.0 for body particle*
   Self-explanatory.

*Molecule template ID for atom_style template does not exist*
   Self-explanatory.

*Molecule template ID for create_atoms does not exist*
   Self-explanatory.

*Molecule template ID for fix deposit does not exist*
   Self-explanatory.

*Molecule template ID for fix gcmc does not exist*
   Self-explanatory.

*Molecule template ID for fix pour does not exist*
   Self-explanatory.

*Molecule template ID for fix rigid/small does not exist*
   Self-explanatory.

*Molecule template ID for fix shake does not exist*
   Self-explanatory.

*Molecule template ID must be alphanumeric or underscore characters*
   Self-explanatory.

*Molecule topology/atom exceeds system topology/atom*
   The number of bonds, angles, etc per-atom in the molecule exceeds the
   system setting.  See the create_box command for how to specify these
   values.

*Molecule topology type exceeds system topology type*
   The number of bond, angle, etc types in the molecule exceeds the
   system setting.  See the create_box command for how to specify these
   values.

*More than one fix deform*
   Only one fix deform can be defined at a time.

*More than one fix freeze*
   Only one of these fixes can be defined, since the granular pair
   potentials access it.

*More than one fix shake*
   Only one fix shake can be defined.

*Mu not allowed when not using semi-grand in fix atom/swap command*
   Self-explanatory.

*Must define angle_style before Angle Coeffs*
   Must use an angle_style command before reading a data file that
   defines Angle Coeffs.

*Must define angle_style before BondAngle Coeffs*
   Must use an angle_style command before reading a data file that
   defines Angle Coeffs.

*Must define angle_style before BondBond Coeffs*
   Must use an angle_style command before reading a data file that
   defines Angle Coeffs.

*Must define bond_style before Bond Coeffs*
   Must use a bond_style command before reading a data file that
   defines Bond Coeffs.

*Must define dihedral_style before AngleAngleTorsion Coeffs*
   Must use a dihedral_style command before reading a data file that
   defines AngleAngleTorsion Coeffs.

*Must define dihedral_style before AngleTorsion Coeffs*
   Must use a dihedral_style command before reading a data file that
   defines AngleTorsion Coeffs.

*Must define dihedral_style before BondBond13 Coeffs*
   Must use a dihedral_style command before reading a data file that
   defines BondBond13 Coeffs.

*Must define dihedral_style before Dihedral Coeffs*
   Must use a dihedral_style command before reading a data file that
   defines Dihedral Coeffs.

*Must define dihedral_style before EndBondTorsion Coeffs*
   Must use a dihedral_style command before reading a data file that
   defines EndBondTorsion Coeffs.

*Must define dihedral_style before MiddleBondTorsion Coeffs*
   Must use a dihedral_style command before reading a data file that
   defines MiddleBondTorsion Coeffs.

*Must define improper_style before AngleAngle Coeffs*
   Must use an improper_style command before reading a data file that
   defines AngleAngle Coeffs.

*Must define improper_style before Improper Coeffs*
   Must use an improper_style command before reading a data file that
   defines Improper Coeffs.

*Must define pair_style before Pair Coeffs*
   Must use a pair_style command before reading a data file that defines
   Pair Coeffs.

*Must define pair_style before PairIJ Coeffs*
   Must use a pair_style command before reading a data file that defines
   PairIJ Coeffs.

*Must have more than one processor partition to temper*
   Cannot use the temper command with only one processor partition.  Use
   the -partition command-line option.

*Must read Atoms before Angles*
   The Atoms section of a data file must come before an Angles section.

*Must read Atoms before Bodies*
   The Atoms section of a data file must come before a Bodies section.

*Must read Atoms before Bonds*
   The Atoms section of a data file must come before a Bonds section.

*Must read Atoms before Dihedrals*
   The Atoms section of a data file must come before a Dihedrals section.

*Must read Atoms before Ellipsoids*
   The Atoms section of a data file must come before a Ellipsoids
   section.

*Must read Atoms before Impropers*
   The Atoms section of a data file must come before an Impropers
   section.

*Must read Atoms before Lines*
   The Atoms section of a data file must come before a Lines section.

*Must read Atoms before Triangles*
   The Atoms section of a data file must come before a Triangles section.

*Must read Atoms before Velocities*
   The Atoms section of a data file must come before a Velocities
   section.

*Must re-specify non-restarted pair style (xxx) after read_restart*
   For pair styles, that do not store their settings in a restart file,
   it must be defined with a new 'pair_style' command after read_restart.

*Must set both respa inner and outer*
   Cannot use just the inner or outer option with respa without using the
   other.

*Must set number of threads via package omp command*
   Because you are using the USER-OMP package, set the number of threads
   via its settings, not by the pair_style snap nthreads setting.

*Must shrink-wrap piston boundary*
   The boundary style of the face where the piston is applied must be of
   type s (shrink-wrapped).

*Must specify a region in fix deposit*
   The region keyword must be specified with this fix.

*Must specify a region in fix pour*
   Self-explanatory.

*Must specify at least 2 types in fix atom/swap command*
   Self-explanatory.

*Must use 'kim_style init' command before simulation box is defined*
   Self-explanatory.

*Must use 'kim_style define' command after simulation box is defined*
   Self-explanatory.

*Must use 'kim_style init' command before 'kim_style define'*
   Self-explanatory.

*Must use 'kspace_modify pressure/scalar no' for rRESPA with kspace_style MSM*
   The kspace scalar pressure option cannot (yet) be used with rRESPA.

*Must use 'kspace_modify pressure/scalar no' for tensor components with kspace_style msm*
   Otherwise MSM will compute only a scalar pressure.  See the kspace_modify
   command for details on this setting.

*Must use 'kspace_modify pressure/scalar no' to obtain per-atom virial with kspace_style MSM*
   The kspace scalar pressure option cannot be used to obtain per-atom virial.

*Must use 'kspace_modify pressure/scalar no' with GPU MSM Pair styles*
   The kspace scalar pressure option is not (yet) compatible with GPU MSM Pair styles.

*Must use 'kspace_modify pressure/scalar no' with kspace_style msm/cg*
   The kspace scalar pressure option is not compatible with kspace_style msm/cg.

*Must use -in switch with multiple partitions*
   A multi-partition simulation cannot read the input script from stdin.
   The -in command-line option must be used to specify a file.

*Must use Kokkos half/thread or full neighbor list with threads or GPUs*
   Using Kokkos half-neighbor lists with threading is not allowed.

*Must use a block or cylinder region with fix pour*
   Self-explanatory.

*Must use a block region with fix pour for 2d simulations*
   Self-explanatory.

*Must use a bond style with TIP4P potential*
   TIP4P potentials assume bond lengths in water are constrained
   by a fix shake command.

*Must use a molecular atom style with fix poems molecule*
   Self-explanatory.

*Must use a z-axis cylinder region with fix pour*
   Self-explanatory.

*Must use an angle style with TIP4P potential*
   TIP4P potentials assume angles in water are constrained by a fix shake
   command.

*Must use atom map style array with Kokkos*
   See the atom_modify map command.

*Must use atom style with molecule IDs with fix bond/swap*
   Self-explanatory.

*Must use pair_style comb or comb3 with fix qeq/comb*
   Self-explanatory.

*Must use variable energy with fix addforce*
   Must define an energy variable when applying a dynamic
   force during minimization.

*Must use variable energy with fix efield*
   You must define an energy when performing a minimization with a
   variable E-field.

*NEB command before simulation box is defined*
   Self-explanatory.

*NEB requires damped dynamics minimizer*
   Use a different minimization style.

*NEB requires use of fix neb*
   Self-explanatory.

*NL ramp in wall/piston only implemented in zlo for now*
   The ramp keyword can only be used for piston applied to face zlo.

*Need nswaptypes mu values in fix atom/swap command*
   Self-explanatory.

*Needed bonus data not in data file*
   Some atom styles require bonus data.  See the read_data doc page for
   details.

*Needed molecular topology not in data file*
   The header of the data file indicated bonds, angles, etc would be
   included, but they are not present.

*Neigh_modify exclude molecule requires atom attribute molecule*
   Self-explanatory.

*Neigh_modify include group != atom_modify first group*
   Self-explanatory.

*Neighbor delay must be 0 or multiple of every setting*
   The delay and every parameters set via the neigh_modify command are
   inconsistent.  If the delay setting is non-zero, then it must be a
   multiple of the every setting.

*Neighbor include group not allowed with ghost neighbors*
   This is a current restriction within LAMMPS.

*Neighbor list overflow, boost neigh_modify one*
   There are too many neighbors of a single atom.  Use the neigh_modify
   command to increase the max number of neighbors allowed for one atom.
   You may also want to boost the page size.

*Neighbor multi not yet enabled for ghost neighbors*
   This is a current restriction within LAMMPS.

*Neighbor multi not yet enabled for granular*
   Self-explanatory.

*Neighbor multi not yet enabled for rRESPA*
   Self-explanatory.

*Neighbor page size must be >= 10x the one atom setting*
   This is required to prevent wasting too much memory.

*New atom IDs exceed maximum allowed ID*
   See the setting for tagint in the src/lmptype.h file.

*New bond exceeded bonds per atom in create_bonds*
See the read_data command for info on using the "extra/bond/per/atom"
keyword to allow for additional bonds to be formed

*New bond exceeded bonds per atom in fix bond/create*
   See the read_data command for info on using the "extra/bond/per/atom"
   keyword to allow for additional bonds to be formed

*New bond exceeded special list size in fix bond/create*
   See the "read_data extra/special/per/atom" command
   (or the "create_box extra/special/per/atom" command)
   for info on how to leave space in the special bonds
   list to allow for additional bonds to be formed.

*Newton bond change after simulation box is defined*
   The newton command cannot be used to change the newton bond value
   after a read_data, read_restart, or create_box command.

*Next command must list all universe and uloop variables*
   This is to insure they stay in sync.

*No Kspace style defined for compute group/group*
   Self-explanatory.

*No OpenMP support compiled in*
   An OpenMP flag is set, but LAMMPS was not built with
   OpenMP support.

*No angle style is defined for compute angle/local*
   Self-explanatory.

*No angles allowed with this atom style*
   Self-explanatory.

*No atoms in data file*
   The header of the data file indicated that atoms would be included,
   but they are not present.

*No basis atoms in lattice*
   Basis atoms must be defined for lattice style user.

*No bodies allowed with this atom style*
   Self-explanatory.  Check data file.

*No bond style is defined for compute bond/local*
   Self-explanatory.

*No bonds allowed with this atom style*
   Self-explanatory.

*No box information in dump. You have to use 'box no'*
   Self-explanatory.

*No count or invalid atom count in molecule file*
   The number of atoms must be specified.

*No dihedral style is defined for compute dihedral/local*
   Self-explanatory.

*No dihedrals allowed with this atom style*
   Self-explanatory.

*No dump custom arguments specified*
   The dump custom command requires that atom quantities be specified to
   output to dump file.

*No dump local arguments specified*
   Self-explanatory.

*No ellipsoids allowed with this atom style*
   Self-explanatory.  Check data file.

*No fix gravity defined for fix pour*
   Gravity is required to use fix pour.

*No improper style is defined for compute improper/local*
   Self-explanatory.

*No impropers allowed with this atom style*
   Self-explanatory.

*No input values for fix ave/spatial*
   Self-explanatory.

*No lines allowed with this atom style*
   Self-explanatory.  Check data file.

*No matching element in ADP potential file*
   The ADP potential file does not contain elements that match the
   requested elements.

*No matching element in EAM potential file*
   The EAM potential file does not contain elements that match the
   requested elements.

*No molecule topology allowed with atom style template*
   The data file cannot specify the number of bonds, angles, etc,
   because this info if inferred from the molecule templates.

*No overlap of box and region for create_atoms*
   Self-explanatory.

*No pair coul/streitz for fix qeq/slater*
   These commands must be used together.

*No pair hbond/dreiding coefficients set*
   Self-explanatory.

*No pair style defined for compute group/group*
   Cannot calculate group interactions without a pair style defined.

*No pair style is defined for compute pair/local*
   Self-explanatory.

*No pair style is defined for compute property/local*
   Self-explanatory.

*No rigid bodies defined*
   The fix specification did not end up defining any rigid bodies.

*No triangles allowed with this atom style*
   Self-explanatory.  Check data file.

*No values in fix ave/chunk command*
   Self-explanatory.

*No values in fix ave/time command*
   Self-explanatory.

*Non digit character between brackets in variable*
   Self-explanatory.

*Non integer # of swaps in temper command*
   Swap frequency in temper command must evenly divide the total # of
   timesteps.

*Non-numeric box dimensions - simulation unstable*
   The box size has apparently blown up.

*Non-zero atom IDs with atom_modify id = no*
   Self-explanatory.

*Non-zero read_data shift z value for 2d simulation*
   Self-explanatory.

*Nprocs not a multiple of N for -reorder*
   Self-explanatory.

*Number of core atoms != number of shell atoms*
   There must be a one-to-one pairing of core and shell atoms.

*Numeric index is out of bounds*
   A command with an argument that specifies an integer or range of
   integers is using a value that is less than 1 or greater than the
   maximum allowed limit.

*One or more Atom IDs is negative*
   Atom IDs must be positive integers.

*One or more atom IDs is too big*
   The limit on atom IDs is set by the SMALLBIG, BIGBIG, SMALLSMALL
   setting in your LAMMPS build.  See the :doc:`Build settings <Build_settings>` doc page for more info.

*One or more atom IDs is zero*
   Either all atoms IDs must be zero or none of them.

*One or more atoms belong to multiple rigid bodies*
   Two or more rigid bodies defined by the fix rigid command cannot
   contain the same atom.

*One or more rigid bodies are a single particle*
   Self-explanatory.

*One or zero atoms in rigid body*
   Any rigid body defined by the fix rigid command must contain 2 or more
   atoms.

*Only 2 types allowed when not using semi-grand in fix atom/swap command*
   Self-explanatory.

*Only one cut-off allowed when requesting all long*
   Self-explanatory.

*Only one cutoff allowed when requesting all long*
   Self-explanatory.

*Only zhi currently implemented for fix append/atoms*
   Self-explanatory.

*Out of range atoms - cannot compute MSM*
   One or more atoms are attempting to map their charge to a MSM grid point
   that is not owned by a processor.  This is likely for one of two
   reasons, both of them bad.  First, it may mean that an atom near the
   boundary of a processor's sub-domain has moved more than 1/2 the
   :doc:`neighbor skin distance <neighbor>` without neighbor lists being
   rebuilt and atoms being migrated to new processors.  This also means
   you may be missing pairwise interactions that need to be computed.
   The solution is to change the re-neighboring criteria via the
   :doc:`neigh_modify <neigh_modify>` command.  The safest settings are
   "delay 0 every 1 check yes".  Second, it may mean that an atom has
   moved far outside a processor's sub-domain or even the entire
   simulation box. This indicates bad physics, e.g. due to highly
   overlapping atoms, too large a timestep, etc.

*Out of range atoms - cannot compute PPPM*
   One or more atoms are attempting to map their charge to a PPPM grid
   point that is not owned by a processor.  This is likely for one of two
   reasons, both of them bad.  First, it may mean that an atom near the
   boundary of a processor's sub-domain has moved more than 1/2 the
   :doc:`neighbor skin distance <neighbor>` without neighbor lists being
   rebuilt and atoms being migrated to new processors.  This also means
   you may be missing pairwise interactions that need to be computed.
   The solution is to change the re-neighboring criteria via the
   :doc:`neigh_modify <neigh_modify>` command.  The safest settings are
   "delay 0 every 1 check yes".  Second, it may mean that an atom has
   moved far outside a processor's sub-domain or even the entire
   simulation box. This indicates bad physics, e.g. due to highly
   overlapping atoms, too large a timestep, etc.

*Out of range atoms - cannot compute PPPMDisp*
   One or more atoms are attempting to map their charge to a PPPM grid
   point that is not owned by a processor.  This is likely for one of two
   reasons, both of them bad.  First, it may mean that an atom near the
   boundary of a processor's sub-domain has moved more than 1/2 the
   :doc:`neighbor skin distance <neighbor>` without neighbor lists being
   rebuilt and atoms being migrated to new processors.  This also means
   you may be missing pairwise interactions that need to be computed.
   The solution is to change the re-neighboring criteria via the
   :doc:`neigh_modify <neigh_modify>` command.  The safest settings are
   "delay 0 every 1 check yes".  Second, it may mean that an atom has
   moved far outside a processor's sub-domain or even the entire
   simulation box. This indicates bad physics, e.g. due to highly
   overlapping atoms, too large a timestep, etc.

*Overflow of allocated fix vector storage*
   This should not normally happen if the fix correctly calculated
   how long the vector will grow to.  Contact the developers.

*Overlapping large/large in pair colloid*
   This potential is infinite when there is an overlap.

*Overlapping small/large in pair colloid*
   This potential is infinite when there is an overlap.

*POEMS fix must come before NPT/NPH fix*
   NPT/NPH fix must be defined in input script after all poems fixes,
   else the fix contribution to the pressure virial is incorrect.

*PPPM can only currently be used with comm_style brick*
   This is a current restriction in LAMMPS.

*PPPM grid is too large*
   The global PPPM grid is larger than OFFSET in one or more dimensions.
   OFFSET is currently set to 4096.  You likely need to decrease the
   requested accuracy.

*PPPM grid stencil extends beyond nearest neighbor processor*
   This is not allowed if the kspace_modify overlap setting is no.

*PPPM order < minimum allowed order*
   The default minimum order is 2.  This can be reset by the
   kspace_modify minorder command.

*PPPM order cannot be < 2 or > than %d*
   This is a limitation of the PPPM implementation in LAMMPS.

*PPPMDisp Coulomb grid is too large*
   The global PPPM grid is larger than OFFSET in one or more dimensions.
   OFFSET is currently set to 4096.  You likely need to decrease the
   requested accuracy.

*PPPMDisp Dispersion grid is too large*
   The global PPPM grid is larger than OFFSET in one or more dimensions.
   OFFSET is currently set to 4096.  You likely need to decrease the
   requested accuracy.

*PPPMDisp can only currently be used with comm_style brick*
   This is a current restriction in LAMMPS.

*PPPMDisp coulomb order cannot be greater than %d*
   This is a limitation of the PPPM implementation in LAMMPS.

*PPPMDisp used but no parameters set, for further information please see the pppm/disp documentation*
   An efficient and accurate usage of the pppm/disp requires settings via the kspace_modify command. Please see the pppm/disp documentation for further instructions.

*PRD command before simulation box is defined*
   The prd command cannot be used before a read_data,
   read_restart, or create_box command.

*PRD nsteps must be multiple of t_event*
   Self-explanatory.

*PRD t_corr must be multiple of t_event*
   Self-explanatory.

*Package command after simulation box is defined*
   The package command cannot be used after a read_data, read_restart, or
   create_box command.

*Package gpu command without GPU package installed*
   The GPU package must be installed via "make yes-gpu" before LAMMPS is
   built.

*Package intel command without USER-INTEL package installed*
   The USER-INTEL package must be installed via "make yes-user-intel"
   before LAMMPS is built.

*Package kokkos command without KOKKOS package enabled*
   The KOKKOS package must be installed via "make yes-kokkos" before
   LAMMPS is built, and the "-k on" must be used to enable the package.

*Package omp command without USER-OMP package installed*
   The USER-OMP package must be installed via "make yes-user-omp" before
   LAMMPS is built.

*Pair body requires atom style body*
   Self-explanatory.

*Pair body requires body style nparticle*
   This pair style is specific to the nparticle body style.

*Pair brownian requires atom style sphere*
   Self-explanatory.

*Pair brownian requires extended particles*
   One of the particles has radius 0.0.

*Pair brownian requires monodisperse particles*
   All particles must be the same finite size.

*Pair brownian/poly requires atom style sphere*
   Self-explanatory.

*Pair brownian/poly requires extended particles*
   One of the particles has radius 0.0.

*Pair brownian/poly requires newton pair off*
   Self-explanatory.

*Pair coeff for hybrid has invalid style*
   Style in pair coeff must have been listed in pair_style command.

*Pair coul/wolf requires atom attribute q*
   The atom style defined does not have this attribute.

*Pair cutoff < Respa interior cutoff*
   One or more pairwise cutoffs are too short to use with the specified
   rRESPA cutoffs.

*Pair dipole/cut requires atom attributes q, mu, torque*
   The atom style defined does not have these attributes.

*Pair dipole/cut/gpu requires atom attributes q, mu, torque*
   The atom style defined does not have this attribute.

*Pair dipole/long requires atom attributes q, mu, torque*
   The atom style defined does not have these attributes.

*Pair dipole/sf/gpu requires atom attributes q, mu, torque*
   The atom style defined does not one or more of these attributes.

*Pair distance < table inner cutoff*
   Two atoms are closer together than the pairwise table allows.

*Pair distance > table outer cutoff*
   Two atoms are further apart than the pairwise table allows.

*Pair dpd requires ghost atoms store velocity*
   Use the comm_modify vel yes command to enable this.

*Pair gayberne epsilon a,b,c coeffs are not all set*
   Each atom type involved in pair_style gayberne must
   have these 3 coefficients set at least once.

*Pair gayberne requires atom style ellipsoid*
   Self-explanatory.

*Pair gayberne requires atoms with same type have same shape*
   Self-explanatory.

*Pair gayberne/gpu requires atom style ellipsoid*
   Self-explanatory.

*Pair gayberne/gpu requires atoms with same type have same shape*
   Self-explanatory.

*Pair granular requires atom attributes radius, rmass*
   The atom style defined does not have these attributes.

*Pair granular requires ghost atoms store velocity*
   Use the comm_modify vel yes command to enable this.

*Pair granular with shear history requires newton pair off*
   This is a current restriction of the implementation of pair
   granular styles with history.

*Pair hybrid single calls do not support per sub-style special bond values*
   Self-explanatory.

*Pair hybrid sub-style does not support single call*
   You are attempting to invoke a single() call on a pair style
   that doesn't support it.

*Pair hybrid sub-style is not used*
   No pair_coeff command used a sub-style specified in the pair_style
   command.

*Pair inner cutoff < Respa interior cutoff*
   One or more pairwise cutoffs are too short to use with the specified
   rRESPA cutoffs.

*Pair inner cutoff >= Pair outer cutoff*
   The specified cutoffs for the pair style are inconsistent.

*Pair line/lj requires atom style line*
   Self-explanatory.

*Pair lj/long/dipole/long requires atom attributes mu, torque*
   The atom style defined does not have these attributes.

*Pair lubricate requires atom style sphere*
   Self-explanatory.

*Pair lubricate requires ghost atoms store velocity*
   Use the comm_modify vel yes command to enable this.

*Pair lubricate requires monodisperse particles*
   All particles must be the same finite size.

*Pair lubricate/poly requires atom style sphere*
   Self-explanatory.

*Pair lubricate/poly requires extended particles*
   One of the particles has radius 0.0.

*Pair lubricate/poly requires ghost atoms store velocity*
   Use the comm_modify vel yes command to enable this.

*Pair lubricate/poly requires newton pair off*
   Self-explanatory.

*Pair lubricateU requires atom style sphere*
   Self-explanatory.

*Pair lubricateU requires ghost atoms store velocity*
   Use the comm_modify vel yes command to enable this.

*Pair lubricateU requires monodisperse particles*
   All particles must be the same finite size.

*Pair lubricateU/poly requires ghost atoms store velocity*
   Use the comm_modify vel yes command to enable this.

*Pair lubricateU/poly requires newton pair off*
   Self-explanatory.

*Pair peri lattice is not identical in x, y, and z*
   The lattice defined by the lattice command must be cubic.

*Pair peri requires a lattice be defined*
   Use the lattice command for this purpose.

*Pair peri requires an atom map, see atom_modify*
   Even for atomic systems, an atom map is required to find Peridynamic
   bonds.  Use the atom_modify command to define one.

*Pair resquared epsilon a,b,c coeffs are not all set*
   Self-explanatory.

*Pair resquared epsilon and sigma coeffs are not all set*
   Self-explanatory.

*Pair resquared requires atom style ellipsoid*
   Self-explanatory.

*Pair resquared requires atoms with same type have same shape*
   Self-explanatory.

*Pair resquared/gpu requires atom style ellipsoid*
   Self-explanatory.

*Pair resquared/gpu requires atoms with same type have same shape*
   Self-explanatory.

*Pair style AIREBO requires atom IDs*
   This is a requirement to use the AIREBO potential.

*Pair style AIREBO requires newton pair on*
   See the newton command.  This is a restriction to use the AIREBO
   potential.

*Pair style BOP requires atom IDs*
   This is a requirement to use the BOP potential.

*Pair style BOP requires newton pair on*
   See the newton command.  This is a restriction to use the BOP
   potential.

*Pair style COMB requires atom IDs*
   This is a requirement to use the AIREBO potential.

*Pair style COMB requires atom attribute q*
   Self-explanatory.

*Pair style COMB requires newton pair on*
   See the newton command.  This is a restriction to use the COMB
   potential.

*Pair style COMB3 requires atom IDs*
   This is a requirement to use the COMB3 potential.

*Pair style COMB3 requires atom attribute q*
   Self-explanatory.

*Pair style COMB3 requires newton pair on*
   See the newton command.  This is a restriction to use the COMB3
   potential.

*Pair style LCBOP requires atom IDs*
   This is a requirement to use the LCBOP potential.

*Pair style LCBOP requires newton pair on*
   See the newton command.  This is a restriction to use the Tersoff
   potential.

*Pair style MEAM requires newton pair on*
   See the newton command.  This is a restriction to use the MEAM
   potential.

*Pair style SNAP requires newton pair on*
   See the newton command.  This is a restriction to use the SNAP
   potential.

*Pair style Stillinger-Weber requires atom IDs*
   This is a requirement to use the SW potential.

*Pair style Stillinger-Weber requires newton pair on*
   See the newton command.  This is a restriction to use the SW
   potential.

*Pair style Tersoff requires atom IDs*
   This is a requirement to use the Tersoff potential.

*Pair style Tersoff requires newton pair on*
   See the newton command.  This is a restriction to use the Tersoff
   potential.

*Pair style Vashishta requires atom IDs*
   This is a requirement to use the Vashishta potential.

*Pair style Vashishta requires newton pair on*
   See the newton command.  This is a restriction to use the Vashishta
   potential.

*Pair style bop requires comm ghost cutoff at least 3x larger than %g*
   Use the communicate ghost command to set this.  See the pair bop
   doc page for more details.

*Pair style born/coul/long requires atom attribute q*
   An atom style that defines this attribute must be used.

*Pair style born/coul/long/gpu requires atom attribute q*
   The atom style defined does not have this attribute.

*Pair style born/coul/wolf requires atom attribute q*
   The atom style defined does not have this attribute.

*Pair style buck/coul/cut requires atom attribute q*
   The atom style defined does not have this attribute.

*Pair style buck/coul/long requires atom attribute q*
   The atom style defined does not have these attributes.

*Pair style buck/coul/long/gpu requires atom attribute q*
   The atom style defined does not have this attribute.

*Pair style buck/long/coul/long requires atom attribute q*
   The atom style defined does not have this attribute.

*Pair style coul/cut requires atom attribute q*
   The atom style defined does not have these attributes.

*Pair style coul/cut/gpu requires atom attribute q*
   The atom style defined does not have this attribute.

*Pair style coul/debye/gpu requires atom attribute q*
   The atom style defined does not have this attribute.

*Pair style coul/dsf requires atom attribute q*
   The atom style defined does not have this attribute.

*Pair style coul/dsf/gpu requires atom attribute q*
   The atom style defined does not have this attribute.

*Pair style coul/long/gpu requires atom attribute q*
   The atom style defined does not have these attributes.

*Pair style coul/streitz requires atom attribute q*
   Self-explanatory.

*Pair style does not have extra field requested by compute pair/local*
   The pair style does not support the pN value requested by the compute
   pair/local command.

*Pair style does not support bond_style quartic*
   The pair style does not have a single() function, so it can
   not be invoked by bond_style quartic.

*Pair style does not support compute group/group*
   The pair_style does not have a single() function, so it cannot be
   invoked by the compute group/group command.

*Pair style does not support compute pair/local*
   The pair style does not have a single() function, so it can
   not be invoked by compute pair/local.

*Pair style does not support compute property/local*
   The pair style does not have a single() function, so it can
   not be invoked by fix bond/swap.

*Pair style does not support fix bond/swap*
   The pair style does not have a single() function, so it can
   not be invoked by fix bond/swap.

*Pair style does not support pair_write*
   The pair style does not have a single() function, so it can
   not be invoked by pair write.

*Pair style does not support rRESPA inner/middle/outer*
   You are attempting to use rRESPA options with a pair style that
   does not support them.

*Pair style granular with history requires atoms have IDs*
   Atoms in the simulation do not have IDs, so history effects
   cannot be tracked by the granular pair potential.

*Pair style hbond/dreiding requires an atom map, see atom_modify*
   Self-explanatory.

*Pair style hbond/dreiding requires atom IDs*
   Self-explanatory.

*Pair style hbond/dreiding requires molecular system*
   Self-explanatory.

*Pair style hbond/dreiding requires newton pair on*
   See the newton command for details.

*Pair style hybrid cannot have hybrid as an argument*
   Self-explanatory.

*Pair style hybrid cannot have none as an argument*
   Self-explanatory.

*Pair style is incompatible with KSpace style*
   If a pair style with a long-range Coulombic component is selected,
   then a kspace style must also be used.

*Pair style is incompatible with TIP4P KSpace style*
   The pair style does not have the requires TIP4P settings.

*Pair style lj/charmm/coul/charmm requires atom attribute q*
   The atom style defined does not have these attributes.

*Pair style lj/charmm/coul/long requires atom attribute q*
   The atom style defined does not have these attributes.

*Pair style lj/charmm/coul/long/gpu requires atom attribute q*
   The atom style defined does not have this attribute.

*Pair style lj/class2/coul/cut requires atom attribute q*
   The atom style defined does not have this attribute.

*Pair style lj/class2/coul/long requires atom attribute q*
   The atom style defined does not have this attribute.

*Pair style lj/class2/coul/long/gpu requires atom attribute q*
   The atom style defined does not have this attribute.

*Pair style lj/cut/coul/cut requires atom attribute q*
   The atom style defined does not have this attribute.

*Pair style lj/cut/coul/cut/gpu requires atom attribute q*
   The atom style defined does not have this attribute.

*Pair style lj/cut/coul/debye/gpu requires atom attribute q*
   The atom style defined does not have this attribute.

*Pair style lj/cut/coul/dsf requires atom attribute q*
   The atom style defined does not have these attributes.

*Pair style lj/cut/coul/dsf/gpu requires atom attribute q*
   The atom style defined does not have this attribute.

*Pair style lj/cut/coul/long requires atom attribute q*
   The atom style defined does not have this attribute.

*Pair style lj/cut/coul/long/gpu requires atom attribute q*
   The atom style defined does not have this attribute.

*Pair style lj/cut/tip4p/cut requires atom IDs*
   This is a requirement to use this potential.

*Pair style lj/cut/tip4p/cut requires atom attribute q*
   The atom style defined does not have this attribute.

*Pair style lj/cut/tip4p/cut requires newton pair on*
   See the newton command.  This is a restriction to use this
   potential.

*Pair style lj/cut/tip4p/long requires atom IDs*
   There are no atom IDs defined in the system and the TIP4P potential
   requires them to find O,H atoms with a water molecule.

*Pair style lj/cut/tip4p/long requires atom attribute q*
   The atom style defined does not have these attributes.

*Pair style lj/cut/tip4p/long requires newton pair on*
   This is because the computation of constraint forces within a water
   molecule adds forces to atoms owned by other processors.

*Pair style lj/gromacs/coul/gromacs requires atom attribute q*
   An atom_style with this attribute is needed.

*Pair style lj/long/dipole/long does not currently support respa*
   This feature is not yet supported.

*Pair style lj/long/tip4p/long requires atom IDs*
   There are no atom IDs defined in the system and the TIP4P potential
   requires them to find O,H atoms with a water molecule.

*Pair style lj/long/tip4p/long requires atom attribute q*
   The atom style defined does not have these attributes.

*Pair style lj/long/tip4p/long requires newton pair on*
   This is because the computation of constraint forces within a water
   molecule adds forces to atoms owned by other processors.

*Pair style lj/sdk/coul/long/gpu requires atom attribute q*
   The atom style defined does not have this attribute.

*Pair style nb3b/harmonic requires atom IDs*
   This is a requirement to use this potential.

*Pair style nb3b/harmonic requires newton pair on*
   See the newton command.  This is a restriction to use this potential.

*Pair style nm/cut/coul/cut requires atom attribute q*
   The atom style defined does not have this attribute.

*Pair style nm/cut/coul/long requires atom attribute q*
   The atom style defined does not have this attribute.

*Pair style peri requires atom style peri*
   Self-explanatory.

*Pair style polymorphic requires atom IDs*
   This is a requirement to use the polymorphic potential.

*Pair style polymorphic requires newton pair on*
   See the newton command.  This is a restriction to use the polymorphic
   potential.

*Pair style reax requires atom IDs*
   This is a requirement to use the ReaxFF potential.

*Pair style reax requires atom attribute q*
   The atom style defined does not have this attribute.

*Pair style reax requires newton pair on*
   This is a requirement to use the ReaxFF potential.

*Pair style requires a KSpace style*
   No kspace style is defined.

*Pair style requires use of kspace_style ewald/disp*
   Self-explanatory.

*Pair style sw/gpu requires atom IDs*
   This is a requirement to use this potential.

*Pair style sw/gpu requires newton pair off*
   See the newton command.  This is a restriction to use this potential.

*Pair style vashishta/gpu requires atom IDs*
   This is a requirement to use this potential.

*Pair style vashishta/gpu requires newton pair off*
   See the newton command.  This is a restriction to use this potential.

*Pair style tersoff/gpu requires atom IDs*
   This is a requirement to use the tersoff/gpu potential.

*Pair style tersoff/gpu requires newton pair off*
   See the newton command.  This is a restriction to use this pair style.

*Pair style tip4p/cut requires atom IDs*
   This is a requirement to use this potential.

*Pair style tip4p/cut requires atom attribute q*
   The atom style defined does not have this attribute.

*Pair style tip4p/cut requires newton pair on*
   See the newton command.  This is a restriction to use this potential.

*Pair style tip4p/long requires atom IDs*
   There are no atom IDs defined in the system and the TIP4P potential
   requires them to find O,H atoms with a water molecule.

*Pair style tip4p/long requires atom attribute q*
   The atom style defined does not have these attributes.

*Pair style tip4p/long requires newton pair on*
   This is because the computation of constraint forces within a water
   molecule adds forces to atoms owned by other processors.

*Pair table cutoffs must all be equal to use with KSpace*
   When using pair style table with a long-range KSpace solver, the
   cutoffs for all atom type pairs must all be the same, since the
   long-range solver starts at that cutoff.

*Pair table parameters did not set N*
   List of pair table parameters must include N setting.

*Pair tersoff/zbl requires metal or real units*
   This is a current restriction of this pair potential.

*Pair tersoff/zbl/kk requires metal or real units*
   This is a current restriction of this pair potential.

*Pair tri/lj requires atom style tri*
   Self-explanatory.

*Pair yukawa/colloid requires atom style sphere*
   Self-explanatory.

*Pair yukawa/colloid requires atoms with same type have same radius*
   Self-explanatory.

*Pair yukawa/colloid/gpu requires atom style sphere*
   Self-explanatory.

*PairKIM only works with 3D problems*
   This is a current limitation.

*Pair_coeff command before pair_style is defined*
   Self-explanatory.

*Pair_coeff command before simulation box is defined*
   The pair_coeff command cannot be used before a read_data,
   read_restart, or create_box command.

*Pair_modify command before pair_style is defined*
   Self-explanatory.

*Pair_modify special setting for pair hybrid incompatible with global special_bonds setting*
   Cannot override a setting of 0.0 or 1.0 or change a setting between
   0.0 and 1.0.

*Pair_write command before pair_style is defined*
   Self-explanatory.

*Particle on or inside fix wall surface*
   Particles must be "exterior" to the wall in order for energy/force to
   be calculated.

*Particle outside surface of region used in fix wall/region*
   Particles must be inside the region for energy/force to be calculated.
   A particle outside the region generates an error.

*Per-atom compute in equal-style variable formula*
   Equal-style variables cannot use per-atom quantities.

*Per-atom energy was not tallied on needed timestep*
   You are using a thermo keyword that requires potentials to
   have tallied energy, but they didn't on this timestep.  See the
   variable doc page for ideas on how to make this work.

*Per-atom fix in equal-style variable formula*
   Equal-style variables cannot use per-atom quantities.

*Per-atom virial was not tallied on needed timestep*
   You are using a thermo keyword that requires potentials to have
   tallied the virial, but they didn't on this timestep.  See the
   variable doc page for ideas on how to make this work.

*Per-processor system is too big*
   The number of owned atoms plus ghost atoms on a single
   processor must fit in 32-bit integer.

*Potential energy ID for fix neb does not exist*
   Self-explanatory.

*Potential energy ID for fix nvt/nph/npt does not exist*
   A compute for potential energy must be defined.

*Potential file has duplicate entry*
   The potential file has more than one entry for the same element.

*Potential file is missing an entry*
   The potential file does not have a needed entry.

*Power by 0 in variable formula*
   Self-explanatory.

*Pressure ID for fix box/relax does not exist*
   The compute ID needed to compute pressure for the fix does not
   exist.

*Pressure ID for fix modify does not exist*
   Self-explanatory.

*Pressure ID for fix npt/nph does not exist*
   Self-explanatory.

*Pressure ID for fix press/berendsen does not exist*
   The compute ID needed to compute pressure for the fix does not
   exist.

*Pressure ID for fix rigid npt/nph does not exist*
   Self-explanatory.

*Pressure ID for thermo does not exist*
   The compute ID needed to compute pressure for thermodynamics does not
   exist.

*Pressure control can not be used with fix nvt*
   Self-explanatory.

*Pressure control can not be used with fix nvt/asphere*
   Self-explanatory.

*Pressure control can not be used with fix nvt/body*
   Self-explanatory.

*Pressure control can not be used with fix nvt/sllod*
   Self-explanatory.

*Pressure control can not be used with fix nvt/sphere*
   Self-explanatory.

*Pressure control must be used with fix nph*
   Self-explanatory.

*Pressure control must be used with fix nph/asphere*
   Self-explanatory.

*Pressure control must be used with fix nph/body*
   Self-explanatory.

*Pressure control must be used with fix nph/small*
   Self-explanatory.

*Pressure control must be used with fix nph/sphere*
   Self-explanatory.

*Pressure control must be used with fix nphug*
   A pressure control keyword (iso, aniso, tri, x, y, or z) must be
   provided.

*Pressure control must be used with fix npt*
   Self-explanatory.

*Pressure control must be used with fix npt/asphere*
   Self-explanatory.

*Pressure control must be used with fix npt/body*
   Self-explanatory.

*Pressure control must be used with fix npt/sphere*
   Self-explanatory.

*Processor count in z must be 1 for 2d simulation*
   Self-explanatory.

*Processor partitions do not match number of allocated processors*
   The total number of processors in all partitions must match the number
   of processors LAMMPS is running on.

*Processors command after simulation box is defined*
   The processors command cannot be used after a read_data, read_restart,
   or create_box command.

*Processors custom grid file is inconsistent*
   The vales in the custom file are not consistent with the number of
   processors you are running on or the Px,Py,Pz settings of the
   processors command.  Or there was not a setting for every processor.

*Processors grid numa and map style are incompatible*
   Using numa for gstyle in the processors command requires using
   cart for the map option.

*Processors part option and grid style are incompatible*
   Cannot use gstyle numa or custom with the part option.

*Processors twogrid requires proc count be a multiple of core count*
   Self-explanatory.

*Pstart and Pstop must have the same value*
   Self-explanatory.

*Python function evaluation failed*
   The Python function did not run successfully and/or did not return a
   value (if it is supposed to return a value).  This is probably due to
   some error condition in the function.

*Python function is not callable*
   The provided Python code was run successfully, but it not
   define a callable function with the required name.

*Python invoke of undefined function*
   Cannot invoke a function that has not been previously defined.

*Python variable does not match Python function*
   This matching is defined by the python-style variable and the python
   command.

*Python variable has no function*
   No python command was used to define the function associated with the
   python-style variable.

*QEQ with 'newton pair off' not supported*
   See the newton command.  This is a restriction to use the QEQ fixes.

*R0 < 0 for fix spring command*
   Equilibrium spring length is invalid.

*RATTLE coordinate constraints are not satisfied up to desired tolerance*
   Self-explanatory.

*RATTLE determinant = 0.0*
   The determinant of the matrix being solved for a single cluster
   specified by the fix rattle command is numerically invalid.

*RATTLE failed*
   Certain constraints were not satisfied.

*RATTLE velocity constraints are not satisfied up to desired tolerance*
   Self-explanatory.

*Read data add offset is too big*
   It cannot be larger than the size of atom IDs, e.g. the maximum 32-bit
   integer.

*Read dump of atom property that isn't allocated*
   Self-explanatory.

*Read rerun dump file timestep > specified stop*
   Self-explanatory.

*Read restart MPI-IO input not allowed with % in filename*
   This is because a % signifies one file per processor and MPI-IO
   creates one large file for all processors.

*Read_data shrink wrap did not assign all atoms correctly*
   This is typically because the box-size specified in the data file is
   large compared to the actual extent of atoms in a shrink-wrapped
   dimension.  When LAMMPS shrink-wraps the box atoms will be lost if the
   processor they are re-assigned to is too far away.  Choose a box
   size closer to the actual extent of the atoms.

*Read_dump command before simulation box is defined*
   The read_dump command cannot be used before a read_data, read_restart,
   or create_box command.

*Read_dump field not found in dump file*
   Self-explanatory.

*Read_dump triclinic status does not match simulation*
   Both the dump snapshot and the current LAMMPS simulation must
   be using either an orthogonal or triclinic box.

*Read_dump xyz fields do not have consistent scaling/wrapping*
   Self-explanatory.

*Reading from MPI-IO filename when MPIIO package is not installed*
   Self-explanatory.

*Reax_defs.h setting for NATDEF is too small*
   Edit the setting in the ReaxFF library and re-compile the
   library and re-build LAMMPS.

*Reax_defs.h setting for NNEIGHMAXDEF is too small*
   Edit the setting in the ReaxFF library and re-compile the
   library and re-build LAMMPS.

*Receiving partition in processors part command is already a receiver*
   Cannot specify a partition to be a receiver twice.

*Region ID for compute chunk/atom does not exist*
   Self-explanatory.

*Region ID for compute reduce/region does not exist*
   Self-explanatory.

*Region ID for compute temp/region does not exist*
   Self-explanatory.

*Region ID for dump custom does not exist*
   Self-explanatory.

*Region ID for fix addforce does not exist*
   Self-explanatory.

*Region ID for fix atom/swap does not exist*
   Self-explanatory.

*Region ID for fix ave/spatial does not exist*
   Self-explanatory.

*Region ID for fix aveforce does not exist*
   Self-explanatory.

*Region ID for fix deposit does not exist*
   Self-explanatory.

*Region ID for fix efield does not exist*
   Self-explanatory.

*Region ID for fix evaporate does not exist*
   Self-explanatory.

*Region ID for fix gcmc does not exist*
   Self-explanatory.

*Region ID for fix heat does not exist*
   Self-explanatory.

*Region ID for fix setforce does not exist*
   Self-explanatory.

*Region ID for fix wall/region does not exist*
   Self-explanatory.

*Region ID for group dynamic does not exist*
   Self-explanatory.

*Region ID in variable formula does not exist*
   Self-explanatory.

*Region cannot have 0 length rotation vector*
   Self-explanatory.

*Region for fix oneway does not exist*
   Self-explanatory.

*Region intersect region ID does not exist*
   Self-explanatory.

*Region union or intersect cannot be dynamic*
   The sub-regions can be dynamic, but not the combined region.

*Region union region ID does not exist*
   One or more of the region IDs specified by the region union command
   does not exist.

*Replacing a fix, but new style != old style*
   A fix ID can be used a 2nd time, but only if the style matches the
   previous fix.  In this case it is assumed you with to reset a fix's
   parameters.  This error may mean you are mistakenly re-using a fix ID
   when you do not intend to.

*Replicate command before simulation box is defined*
   The replicate command cannot be used before a read_data, read_restart,
   or create_box command.

*Replicate did not assign all atoms correctly*
   Atoms replicated by the replicate command were not assigned correctly
   to processors.  This is likely due to some atom coordinates being
   outside a non-periodic simulation box.

*Replicated system atom IDs are too big*
   See the setting for tagint in the src/lmptype.h file.

*Replicated system is too big*
   See the setting for bigint in the src/lmptype.h file.

*Required border comm not yet implemented with Kokkos*
   There are various limitations in the communication options supported
   by Kokkos.

*Rerun command before simulation box is defined*
   The rerun command cannot be used before a read_data, read_restart, or
   create_box command.

*Rerun dump file does not contain requested snapshot*
   Self-explanatory.

*Resetting timestep size is not allowed with fix move*
   This is because fix move is moving atoms based on elapsed time.

*Respa inner cutoffs are invalid*
   The first cutoff must be <= the second cutoff.

*Respa levels must be >= 1*
   Self-explanatory.

*Respa middle cutoffs are invalid*
   The first cutoff must be <= the second cutoff.

*Restart file MPI-IO output not allowed with % in filename*
   This is because a % signifies one file per processor and MPI-IO
   creates one large file for all processors.

*Restart file byte ordering is not recognized*
   The file does not appear to be a LAMMPS restart file since it doesn't
   contain a recognized byte-ordering flag at the beginning.

*Restart file byte ordering is swapped*
   The file was written on a machine with different byte-ordering than
   the machine you are reading it on.  Convert it to a text data file
   instead, on the machine you wrote it on.

*Restart file incompatible with current version*
   This is probably because you are trying to read a file created with a
   version of LAMMPS that is too old compared to the current version.
   Use your older version of LAMMPS and convert the restart file
   to a data file.

*Restart file is a MPI-IO file*
   The file is inconsistent with the filename you specified for it.

*Restart file is a multi-proc file*
   The file is inconsistent with the filename you specified for it.

*Restart file is not a MPI-IO file*
   The file is inconsistent with the filename you specified for it.

*Restart file is not a multi-proc file*
   The file is inconsistent with the filename you specified for it.

*Restart variable returned a bad timestep*
   The variable must return a timestep greater than the current timestep.

*Restrain atoms %d %d %d %d missing on proc %d at step %ld*
   The 4 atoms in a restrain dihedral specified by the fix restrain
   command are not all accessible to a processor.  This probably means an
   atom has moved too far.

*Restrain atoms %d %d %d missing on proc %d at step %ld*
   The 3 atoms in a restrain angle specified by the fix restrain
   command are not all accessible to a processor.  This probably means an
   atom has moved too far.

*Restrain atoms %d %d missing on proc %d at step %ld*
   The 2 atoms in a restrain bond specified by the fix restrain
   command are not all accessible to a processor.  This probably means an
   atom has moved too far.

*Reuse of compute ID*
   A compute ID cannot be used twice.

*Reuse of dump ID*
   A dump ID cannot be used twice.

*Reuse of molecule template ID*
   The template IDs must be unique.

*Reuse of region ID*
   A region ID cannot be used twice.

*Rigid body atoms %d %d missing on proc %d at step %ld*
   This means that an atom cannot find the atom that owns the rigid body
   it is part of, or vice versa.  The solution is to use the communicate
   cutoff command to insure ghost atoms are acquired from far enough away
   to encompass the max distance printed when the fix rigid/small command
   was invoked.

*Rigid body has degenerate moment of inertia*
   Fix poems will only work with bodies (collections of atoms) that have
   non-zero principal moments of inertia.  This means they must be 3 or
   more non-collinear atoms, even with joint atoms removed.

*Rigid fix must come before NPT/NPH fix*
   NPT/NPH fix must be defined in input script after all rigid fixes,
   else the rigid fix contribution to the pressure virial is
   incorrect.

*Rmask function in equal-style variable formula*
   Rmask is per-atom operation.

*Run command before simulation box is defined*
   The run command cannot be used before a read_data, read_restart, or
   create_box command.

*Run command start value is after start of run*
   Self-explanatory.

*Run command stop value is before end of run*
   Self-explanatory.

*Run_style command before simulation box is defined*
   The run_style command cannot be used before a read_data,
   read_restart, or create_box command.

*SRD bin size for fix srd differs from user request*
   Fix SRD had to adjust the bin size to fit the simulation box.  See the
   cubic keyword if you want this message to be an error vs warning.

*SRD bins for fix srd are not cubic enough*
   The bin shape is not within tolerance of cubic.  See the cubic
   keyword if you want this message to be an error vs warning.

*SRD particle %d started inside big particle %d on step %ld bounce %d*
   See the inside keyword if you want this message to be an error vs
   warning.

*SRD particle %d started inside wall %d on step %ld bounce %d*
   See the inside keyword if you want this message to be an error vs
   warning.

*Same dimension twice in fix ave/spatial*
   Self-explanatory.

*Sending partition in processors part command is already a sender*
   Cannot specify a partition to be a sender twice.

*Set command before simulation box is defined*
   The set command cannot be used before a read_data, read_restart,
   or create_box command.

*Set command floating point vector does not exist*
   Self-explanatory.

*Set command integer vector does not exist*
   Self-explanatory.

*Set command with no atoms existing*
   No atoms are yet defined so the set command cannot be used.

*Set region ID does not exist*
   Region ID specified in set command does not exist.

*Shake angles have different bond types*
   All 3-atom angle-constrained SHAKE clusters specified by the fix shake
   command that are the same angle type, must also have the same bond
   types for the 2 bonds in the angle.

*Shake atoms %d %d %d %d missing on proc %d at step %ld*
   The 4 atoms in a single shake cluster specified by the fix shake
   command are not all accessible to a processor.  This probably means
   an atom has moved too far.

*Shake atoms %d %d %d missing on proc %d at step %ld*
   The 3 atoms in a single shake cluster specified by the fix shake
   command are not all accessible to a processor.  This probably means
   an atom has moved too far.

*Shake atoms %d %d missing on proc %d at step %ld*
   The 2 atoms in a single shake cluster specified by the fix shake
   command are not all accessible to a processor.  This probably means
   an atom has moved too far.

*Shake cluster of more than 4 atoms*
   A single cluster specified by the fix shake command can have no more
   than 4 atoms.

*Shake clusters are connected*
   A single cluster specified by the fix shake command must have a single
   central atom with up to 3 other atoms bonded to it.

*Shake determinant = 0.0*
   The determinant of the matrix being solved for a single cluster
   specified by the fix shake command is numerically invalid.

*Shake fix must come before NPT/NPH fix*
   NPT fix must be defined in input script after SHAKE fix, else the
   SHAKE fix contribution to the pressure virial is incorrect.

*Shear history overflow, boost neigh_modify one*
   There are too many neighbors of a single atom.  Use the neigh_modify
   command to increase the max number of neighbors allowed for one atom.
   You may also want to boost the page size.

*Small to big integers are not sized correctly*
   This error occurs when the sizes of smallint, imageint, tagint, bigint,
   as defined in src/lmptype.h are not what is expected.  Contact
   the developers if this occurs.

*Smallint setting in lmptype.h is invalid*
   It has to be the size of an integer.

*Smallint setting in lmptype.h is not compatible*
   Smallint stored in restart file is not consistent with LAMMPS version
   you are running.

*Special list size exceeded in fix bond/create*
   See the "read_data extra/special/per/atom" command
   (or the "create_box extra/special/per/atom" command)
   for info on how to leave space in the special bonds
   list to allow for additional bonds to be formed.

*Species XXX is not supported by this KIM Simulator Model*
   The kim_style define command was referencing a species that is not
   present in the requested KIM Simulator Model.

*Specified processors != physical processors*
   The 3d grid of processors defined by the processors command does not
   match the number of processors LAMMPS is being run on.

*Specified target stress must be uniaxial or hydrostatic*
   Self-explanatory.

*Sqrt of negative value in variable formula*
   Self-explanatory.

*Subsequent read data induced too many angles per atom*
   See the extra/angle/per/atom keyword for the create_box
   or the read_data command to set this limit larger

*Subsequent read data induced too many bonds per atom*
   See the extra/bond/per/atom keyword for the create_box
   or the read_data command to set this limit larger

*Subsequent read data induced too many dihedrals per atom*
   See the extra/dihedral/per/atom keyword for the create_box
   or the read_data command to set this limit larger

*Subsequent read data induced too many impropers per atom*
   See the extra/improper/per/atom keyword for the create_box
   or the read_data command to set this limit larger

*Substitution for illegal variable*
   Input script line contained a variable that could not be substituted
   for.

*Support for writing images in JPEG format not included*
   LAMMPS was not built with the -DLAMMPS_JPEG switch in the Makefile.

*Support for writing images in PNG format not included*
   LAMMPS was not built with the -DLAMMPS_PNG switch in the Makefile.

*Support for writing movies not included*
   LAMMPS was not built with the -DLAMMPS_FFMPEG switch in the Makefile

*System in data file is too big*
   See the setting for bigint in the src/lmptype.h file.

*System is not charge neutral, net charge = %g*
   The total charge on all atoms on the system is not 0.0.
   For some KSpace solvers this is an error.

*TAD nsteps must be multiple of t_event*
   Self-explanatory.

*TIP4P hydrogen has incorrect atom type*
   The TIP4P pairwise computation found an H atom whose type does not
   agree with the specified H type.

*TIP4P hydrogen is missing*
   The TIP4P pairwise computation failed to find the correct H atom
   within a water molecule.

*TMD target file did not list all group atoms*
   The target file for the fix tmd command did not list all atoms in the
   fix group.

*Tad command before simulation box is defined*
   Self-explanatory.

*Tagint setting in lmptype.h is invalid*
   Tagint must be as large or larger than smallint.

*Tagint setting in lmptype.h is not compatible*
   Format of tagint stored in restart file is not consistent with LAMMPS
   version you are running.  See the settings in src/lmptype.h

*Target pressure for fix rigid/nph cannot be < 0.0*
   Self-explanatory.

*Target pressure for fix rigid/npt/small cannot be < 0.0*
   Self-explanatory.

*Target temperature for fix nvt/npt/nph cannot be 0.0*
   Self-explanatory.

*Target temperature for fix rigid/npt cannot be 0.0*
   Self-explanatory.

*Target temperature for fix rigid/npt/small cannot be 0.0*
   Self-explanatory.

*Target temperature for fix rigid/nvt cannot be 0.0*
   Self-explanatory.

*Target temperature for fix rigid/nvt/small cannot be 0.0*
   Self-explanatory.

*Temper command before simulation box is defined*
   The temper command cannot be used before a read_data, read_restart, or
   create_box command.

*Temperature ID for fix bond/swap does not exist*
   Self-explanatory.

*Temperature ID for fix box/relax does not exist*
   Self-explanatory.

*Temperature ID for fix nvt/npt does not exist*
   Self-explanatory.

*Temperature ID for fix press/berendsen does not exist*
   Self-explanatory.

*Temperature ID for fix rigid nvt/npt/nph does not exist*
   Self-explanatory.

*Temperature ID for fix temp/berendsen does not exist*
   Self-explanatory.

*Temperature ID for fix temp/csld does not exist*
   Self-explanatory.

*Temperature ID for fix temp/csvr does not exist*
   Self-explanatory.

*Temperature ID for fix temp/rescale does not exist*
   Self-explanatory.

*Temperature compute degrees of freedom < 0*
   This should not happen if you are calculating the temperature
   on a valid set of atoms.

*Temperature control can not be used with fix nph*
   Self-explanatory.

*Temperature control can not be used with fix nph/asphere*
   Self-explanatory.

*Temperature control can not be used with fix nph/body*
   Self-explanatory.

*Temperature control can not be used with fix nph/sphere*
   Self-explanatory.

*Temperature control must be used with fix nphug*
   The temp keyword must be provided.

*Temperature control must be used with fix npt*
   Self-explanatory.

*Temperature control must be used with fix npt/asphere*
   Self-explanatory.

*Temperature control must be used with fix npt/body*
   Self-explanatory.

*Temperature control must be used with fix npt/sphere*
   Self-explanatory.

*Temperature control must be used with fix nvt*
   Self-explanatory.

*Temperature control must be used with fix nvt/asphere*
   Self-explanatory.

*Temperature control must be used with fix nvt/body*
   Self-explanatory.

*Temperature control must be used with fix nvt/sllod*
   Self-explanatory.

*Temperature control must be used with fix nvt/sphere*
   Self-explanatory.

*Temperature control must not be used with fix nph/small*
   Self-explanatory.

*Temperature for fix nvt/sllod does not have a bias*
   The specified compute must compute temperature with a bias.

*Tempering could not find thermo_pe compute*
   This compute is created by the thermo command.  It must have been
   explicitly deleted by a uncompute command.

*Tempering fix ID is not defined*
   The fix ID specified by the temper command does not exist.

*Tempering temperature fix is not valid*
   The fix specified by the temper command is not one that controls
   temperature (nvt or langevin).

*Test_descriptor_string already allocated*
   This is an internal error.  Contact the developers.

*The package gpu command is required for gpu styles*
   Self-explanatory.

*Thermo and fix not computed at compatible times*
   Fixes generate values on specific timesteps.  The thermo output
   does not match these timesteps.

*Thermo compute array is accessed out-of-range*
   Self-explanatory.

*Thermo compute does not compute array*
   Self-explanatory.

*Thermo compute does not compute scalar*
   Self-explanatory.

*Thermo compute does not compute vector*
   Self-explanatory.

*Thermo compute vector is accessed out-of-range*
   Self-explanatory.

*Thermo custom variable cannot be indexed*
   Self-explanatory.

*Thermo custom variable is not equal-style variable*
   Only equal-style variables can be output with thermodynamics, not
   atom-style variables.

*Thermo every variable returned a bad timestep*
   The variable must return a timestep greater than the current timestep.

*Thermo fix array is accessed out-of-range*
   Self-explanatory.

*Thermo fix does not compute array*
   Self-explanatory.

*Thermo fix does not compute scalar*
   Self-explanatory.

*Thermo fix does not compute vector*
   Self-explanatory.

*Thermo fix vector is accessed out-of-range*
   Self-explanatory.

*Thermo keyword in variable requires thermo to use/init pe*
   You are using a thermo keyword in a variable that requires
   potential energy to be calculated, but your thermo output
   does not use it.  Add it to your thermo output.

*Thermo keyword in variable requires thermo to use/init press*
   You are using a thermo keyword in a variable that requires pressure to
   be calculated, but your thermo output does not use it.  Add it to your
   thermo output.

*Thermo keyword in variable requires thermo to use/init temp*
   You are using a thermo keyword in a variable that requires temperature
   to be calculated, but your thermo output does not use it.  Add it to
   your thermo output.

*Thermo style does not use press*
   Cannot use thermo_modify to set this parameter since the thermo_style
   is not computing this quantity.

*Thermo style does not use temp*
   Cannot use thermo_modify to set this parameter since the thermo_style
   is not computing this quantity.

*Thermo_modify every variable returned a bad timestep*
   The returned timestep is less than or equal to the current timestep.

*Thermo_modify int format does not contain d character*
   Self-explanatory.

*Thermo_modify pressure ID does not compute pressure*
   The specified compute ID does not compute pressure.

*Thermo_modify temperature ID does not compute temperature*
   The specified compute ID does not compute temperature.

*Thermo_style command before simulation box is defined*
   The thermo_style command cannot be used before a read_data,
   read_restart, or create_box command.

*This variable thermo keyword cannot be used between runs*
   Keywords that refer to time (such as cpu, elapsed) do not
   make sense in between runs.

*Threshold for an atom property that isn't allocated*
   A dump threshold has been requested on a quantity that is
   not defined by the atom style used in this simulation.

*Timestep must be >= 0*
   Specified timestep is invalid.

*Too big a problem to use velocity create loop all*
   The system size must fit in a 32-bit integer to use this option.

*Too big a timestep for dump dcd*
   The timestep must fit in a 32-bit integer to use this dump style.

*Too big a timestep for dump xtc*
   The timestep must fit in a 32-bit integer to use this dump style.

*Too few bits for lookup table*
   Table size specified via pair_modify command does not work with your
   machine's floating point representation.

*Too few lines in %s section of data file*
   Self-explanatory.

*Too few values in body lines in data file*
   Self-explanatory.

*Too few values in body section of molecule file*
   Self-explanatory.

*Too many -pk arguments in command line*
   The string formed by concatenating the arguments is too long.  Use a
   package command in the input script instead.

*Too many MSM grid levels*
   The max number of MSM grid levels is hardwired to 10.

*Too many args in variable function*
   More args are used than any variable function allows.

*Too many atom pairs for pair bop*
   The number of atomic pairs exceeds the expected number.  Check your
   atomic structure to ensure that it is realistic.

*Too many atom sorting bins*
   This is likely due to an immense simulation box that has blown up
   to a large size.

*Too many atom triplets for pair bop*
   The number of three atom groups for angle determinations exceeds the
   expected number.  Check your atomic structure to ensure that it is
   realistic.

*Too many atoms for dump dcd*
   The system size must fit in a 32-bit integer to use this dump
   style.

*Too many atoms for dump xtc*
   The system size must fit in a 32-bit integer to use this dump
   style.

*Too many atoms to dump sort*
   Cannot sort when running with more than 2\^31 atoms.

*Too many exponent bits for lookup table*
   Table size specified via pair_modify command does not work with your
   machine's floating point representation.

*Too many groups*
   The maximum number of atom groups (including the "all" group) is
   given by MAX_GROUP in group.cpp and is 32.

*Too many iterations*
   You must use a number of iterations that fit in a 32-bit integer
   for minimization.

*Too many lines in one body in data file - boost MAXBODY*
   MAXBODY is a setting at the top of the src/read_data.cpp file.
   Set it larger and re-compile the code.

*Too many local+ghost atoms for neighbor list*
   The number of nlocal + nghost atoms on a processor
   is limited by the size of a 32-bit integer with 2 bits
   removed for masking 1-2, 1-3, 1-4 neighbors.

*Too many mantissa bits for lookup table*
   Table size specified via pair_modify command does not work with your
   machine's floating point representation.

*Too many masses for fix shake*
   The fix shake command cannot list more masses than there are atom
   types.

*Too many molecules for fix poems*
   The limit is 2\^31 = ~2 billion molecules.

*Too many molecules for fix rigid*
   The limit is 2\^31 = ~2 billion molecules.

*Too many neighbor bins*
   This is likely due to an immense simulation box that has blown up
   to a large size.

*Too many timesteps*
   The cumulative timesteps must fit in a 64-bit integer.

*Too many timesteps for NEB*
   You must use a number of timesteps that fit in a 32-bit integer
   for NEB.

*Too many total atoms*
   See the setting for bigint in the src/lmptype.h file.

*Too many total bits for bitmapped lookup table*
   Table size specified via pair_modify command is too large.  Note that
   a value of N generates a 2\^N size table.

*Too many values in body lines in data file*
   Self-explanatory.

*Too many values in body section of molecule file*
   Self-explanatory.

*Too much buffered per-proc info for dump*
   The size of the buffered string must fit in a 32-bit integer for a
   dump.

*Too much per-proc info for dump*
   Number of local atoms times number of columns must fit in a 32-bit
   integer for dump.

*Tree structure in joint connections*
   Fix poems cannot (yet) work with coupled bodies whose joints connect
   the bodies in a tree structure.

*Triclinic box skew is too large*
   The displacement in a skewed direction must be less than half the box
   length in that dimension.  E.g. the xy tilt must be between -half and
   +half of the x box length.  This constraint can be relaxed by using
   the box tilt command.

*Tried to convert a double to int, but input_double > INT_MAX*
   Self-explanatory.

*Trying to build an occasional neighbor list before initialization completed*
   This is not allowed.  Source code caller needs to be modified.

*Two fix ave commands using same compute chunk/atom command in incompatible ways*
   They are both attempting to "lock" the chunk/atom command so that the
   chunk assignments persist for some number of timesteps, but are doing
   it in different ways.

*Two groups cannot be the same in fix spring couple*
   Self-explanatory.

*Unable to initialize accelerator for use*
   There was a problem initializing an accelerator for the gpu package

*Unbalanced quotes in input line*
   No matching end double quote was found following a leading double
   quote.

*Unexpected end of -reorder file*
   Self-explanatory.

*Unexpected empty line in AngleCoeffs section*
   Read a blank line where there should be coefficient data.

*Unexpected empty line in BondCoeffs section*
   Read a blank line where there should be coefficient data.

*Unexpected empty line in DihedralCoeffs section*
   Read a blank line where there should be coefficient data.

*Unexpected empty line in ImproperCoeffs section*
   Read a blank line where there should be coefficient data.

*Unexpected empty line in PairCoeffs section*
   Read a blank line where there should be coefficient data.

*Unexpected end of custom file*
   Self-explanatory.

*Unexpected end of data file*
   LAMMPS hit the end of the data file while attempting to read a
   section.  Something is wrong with the format of the data file.

*Unexpected end of dump file*
   A read operation from the file failed.

*Unexpected end of fix rigid file*
   A read operation from the file failed.

*Unexpected end of fix rigid/small file*
   A read operation from the file failed.

*Unexpected end of molecule file*
   Self-explanatory.

*Unexpected end of neb file*
   A read operation from the file failed.

*Units command after simulation box is defined*
   The units command cannot be used after a read_data, read_restart, or
   create_box command.

*Universe/uloop variable count < # of partitions*
   A universe or uloop style variable must specify a number of values >= to the
   number of processor partitions.

*Unrecognized angle style*
   The choice of angle style is unknown.

*Unrecognized atom style*
   The choice of atom style is unknown.

*Unrecognized body style*
   The choice of body style is unknown.

*Unrecognized bond style*
   The choice of bond style is unknown.

*Unknown category for info is_active()*
   Self-explanatory.

*Unknown category for info is_available()*
   Self-explanatory.

*Unknown category for info is_defined()*
   Self-explanatory.

*Unrecognized command: %s*
   The command is not known to LAMMPS.  Check the input script.

*Unrecognized compute style*
   The choice of compute style is unknown.

*Unrecognized dihedral style*
   The choice of dihedral style is unknown.

*Unrecognized dump reader style*
   The choice of dump reader style via the format keyword is unknown.

*Unrecognized dump style*
   The choice of dump style is unknown.

*Unknown error in GPU library*
   Self-explanatory.

*Unrecognized fix style*
   The choice of fix style is unknown.

*Unknown identifier in data file: %s*
   A section of the data file cannot be read by LAMMPS.

*Unrecognized improper style*
   The choice of improper style is unknown.

*Unknown keyword in thermo_style custom command*
   One or more specified keywords are not recognized.

*Unrecognized kspace style*
   The choice of kspace style is unknown.

*Unknown name for info newton category*
   Self-explanatory.

*Unknown name for info package category*
   Self-explanatory.

*Unknown name for info pair category*
   Self-explanatory.

*Unrecognized pair style*
   The choice of pair style is unknown.

*Unknown pair_modify hybrid sub-style*
   The choice of sub-style is unknown.

*Unrecognized region style*
   The choice of region style is unknown.

*Unknown section in molecule file*
   Self-explanatory.

*Unknown table style in angle style table*
   Self-explanatory.

*Unknown table style in bond style table*
   Self-explanatory.

*Unknown table style in pair_style command*
   Style of table is invalid for use with pair_style table command.

*Unknown unit_style*
   Self-explanatory. Check the input script or data file.

*Unrecognized lattice type in MEAM file 1*
   The lattice type in an entry of the MEAM library file is not
   valid.

*Unrecognized lattice type in MEAM file 2*
   The lattice type in an entry of the MEAM parameter file is not
   valid.

*Unrecognized pair style in compute pair command*
   Self-explanatory.

*Unsupported mixing rule in kspace_style ewald/disp*
   Only geometric mixing is supported.

*Unsupported order in kspace_style ewald/disp*
   Only 1/r\^6 dispersion or dipole terms are supported.

*Unsupported order in kspace_style pppm/disp, pair_style %s*
   Only pair styles with 1/r and 1/r\^6 dependence are currently supported.

*Use cutoff keyword to set cutoff in single mode*
   Mode is single so cutoff/multi keyword cannot be used.

*Use cutoff/multi keyword to set cutoff in multi mode*
   Mode is multi so cutoff keyword cannot be used.

*Using fix nvt/sllod with inconsistent fix deform remap option*
   Fix nvt/sllod requires that deforming atoms have a velocity profile
   provided by "remap v" as a fix deform option.

*Using fix nvt/sllod with no fix deform defined*
   Self-explanatory.

*Using fix srd with inconsistent fix deform remap option*
   When shearing the box in an SRD simulation, the remap v option for fix
   deform needs to be used.

*Using pair lubricate with inconsistent fix deform remap option*
   Must use remap v option with fix deform with this pair style.

*Using pair lubricate/poly with inconsistent fix deform remap option*
   If fix deform is used, the remap v option is required.

*Using suffix gpu without GPU package installed*
   Self-explanatory.

*Using suffix intel without USER-INTEL package installed*
   Self-explanatory.

*Using suffix kk without KOKKOS package enabled*
   Self-explanatory.

*Using suffix omp without USER-OMP package installed*
   Self-explanatory.

*Using update dipole flag requires atom attribute mu*
   Self-explanatory.

*Using update dipole flag requires atom style sphere*
   Self-explanatory.

*Variable ID in variable formula does not exist*
   Self-explanatory.

*Variable atom ID is too large*
   Specified ID is larger than the maximum allowed atom ID.

*Variable evaluation before simulation box is defined*
   Cannot evaluate a compute or fix or atom-based value in a variable
   before the simulation has been setup.

*Variable evaluation in fix wall gave bad value*
   The returned value for epsilon or sigma < 0.0.

*Variable evaluation in region gave bad value*
   Variable returned a radius < 0.0.

*Variable for compute ti is invalid style*
   Self-explanatory.

*Variable for create_atoms is invalid style*
   The variables must be equal-style variables.

*Variable for displace_atoms is invalid style*
   It must be an equal-style or atom-style variable.

*Variable for dump every is invalid style*
   Only equal-style variables can be used.

*Variable for dump image center is invalid style*
   Must be an equal-style variable.

*Variable for dump image persp is invalid style*
   Must be an equal-style variable.

*Variable for dump image phi is invalid style*
   Must be an equal-style variable.

*Variable for dump image theta is invalid style*
   Must be an equal-style variable.

*Variable for dump image zoom is invalid style*
   Must be an equal-style variable.

*Variable for fix adapt is invalid style*
   Only equal-style variables can be used.

*Variable for fix addforce is invalid style*
   Self-explanatory.

*Variable for fix aveforce is invalid style*
   Only equal-style variables can be used.

*Variable for fix deform is invalid style*
   The variable must be an equal-style variable.

*Variable for fix efield is invalid style*
   The variable must be an equal- or atom-style variable.

*Variable for fix gravity is invalid style*
   Only equal-style variables can be used.

*Variable for fix heat is invalid style*
   Only equal-style or atom-style variables can be used.

*Variable for fix indent is invalid style*
   Only equal-style variables can be used.

*Variable for fix indent is not equal style*
   Only equal-style variables can be used.

*Variable for fix langevin is invalid style*
   It must be an equal-style variable.

*Variable for fix move is invalid style*
   Only equal-style variables can be used.

*Variable for fix setforce is invalid style*
   Only equal-style variables can be used.

*Variable for fix temp/berendsen is invalid style*
   Only equal-style variables can be used.

*Variable for fix temp/csld is invalid style*
   Only equal-style variables can be used.

*Variable for fix temp/csvr is invalid style*
   Only equal-style variables can be used.

*Variable for fix temp/rescale is invalid style*
   Only equal-style variables can be used.

*Variable for fix wall is invalid style*
   Only equal-style variables can be used.

*Variable for fix wall/reflect is invalid style*
   Only equal-style variables can be used.

*Variable for fix wall/srd is invalid style*
   Only equal-style variables can be used.

*Variable for group dynamic is invalid style*
   The variable must be an atom-style variable.

*Variable for group is invalid style*
   Only atom-style variables can be used.

*Variable for region cylinder is invalid style*
   Only equal-style variables are allowed.

*Variable for region is invalid style*
   Only equal-style variables can be used.

*Variable for region is not equal style*
   Self-explanatory.

*Variable for region sphere is invalid style*
   Only equal-style variables are allowed.

*Variable for restart is invalid style*
   Only equal-style variables can be used.

*Variable for set command is invalid style*
   Only atom-style variables can be used.

*Variable for thermo every is invalid style*
   Only equal-style variables can be used.

*Variable for velocity set is invalid style*
   Only atom-style variables can be used.

*Variable for voronoi radius is not atom style*
   Self-explanatory.

*Variable formula compute array is accessed out-of-range*
   Self-explanatory.

*Variable formula compute vector is accessed out-of-range*
   Self-explanatory.

*Variable formula fix array is accessed out-of-range*
   Self-explanatory.

*Variable formula fix vector is accessed out-of-range*
   Self-explanatory.

*Variable has circular dependency*
   A circular dependency is when variable "a" in used by variable "b" and
   variable "b" is also used by variable "a".  Circular dependencies with
   longer chains of dependence are also not allowed.

*Variable name between brackets must be alphanumeric or underscore characters*
   Self-explanatory.

*Variable name for compute chunk/atom does not exist*
   Self-explanatory.

*Variable name for compute reduce does not exist*
   Self-explanatory.

*Variable name for compute ti does not exist*
   Self-explanatory.

*Variable name for create_atoms does not exist*
   Self-explanatory.

*Variable name for displace_atoms does not exist*
   Self-explanatory.

*Variable name for dump every does not exist*
   Self-explanatory.

*Variable name for dump image center does not exist*
   Self-explanatory.

*Variable name for dump image persp does not exist*
   Self-explanatory.

*Variable name for dump image phi does not exist*
   Self-explanatory.

*Variable name for dump image theta does not exist*
   Self-explanatory.

*Variable name for dump image zoom does not exist*
   Self-explanatory.

*Variable name for fix adapt does not exist*
   Self-explanatory.

*Variable name for fix addforce does not exist*
   Self-explanatory.

*Variable name for fix ave/atom does not exist*
   Self-explanatory.

*Variable name for fix ave/chunk does not exist*
   Self-explanatory.

*Variable name for fix ave/correlate does not exist*
   Self-explanatory.

*Variable name for fix ave/histo does not exist*
   Self-explanatory.

*Variable name for fix ave/spatial does not exist*
   Self-explanatory.

*Variable name for fix ave/time does not exist*
   Self-explanatory.

*Variable name for fix aveforce does not exist*
   Self-explanatory.

*Variable name for fix deform does not exist*
   Self-explanatory.

*Variable name for fix efield does not exist*
   Self-explanatory.

*Variable name for fix gravity does not exist*
   Self-explanatory.

*Variable name for fix heat does not exist*
   Self-explanatory.

*Variable name for fix indent does not exist*
   Self-explanatory.

*Variable name for fix langevin does not exist*
   Self-explanatory.

*Variable name for fix move does not exist*
   Self-explanatory.

*Variable name for fix setforce does not exist*
   Self-explanatory.

*Variable name for fix store/state does not exist*
   Self-explanatory.

*Variable name for fix temp/berendsen does not exist*
   Self-explanatory.

*Variable name for fix temp/csld does not exist*
   Self-explanatory.

*Variable name for fix temp/csvr does not exist*
   Self-explanatory.

*Variable name for fix temp/rescale does not exist*
   Self-explanatory.

*Variable name for fix vector does not exist*
   Self-explanatory.

*Variable name for fix wall does not exist*
   Self-explanatory.

*Variable name for fix wall/reflect does not exist*
   Self-explanatory.

*Variable name for fix wall/srd does not exist*
   Self-explanatory.

*Variable name for group does not exist*
   Self-explanatory.

*Variable name for group dynamic does not exist*
   Self-explanatory.

*Variable name for region cylinder does not exist*
   Self-explanatory.

*Variable name for region does not exist*
   Self-explanatory.

*Variable name for region sphere does not exist*
   Self-explanatory.

*Variable name for restart does not exist*
   Self-explanatory.

*Variable name for set command does not exist*
   Self-explanatory.

*Variable name for thermo every does not exist*
   Self-explanatory.

*Variable name for velocity set does not exist*
   Self-explanatory.

*Variable name for voronoi radius does not exist*
   Self-explanatory.

*Variable name must be alphanumeric or underscore characters*
   Self-explanatory.

*Variable uses atom property that isn't allocated*
   Self-explanatory.

*Velocity command before simulation box is defined*
   The velocity command cannot be used before a read_data, read_restart,
   or create_box command.

*Velocity command with no atoms existing*
   A velocity command has been used, but no atoms yet exist.

*Velocity ramp in z for a 2d problem*
   Self-explanatory.

*Velocity rigid used with non-rigid fix-ID*
   Self-explanatory.

*Velocity temperature ID does calculate a velocity bias*
   The specified compute must compute a bias for temperature.

*Velocity temperature ID does not compute temperature*
   The compute ID given to the velocity command must compute
   temperature.

*Verlet/split can only currently be used with comm_style brick*
   This is a current restriction in LAMMPS.

*Verlet/split does not yet support TIP4P*
   This is a current limitation.

*Verlet/split requires 2 partitions*
   See the -partition command-line switch.

*Verlet/split requires Rspace partition layout be multiple of Kspace partition layout in each dim*
   This is controlled by the processors command.

*Verlet/split requires Rspace partition size be multiple of Kspace partition size*
   This is so there is an equal number of Rspace processors for every
   Kspace processor.

*Virial was not tallied on needed timestep*
   You are using a thermo keyword that requires potentials to
   have tallied the virial, but they didn't on this timestep.  See the
   variable doc page for ideas on how to make this work.

*Voro++ error: narea and neigh have a different size*
   This error is returned by the Voro++ library.

*Wall defined twice in fix wall command*
   Self-explanatory.

*Wall defined twice in fix wall/reflect command*
   Self-explanatory.

*Wall defined twice in fix wall/srd command*
   Self-explanatory.

*Water H epsilon must be 0.0 for pair style lj/cut/tip4p/cut*
   This is because LAMMPS does not compute the Lennard-Jones interactions
   with these particles for efficiency reasons.

*Water H epsilon must be 0.0 for pair style lj/cut/tip4p/long*
   This is because LAMMPS does not compute the Lennard-Jones interactions
   with these particles for efficiency reasons.

*Water H epsilon must be 0.0 for pair style lj/long/tip4p/long*
   This is because LAMMPS does not compute the Lennard-Jones interactions
   with these particles for efficiency reasons.

*World variable count doesn't match # of partitions*
   A world-style variable must specify a number of values equal to the
   number of processor partitions.

*Write_data command before simulation box is defined*
   Self-explanatory.

*Write_restart command before simulation box is defined*
   The write_restart command cannot be used before a read_data,
   read_restart, or create_box command.

*Writing to MPI-IO filename when MPIIO package is not installed*
   Self-explanatory.

*Zero length rotation vector with displace_atoms*
   Self-explanatory.

*Zero length rotation vector with fix move*
   Self-explanatory.

*Zero-length lattice orient vector*
   Self-explanatory.
