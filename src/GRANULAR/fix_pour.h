/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(pour,FixPour)

#else

#ifndef LMP_FIX_POUR_H
#define LMP_FIX_POUR_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPour : public Fix {
 public:
  FixPour(class LAMMPS *, int, char **);
  ~FixPour();
  int setmask();
  void init();
  void pre_exchange();
  void reset_dt();
  void *extract(const char *, int &);

 private:
  int ninsert,ntype,seed;
  int iregion,mode,idnext,dstyle,npoly,rigidflag,shakeflag;
  double radius_one,radius_max;
  double radius_lo,radius_hi;
  double *radius_poly,*frac_poly;
  double density_lo,density_hi;
  double volfrac;
  int maxattempt;
  int region_style;
  double rate;
  double vxlo,vxhi,vylo,vyhi,vy,vz;
  double xlo,xhi,ylo,yhi,zlo,zhi;
  double xc,yc,rc;
  double grav;
  char *idrigid,*idshake;

  class Molecule *onemol;
  int natom;
  double **coords;
  imageint *imageflags;
  class Fix *fixrigid,*fixshake;
  double oneradius;

  int me,nprocs;
  int *recvcounts,*displs;
  int nfreq,nfirst,ninserted,nper;
  double lo_current,hi_current;
  tagint maxtag_all,maxmol_all;
  class RanPark *random,*random2;

  void find_maxid();
  int overlap(int);
  int outside(int, double, double, double);
  void xyz_random(double, double *);
  double radius_sample();
  void options(int, char **);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix pour requires atom attributes radius, rmass

UNDOCUMENTED

E: Invalid atom type in fix pour command

UNDOCUMENTED

E: Must specify a region in fix pour

UNDOCUMENTED

E: Fix pour region does not support a bounding box

UNDOCUMENTED

E: Fix pour region cannot be dynamic

UNDOCUMENTED

E: Insertion region extends outside simulation box

UNDOCUMENTED

E: Must use a z-axis cylinder with fix pour

UNDOCUMENTED

E: Must use a block or cylinder region with fix pour

UNDOCUMENTED

E: Must use a block region with fix pour for 2d simulations

UNDOCUMENTED

E: Cannot use fix_pour unless atoms have IDs

UNDOCUMENTED

E: Fix pour molecule must have coordinates

UNDOCUMENTED

E: Fix pour molecule must have atom types

UNDOCUMENTED

E: Invalid atom type in fix pour mol command

UNDOCUMENTED

E: Fix pour molecule template ID must be same as atom style template ID

UNDOCUMENTED

E: Cannot use fix pour rigid and not molecule

UNDOCUMENTED

E: Cannot use fix pour shake and not molecule

UNDOCUMENTED

E: Cannot use fix pour rigid and shake

UNDOCUMENTED

E: No fix gravity defined for fix pour

UNDOCUMENTED

E: Cannot use fix pour with triclinic box

UNDOCUMENTED

E: Gravity must point in -z to use with fix pour in 3d

UNDOCUMENTED

E: Gravity must point in -y to use with fix pour in 2d

UNDOCUMENTED

E: Gravity changed since fix pour was created

UNDOCUMENTED

E: Fix pour rigid fix does not exist

UNDOCUMENTED

E: Fix pour and fix rigid/small not using same molecule template ID

UNDOCUMENTED

E: Fix pour shake fix does not exist

UNDOCUMENTED

E: Fix pour and fix shake not using same molecule template ID

UNDOCUMENTED

W: Less insertions than requested

UNDOCUMENTED

E: Too many total atoms

UNDOCUMENTED

E: New atom IDs exceed maximum allowed ID

UNDOCUMENTED

E: Fix pour region ID does not exist

UNDOCUMENTED

E: Molecule template ID for fix pour does not exist

UNDOCUMENTED

W: Molecule template for fix pour has multiple molecules

UNDOCUMENTED

E: Fix pour polydisperse fractions do not sum to 1.0

UNDOCUMENTED

E: Cannot change timestep with fix pour

UNDOCUMENTED

U: Cannot read_data after simulation box is defined

The read_data command cannot be used after a read_data,
read_restart, or create_box command.

U: Cannot run 2d simulation with nonperiodic Z dimension

Use the boundary command to make the z dimension periodic in order to
run a 2d simulation.

U: Fix ID for read_data does not exist

Self-explanatory.

U: Must read Atoms before Velocities

The Atoms section of a data file must come before a Velocities
section.

U: Invalid data file section: Bonds

Atom style does not allow bonds.

U: Must read Atoms before Bonds

The Atoms section of a data file must come before a Bonds section.

U: Invalid data file section: Angles

Atom style does not allow angles.

U: Must read Atoms before Angles

The Atoms section of a data file must come before an Angles section.

U: Invalid data file section: Dihedrals

Atom style does not allow dihedrals.

U: Must read Atoms before Dihedrals

The Atoms section of a data file must come before a Dihedrals section.

U: Invalid data file section: Impropers

Atom style does not allow impropers.

U: Must read Atoms before Impropers

The Atoms section of a data file must come before an Impropers
section.

U: Invalid data file section: Ellipsoids

Atom style does not allow ellipsoids.

U: Must read Atoms before Ellipsoids

The Atoms section of a data file must come before a Ellipsoids
section.

U: Invalid data file section: Lines

Atom style does not allow lines.

U: Must read Atoms before Lines

The Atoms section of a data file must come before a Lines section.

U: Invalid data file section: Triangles

Atom style does not allow triangles.

U: Must read Atoms before Triangles

The Atoms section of a data file must come before a Triangles section.

U: Invalid data file section: Bodies

Atom style does not allow bodies.

U: Must read Atoms before Bodies

The Atoms section of a data file must come before a Bodies section.

U: Must define pair_style before Pair Coeffs

Must use a pair_style command before reading a data file that defines
Pair Coeffs.

U: Must define pair_style before PairIJ Coeffs

UNDOCUMENTED

U: Invalid data file section: Bond Coeffs

Atom style does not allow bonds.

U: Must define bond_style before Bond Coeffs

Must use a bond_style command before reading a data file that
defines Bond Coeffs.

U: Invalid data file section: Angle Coeffs

Atom style does not allow angles.

U: Must define angle_style before Angle Coeffs

Must use an angle_style command before reading a data file that
defines Angle Coeffs.

U: Invalid data file section: Dihedral Coeffs

Atom style does not allow dihedrals.

U: Must define dihedral_style before Dihedral Coeffs

Must use a dihedral_style command before reading a data file that
defines Dihedral Coeffs.

U: Invalid data file section: Improper Coeffs

Atom style does not allow impropers.

U: Must define improper_style before Improper Coeffs

Must use an improper_style command before reading a data file that
defines Improper Coeffs.

U: Invalid data file section: BondBond Coeffs

Atom style does not allow angles.

U: Must define angle_style before BondBond Coeffs

Must use an angle_style command before reading a data file that
defines Angle Coeffs.

U: Invalid data file section: BondAngle Coeffs

Atom style does not allow angles.

U: Must define angle_style before BondAngle Coeffs

Must use an angle_style command before reading a data file that
defines Angle Coeffs.

U: Invalid data file section: MiddleBondTorsion Coeffs

Atom style does not allow dihedrals.

U: Must define dihedral_style before MiddleBondTorsion Coeffs

Must use a dihedral_style command before reading a data file that
defines MiddleBondTorsion Coeffs.

U: Invalid data file section: EndBondTorsion Coeffs

Atom style does not allow dihedrals.

U: Must define dihedral_style before EndBondTorsion Coeffs

Must use a dihedral_style command before reading a data file that
defines EndBondTorsion Coeffs.

U: Invalid data file section: AngleTorsion Coeffs

Atom style does not allow dihedrals.

U: Must define dihedral_style before AngleTorsion Coeffs

Must use a dihedral_style command before reading a data file that
defines AngleTorsion Coeffs.

U: Invalid data file section: AngleAngleTorsion Coeffs

Atom style does not allow dihedrals.

U: Must define dihedral_style before AngleAngleTorsion Coeffs

Must use a dihedral_style command before reading a data file that
defines AngleAngleTorsion Coeffs.

U: Invalid data file section: BondBond13 Coeffs

Atom style does not allow dihedrals.

U: Must define dihedral_style before BondBond13 Coeffs

Must use a dihedral_style command before reading a data file that
defines BondBond13 Coeffs.

U: Invalid data file section: AngleAngle Coeffs

Atom style does not allow impropers.

U: Must define improper_style before AngleAngle Coeffs

Must use an improper_style command before reading a data file that
defines AngleAngle Coeffs.

U: Unknown identifier in data file: %s

A section of the data file cannot be read by LAMMPS.

U: No atoms in data file

The header of the data file indicated that atoms would be included,
but they were not present.

U: Needed molecular topology not in data file

UNDOCUMENTED

U: Needed bonus data not in data file

Some atom styles require bonus data.  See the read_data doc page for
details.

U: Unexpected end of data file

LAMMPS hit the end of the data file while attempting to read a
section.  Something is wrong with the format of the data file.

U: No ellipsoids allowed with this atom style

Self-explanatory.  Check data file.

U: No lines allowed with this atom style

Self-explanatory.  Check data file.

U: No triangles allowed with this atom style

Self-explanatory.  Check data file.

U: No bodies allowed with this atom style

Self-explanatory.  Check data file.

U: System in data file is too big

See the setting for bigint in the src/lmptype.h file.

U: No bonds allowed with this atom style

Self-explanatory.  Check data file.

U: No angles allowed with this atom style

Self-explanatory.  Check data file.

U: No dihedrals allowed with this atom style

Self-explanatory.  Check data file.

U: No impropers allowed with this atom style

Self-explanatory.  Check data file.

U: Bonds defined but no bond types

The data file header lists bonds but no bond types.

U: Angles defined but no angle types

The data file header lists angles but no angle types.

U: Dihedrals defined but no dihedral types

The data file header lists dihedrals but no dihedral types.

U: Impropers defined but no improper types

The data file header lists improper but no improper types.

U: No molecule topology allowed with atom style template

UNDOCUMENTED

U: Did not assign all atoms correctly

Atoms read in from a data file were not assigned correctly to
processors.  This is likely due to some atom coordinates being
outside a non-periodic simulation box.

U: Bonds assigned incorrectly

Bonds read in from the data file were not assigned correctly to atoms.
This means there is something invalid about the topology definitions.

U: Angles assigned incorrectly

Angles read in from the data file were not assigned correctly to
atoms.  This means there is something invalid about the topology
definitions.

U: Dihedrals assigned incorrectly

Dihedrals read in from the data file were not assigned correctly to
atoms.  This means there is something invalid about the topology
definitions.

U: Impropers assigned incorrectly

Impropers read in from the data file were not assigned correctly to
atoms.  This means there is something invalid about the topology
definitions.

U: Too many lines in one body in data file - boost MAXBODY

MAXBODY is a setting at the top of the src/read_data.cpp file.
Set it larger and re-compile the code.

U: Cannot open gzipped file

LAMMPS was compiled without support for reading and writing gzipped
files through a pipeline to the gzip program with -DLAMMPS_GZIP.

U: Cannot open file %s

The specified file cannot be opened.  Check that the path and name are
correct. If the file is a compressed file, also check that the gzip
executable can be found and run.

U: Invalid atom ID in Atoms section of data file

Atom IDs must be positive integers.

U: Molecular data file has too many atoms

These kids of data files are currently limited to a number
of atoms that fits in a 32-bit integer.

U: Needed topology not in data file

The header of the data file indicated that bonds or angles or
dihedrals or impropers would be included, but they were not present.

*/
