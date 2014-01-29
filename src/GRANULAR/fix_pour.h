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

E: Cannot read_data after simulation box is defined

The read_data command cannot be used after a read_data,
read_restart, or create_box command.

E: Cannot run 2d simulation with nonperiodic Z dimension

Use the boundary command to make the z dimension periodic in order to
run a 2d simulation.

E: Fix ID for read_data does not exist

Self-explanatory.

E: Must read Atoms before Velocities

The Atoms section of a data file must come before a Velocities
section.

E: Invalid data file section: Bonds

Atom style does not allow bonds.

E: Must read Atoms before Bonds

The Atoms section of a data file must come before a Bonds section.

E: Invalid data file section: Angles

Atom style does not allow angles.

E: Must read Atoms before Angles

The Atoms section of a data file must come before an Angles section.

E: Invalid data file section: Dihedrals

Atom style does not allow dihedrals.

E: Must read Atoms before Dihedrals

The Atoms section of a data file must come before a Dihedrals section.

E: Invalid data file section: Impropers

Atom style does not allow impropers.

E: Must read Atoms before Impropers

The Atoms section of a data file must come before an Impropers
section.

E: Invalid data file section: Ellipsoids

Atom style does not allow ellipsoids.

E: Must read Atoms before Ellipsoids

The Atoms section of a data file must come before a Ellipsoids
section.

E: Invalid data file section: Lines

Atom style does not allow lines.

E: Must read Atoms before Lines

The Atoms section of a data file must come before a Lines section.

E: Invalid data file section: Triangles

Atom style does not allow triangles.

E: Must read Atoms before Triangles

The Atoms section of a data file must come before a Triangles section.

E: Invalid data file section: Bodies

Atom style does not allow bodies.

E: Must read Atoms before Bodies

The Atoms section of a data file must come before a Bodies section.

E: Must define pair_style before Pair Coeffs

Must use a pair_style command before reading a data file that defines
Pair Coeffs.

E: Must define pair_style before PairIJ Coeffs

UNDOCUMENTED

E: Invalid data file section: Bond Coeffs

Atom style does not allow bonds.

E: Must define bond_style before Bond Coeffs

Must use a bond_style command before reading a data file that
defines Bond Coeffs.

E: Invalid data file section: Angle Coeffs

Atom style does not allow angles.

E: Must define angle_style before Angle Coeffs

Must use an angle_style command before reading a data file that
defines Angle Coeffs.

E: Invalid data file section: Dihedral Coeffs

Atom style does not allow dihedrals.

E: Must define dihedral_style before Dihedral Coeffs

Must use a dihedral_style command before reading a data file that
defines Dihedral Coeffs.

E: Invalid data file section: Improper Coeffs

Atom style does not allow impropers.

E: Must define improper_style before Improper Coeffs

Must use an improper_style command before reading a data file that
defines Improper Coeffs.

E: Invalid data file section: BondBond Coeffs

Atom style does not allow angles.

E: Must define angle_style before BondBond Coeffs

Must use an angle_style command before reading a data file that
defines Angle Coeffs.

E: Invalid data file section: BondAngle Coeffs

Atom style does not allow angles.

E: Must define angle_style before BondAngle Coeffs

Must use an angle_style command before reading a data file that
defines Angle Coeffs.

E: Invalid data file section: MiddleBondTorsion Coeffs

Atom style does not allow dihedrals.

E: Must define dihedral_style before MiddleBondTorsion Coeffs

Must use a dihedral_style command before reading a data file that
defines MiddleBondTorsion Coeffs.

E: Invalid data file section: EndBondTorsion Coeffs

Atom style does not allow dihedrals.

E: Must define dihedral_style before EndBondTorsion Coeffs

Must use a dihedral_style command before reading a data file that
defines EndBondTorsion Coeffs.

E: Invalid data file section: AngleTorsion Coeffs

Atom style does not allow dihedrals.

E: Must define dihedral_style before AngleTorsion Coeffs

Must use a dihedral_style command before reading a data file that
defines AngleTorsion Coeffs.

E: Invalid data file section: AngleAngleTorsion Coeffs

Atom style does not allow dihedrals.

E: Must define dihedral_style before AngleAngleTorsion Coeffs

Must use a dihedral_style command before reading a data file that
defines AngleAngleTorsion Coeffs.

E: Invalid data file section: BondBond13 Coeffs

Atom style does not allow dihedrals.

E: Must define dihedral_style before BondBond13 Coeffs

Must use a dihedral_style command before reading a data file that
defines BondBond13 Coeffs.

E: Invalid data file section: AngleAngle Coeffs

Atom style does not allow impropers.

E: Must define improper_style before AngleAngle Coeffs

Must use an improper_style command before reading a data file that
defines AngleAngle Coeffs.

E: Unknown identifier in data file: %s

A section of the data file cannot be read by LAMMPS.

E: No atoms in data file

The header of the data file indicated that atoms would be included,
but they were not present.

E: Needed molecular topology not in data file

UNDOCUMENTED

E: Needed bonus data not in data file

Some atom styles require bonus data.  See the read_data doc page for
details.

E: Unexpected end of data file

LAMMPS hit the end of the data file while attempting to read a
section.  Something is wrong with the format of the data file.

E: No ellipsoids allowed with this atom style

Self-explanatory.  Check data file.

E: No lines allowed with this atom style

Self-explanatory.  Check data file.

E: No triangles allowed with this atom style

Self-explanatory.  Check data file.

E: No bodies allowed with this atom style

Self-explanatory.  Check data file.

E: System in data file is too big

See the setting for bigint in the src/lmptype.h file.

E: No bonds allowed with this atom style

Self-explanatory.  Check data file.

E: No angles allowed with this atom style

Self-explanatory.  Check data file.

E: No dihedrals allowed with this atom style

Self-explanatory.  Check data file.

E: No impropers allowed with this atom style

Self-explanatory.  Check data file.

E: Bonds defined but no bond types

The data file header lists bonds but no bond types.

E: Angles defined but no angle types

The data file header lists angles but no angle types.

E: Dihedrals defined but no dihedral types

The data file header lists dihedrals but no dihedral types.

E: Impropers defined but no improper types

The data file header lists improper but no improper types.

E: No molecule topology allowed with atom style template

UNDOCUMENTED

E: Did not assign all atoms correctly

Atoms read in from a data file were not assigned correctly to
processors.  This is likely due to some atom coordinates being
outside a non-periodic simulation box.

E: Bonds assigned incorrectly

Bonds read in from the data file were not assigned correctly to atoms.
This means there is something invalid about the topology definitions.

E: Angles assigned incorrectly

Angles read in from the data file were not assigned correctly to
atoms.  This means there is something invalid about the topology
definitions.

E: Dihedrals assigned incorrectly

Dihedrals read in from the data file were not assigned correctly to
atoms.  This means there is something invalid about the topology
definitions.

E: Impropers assigned incorrectly

Impropers read in from the data file were not assigned correctly to
atoms.  This means there is something invalid about the topology
definitions.

E: Too many lines in one body in data file - boost MAXBODY

MAXBODY is a setting at the top of the src/read_data.cpp file.
Set it larger and re-compile the code.

E: Cannot open gzipped file

LAMMPS was compiled without support for reading and writing gzipped
files through a pipeline to the gzip program with -DLAMMPS_GZIP.

E: Cannot open file %s

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
