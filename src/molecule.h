/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_ONE_MOLECULE_H
#define LMP_ONE_MOLECULE_H

#include "pointers.h"

namespace LAMMPS_NS {

class Molecule : protected Pointers {
 public:
  char *id;    // template id of this molecule, same for all molecules in set
  int nset;    // if first in set, # of molecules in this set
               // else 0 if not first in set
  int last;    // 1 if last molecule in set, else 0

  // number of atoms,bonds,etc in molecule
  // nibody,ndbody = # of integer/double fields in body

  int natoms;
  int nbonds, nangles, ndihedrals, nimpropers;
  int ntypes, nmolecules, nfragments;
  int nbondtypes, nangletypes, ndihedraltypes, nimpropertypes;
  int nibody, ndbody;

  // max bond,angle,etc per atom

  int bond_per_atom, angle_per_atom, dihedral_per_atom, improper_per_atom;
  int maxspecial;

  // 1 if attribute defined in file, 0 if not

  int xflag, typeflag, moleculeflag, fragmentflag, qflag, radiusflag, rmassflag;
  int bondflag, angleflag, dihedralflag, improperflag;
  int nspecialflag, specialflag;
  int shakeflag, shakeflagflag, shakeatomflag, shaketypeflag;
  int bodyflag, ibodyflag, dbodyflag;

  // 1 if attribute defined or computed, 0 if not

  int centerflag, massflag, comflag, inertiaflag;

  // 1 if molecule fields require atom IDs

  int tag_require;

  // attributes

  double **x;          // displacement of each atom from origin
  int *type;           // type of each atom
  tagint *molecule;    // molecule of each atom
  double *q;           // charge on each atom
  double *radius;      // radius of each atom
  double *rmass;       // mass of each atom

  int *num_bond;    // bonds, angles, dihedrals, impropers for each atom
  int **bond_type;
  tagint **bond_atom;

  int *num_angle;
  int **angle_type;
  tagint **angle_atom1, **angle_atom2, **angle_atom3;

  int *num_dihedral;
  int **dihedral_type;
  tagint **dihedral_atom1, **dihedral_atom2, **dihedral_atom3, **dihedral_atom4;

  int *num_improper;
  int **improper_type;
  tagint **improper_atom1, **improper_atom2, **improper_atom3, **improper_atom4;

  int **nspecial;
  tagint **special;

  int *shake_flag;
  tagint **shake_atom;
  int **shake_type;

  class AtomVecBody *avec_body;
  int *ibodyparams;    // integer and double body params
  double *dbodyparams;

  // fragment info

  int **fragmentmask;    // nfragments by natoms
  std::vector<std::string> fragmentnames;

  double center[3];              // geometric center of molecule
  double masstotal;              // total mass of molecule
  double com[3];                 // center of mass of molecule
  double itensor[6];             // moments of inertia of molecule
  double inertia[3];             // principal moments of inertia of molecule
  double ex[3], ey[3], ez[3];    // principal axes of molecule in space coords
  double quat[4];                // quaternion for orientation of molecule

  double maxradius;    // max radius of any atom in molecule
  double molradius;    // radius of molecule from geometric center
                       // including finite-size particle radii
  int comatom;         // index (1-Natom) of atom closest to COM
  double maxextent;    // furthest any atom in molecule is from comatom

  double **dx;        // displacement of each atom relative to center
  double **dxcom;     // displacement of each atom relative to COM
  double **dxbody;    // displacement of each atom relative to COM
                      // in body frame (diagonalized interia tensor)

  double *quat_external;    // orientation imposed by external class
                            // e.g. FixPour or CreateAtoms

  Molecule(class LAMMPS *, int, char **, int &);
  ~Molecule();
  void compute_center();
  void compute_mass();
  void compute_com();
  void compute_inertia();
  int findfragment(const char *);
  void check_attributes(int);

 private:
  int me;
  FILE *fp;
  int *count;
  int toffset, boffset, aoffset, doffset, ioffset;
  int autospecial;
  double sizescale;

  void read(int);
  void coords(char *);
  void types(char *);
  void molecules(char *);
  void fragments(char *);
  void charges(char *);
  void diameters(char *);
  void masses(char *);
  void bonds(int, char *);
  void angles(int, char *);
  void dihedrals(int, char *);
  void impropers(int, char *);
  void nspecial_read(int, char *);
  void special_read(char *);
  void special_generate();
  void shakeflag_read(char *);
  void shakeatom_read(char *);
  void shaketype_read(char *);
  void body(int, int, char *);

  void initialize();
  void allocate();
  void deallocate();

  void readline(char *);
  std::string parse_keyword(int, char *);
  void skip_lines(int, char *, const std::string &);

  // void print();
};

}    // namespace LAMMPS_NS

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Molecule template ID must be alphanumeric or underscore characters

Self-explanatory.

E: Insufficient Jacobi rotations for rigid molecule

Eigensolve for rigid body was not sufficiently accurate.

E: Unexpected end of molecule file

Self-explanatory.

E: Molecule file z center-of-mass must be 0.0 for 2d

Self-explanatory.

E: Molecule file requires atom style body

Self-explanatory.

E: Invalid header in molecule file

UNDOCUMENTED

E: No count or invalid atom count in molecule file

The number of atoms must be specified.

E: Invalid bond count in molecule file

Self-explanatory.

E: Invalid angle count in molecule file

Self-explanatory.

E: Invalid dihedral count in molecule file

Self-explanatory.

E: Invalid improper count in molecule file

Self-explanatory.

E: Molecule file has bonds but no nbonds setting

Self-explanatory.

E: Molecule file has angles but no nangles setting

Self-explanatory.

E: Molecule file has dihedrals but no ndihedrals setting

Self-explanatory.

E: Molecule file has impropers but no nimpropers setting

Self-explanatory.

E: Molecule file has fragments but no nfragments setting

Self-explanatory.

E: Molecule file shake flags not before shake atoms

The order of the two sections is important.

E: Molecule file shake flags not before shake bonds

The order of the two sections is important.

E: Molecule file has body params but no setting for them

Self-explanatory.

E: Unknown section in molecule file

Self-explanatory.

E: Molecule file needs both Special Bond sections

Self-explanatory.

E: Molecule file has special flags but no bonds

Self-explanatory.

E: Molecule file shake info is incomplete

All 3 SHAKE sections are needed.

E: Molecule file has no Body Integers section

Self-explanatory.

E: Molecule file has no Body Doubles section

Self-explanatory.

E: Molecule file has no Fragments section

Self-explanatory.

E: Cannot auto-generate special bonds before simulation box is defined

UNDOCUMENTED

E: Molecule natoms must be 1 for body particle

Self-explanatory.

E: Molecule sizescale must be 1.0 for body particle

Self-explanatory.

E: Invalid Coords section in molecule file

Self-explanatory.

E: Molecule file z coord must be 0.0 for 2d

Self-explanatory.

E: Invalid Types section in molecule file

Self-explanatory.

E: Invalid atom type in molecule file

Atom types must range from 1 to specified # of types.

E: Invalid Charges section in molecule file

Self-explanatory.

E: Invalid Diameters section in molecule file

Self-explanatory.

E: Invalid atom diameter in molecule file

Diameters must be >= 0.0.

E: Invalid Masses section in molecule file

Self-explanatory.

E: Invalid atom mass in molecule file

Masses must be > 0.0.

E: Invalid Bonds section in molecule file

Self-explanatory.

E: Invalid atom ID in Bonds section of molecule file

Self-explanatory.

E: Invalid bond type in Bonds section of molecule file

Self-explanatory.

E: Invalid Angles section in molecule file

Self-explanatory.

E: Invalid atom ID in Angles section of molecule file

Self-explanatory.

E: Invalid angle type in Angles section of molecule file

Self-explanatory.

E: Invalid Dihedrals section in molecule file

Self-explanatory.

E: Invalid atom ID in dihedrals section of molecule file

Self-explanatory.

E: Invalid dihedral type in dihedrals section of molecule file

Self-explanatory.

E: Invalid Impropers section in molecule file

Self-explanatory.

E: Invalid atom ID in impropers section of molecule file

Self-explanatory.

E: Invalid improper type in impropers section of molecule file

Self-explanatory.

E: Invalid molecule ID in molecule file

Molecule ID must be a non-zero positive integer.

E: Invalid Molecules section in molecule file

Self-explanatory.

E: Invalid atom ID in Fragments section of molecule file

Self-explanatory.

E: Invalid Special Bond Counts section in molecule file

Self-explanatory.

E: Molecule file special list does not match special count

The number of values in an atom's special list does not match count.

E: Invalid special atom index in molecule file

Self-explanatory.

E: Molecule auto special bond generation overflow

Counts exceed maxspecial setting for other atoms in system.

E: Invalid Shake Flags section in molecule file

UNDOCUMENTED

E: Invalid shake flag in molecule file

Self-explanatory.

E: Invalid shake atom in molecule file

Self-explanatory.

E: Invalid shake type data in molecule file

UNDOCUMENTED

E: Invalid shake bond type in molecule file

Self-explanatory.

E: Invalid shake angle type in molecule file

Self-explanatory.

E: Too few values in body section of molecule file

Self-explanatory.

E: Too many values in body section of molecule file

Self-explanatory.

W: Molecule attributes do not match system attributes

An attribute is specified (e.g. diameter, charge) that is
not defined for the specified atom style.

E: Molecule topology type exceeds system topology type

The number of bond, angle, etc types in the molecule exceeds the
system setting.  See the create_box command for how to specify these
values.

E: Molecule topology/atom exceeds system topology/atom

The number of bonds, angles, etc per-atom in the molecule exceeds the
system setting.  See the create_box command for how to specify these
values.

W: Molecule has bond topology but no special bond settings

This means the bonded atoms will not be excluded in pair-wise
interactions.

E: Cannot open molecule file %s

The specified file cannot be opened.  Check that the path and name are
correct.

*/
