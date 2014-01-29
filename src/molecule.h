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

#ifndef LMP_ONE_MOLECULE_H
#define LMP_ONE_MOLEUCULE_H

#include "pointers.h"

namespace LAMMPS_NS {

class Molecule : protected Pointers {
 public:
  char *id;   // template id of this molecule, same for all molecules in set
  int nset;   // if first in set, # of molecules in this set
              // else 0 if not first in set

  // number of atoms,bonds,etc in molecule

  int natoms;
  int nbonds,nangles,ndihedrals,nimpropers;
  int ntypes;
  int nbondtypes,nangletypes,ndihedraltypes,nimpropertypes;

  // max bond,angle,etc per atom

  int bond_per_atom,angle_per_atom,dihedral_per_atom,improper_per_atom;
  int maxspecial;

  // 1 if attribute defined in file, 0 if not

  int xflag,typeflag,qflag,radiusflag,rmassflag;
  int bondflag,angleflag,dihedralflag,improperflag;
  int nspecialflag,specialflag;
  int shakeflag,shakeflagflag,shakeatomflag,shaketypeflag;

  // 1 if attribute defined or computed, 0 if not

  int centerflag,massflag,comflag,inertiaflag;

  // 1 if molecule fields require atom IDs

  int tag_require;

  // attributes

  double **x;          // displacement of each atom from origin
  int *type;           // type of each atom
  double *q;           // charge on each atom
  double *radius;      // radius of each atom
  double *rmass;       // mass of each atom

  int *num_bond;       // bonds, angles, dihedrals, impropers for each atom
  int **bond_type;
  tagint **bond_atom;
  
  int *num_angle;
  int **angle_type;
  tagint **angle_atom1,**angle_atom2,**angle_atom3;
  
  int *num_dihedral;
  int **dihedral_type;
  tagint **dihedral_atom1,**dihedral_atom2,**dihedral_atom3,**dihedral_atom4;
  
  int *num_improper;
  int **improper_type;
  tagint **improper_atom1,**improper_atom2,**improper_atom3,**improper_atom4;

  int **nspecial;
  tagint **special;

  int *shake_flag;
  tagint **shake_atom;
  int **shake_type;

  double center[3];         // geometric center of molecule
  double masstotal;         // total mass of molecule
  double com[3];            // center of mass of molecule
  double itensor[6];        // moments of inertia of molecule
  double inertia[3];        // principal moments of inertia of molecule
  double ex[3],ey[3],ez[3]; // principal axes of molecule in space coords
  double quat[4];           // quaternion for orientation of molecule

  double molradius;    // radius of molecule from COM,
                       // including finite-size particle radii 
  int comatom;         // index (1-Natom) of atom closest to COM
  double maxextent;    // furthest any atom in molecule is from comatom

  double **dx;         // displacement of each atom relative to center
  double **dxcom;      // displacement of each atom relative to COM
  double **dxbody;     // displacement of each atom relative to COM
                       // in body frame (diagonalized interia tensor)

  Molecule(class LAMMPS *, char *, char *);
  ~Molecule();
  void compute_center();
  void compute_mass();
  void compute_com();
  void compute_inertia();
  void check_attributes(int);

 private:
  int me;
  FILE *fp;
  int *count;

  void read(int);
  void coords(char *);
  void types(char *);
  void charges(char *);
  void diameters(char *);
  void masses(char *);
  void bonds(int, char *);
  void angles(int, char *);
  void dihedrals(int, char *);
  void impropers(int, char *);
  void nspecial_read(int, char *);
  void special_read(char *);
  void shakeflag_read(char *);
  void shakeatom_read(char *);
  void shaketype_read(char *);

  void initialize();
  void allocate();
  void deallocate();

  void open(char *);
  void readline(char *);
  void parse_keyword(int, char *, char *);
  void skip_lines(int, char *);
  int parse(char *, char **, int);

  // void print();
};

}

#endif

/* ERROR/WARNING messages:

E: Molecule template ID must be alphanumeric or underscore characters

UNDOCUMENTED

E: Insufficient Jacobi rotations for rigid molecule

UNDOCUMENTED

E: Unexpected end of molecule file

UNDOCUMENTED

E: Molecule file z center-of-mass must be 0.0 for 2d

UNDOCUMENTED

E: No atom count in molecule file

UNDOCUMENTED

E: Molecule file has bonds but no nbonds setting

UNDOCUMENTED

E: Molecule file has angles but no nangles setting

UNDOCUMENTED

E: Molecule file has dihedrals but no ndihedrals setting

UNDOCUMENTED

E: Molecule file has impropers but no nimpropers setting

UNDOCUMENTED

E: Molecule file shake flags not before shake atoms

UNDOCUMENTED

E: Molecule file shake flags not before shake bonds

UNDOCUMENTED

E: Unknown section in molecule file

UNDOCUMENTED

E: Molecule file needs both Special Bond sections

UNDOCUMENTED

E: Molecule file has special flags but no bonds

UNDOCUMENTED

E: Molecule file shake info is incomplete

UNDOCUMENTED

E: Molecule file z coord must be 0.0 for 2d

UNDOCUMENTED

E: Invalid atom type in molecule file

UNDOCUMENTED

E: Invalid atom diameter in molecule file

UNDOCUMENTED

E: Invalid atom mass in molecule file

UNDOCUMENTED

E: Invalid atom ID in Bonds section of molecule file

UNDOCUMENTED

E: Invalid bond type in Bonds section of molecule file

UNDOCUMENTED

E: Invalid atom ID in Angles section of molecule file

UNDOCUMENTED

E: Invalid angle type in Angles section of molecule file

UNDOCUMENTED

E: Invalid atom ID in dihedrals section of molecule file

UNDOCUMENTED

E: Invalid dihedral type in dihedrals section of molecule file

UNDOCUMENTED

E: Invalid atom ID in impropers section of molecule file

UNDOCUMENTED

E: Invalid improper type in impropers section of molecule file

UNDOCUMENTED

E: Molecule file special list does not match special count

UNDOCUMENTED

E: Invalid special atom index in molecule file

UNDOCUMENTED

E: Invalid shake flag in molecule file

UNDOCUMENTED

E: Invalid shake atom in molecule file

UNDOCUMENTED

E: Invalid shake bond type in molecule file

UNDOCUMENTED

E: Invalid shake angle type in molecule file

UNDOCUMENTED

W: Molecule attributes do not match system attributes

UNDOCUMENTED

E: Molecule topology type exceeds system topology type

UNDOCUMENTED

E: Molecule toplogy/atom exceeds system topology/atom

UNDOCUMENTED

W: Molecule has bond topology but no special bond settings

UNDOCUMENTED

E: Cannot open molecule file %s

UNDOCUMENTED

*/
