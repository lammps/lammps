/********************************
*
* Header file for Msi2LMP2 conversion program.
*
* Msi2lmp3
*
* This is the header file for the third version of a program
* that generates a LAMMPS data file based on the information
* in an MSI car file (atom coordinates) and mdf file (molecular
* topology). A key part of the program looks up forcefield parameters
* from an MSI frc file.
*
* The first version was written by Steve Lustig at Dupont, but
* required using Discover to derive internal coordinates and
* forcefield parameters
*
* The second version was written by Michael Peachey while an
* intern in the Cray Chemistry Applications Group managed
* by John Carpenter. This version derived internal coordinates
* from the mdf file and looked up parameters in the frc file
* thus eliminating the need for Discover.
*
* The third version was written by John Carpenter to optimize
* the performance of the program for large molecular systems
* (the original code for derving atom numbers was quadratic in time)
* and to make the program fully dynamic. The second version used
* fixed dimension arrays for the internal coordinates.
*
* John Carpenter can be contacted by sending email to
* jec374@earthlink.net
*
* November 2000
*
*/

# include <string.h>
# include <stddef.h>
# include <stdio.h>
# include <stdlib.h>
# include <math.h>

#ifdef   MAIN
#define  _EX
#define  _ARG(arg)    = (arg)
#else
#define  _EX extern
#define  _ARG(arg)
#endif

#define MAX_ATOM_TYPES 100
#define MAX_BOND_TYPES 200
#define MAX_ANGLE_TYPES 300
#define MAX_DIHEDRAL_TYPES 400
#define MAX_OOP_TYPES 400
#define MAX_ANGLEANGLE_TYPES 400
#define MAX_TYPES 12000
#define MAX_MEMBERS  5
#define MAX_LINE_LENGTH   200
#define MAX_PARAMS        8
#define PI_180  0.01745329252
#define MAX_CONNECTIONS 6
#define MAX_STRING 50



struct ResidueList {
  int start;
  int end;
  char name[16];
};

struct MoleculeList {
  int start;
  int end;
  int no_residues;
  struct ResidueList *residue;
};

/* Internal coodinate Lists */

struct BondList {
  int type;
  int members[2];
};

struct AngleList {
  int type;
  int members[3];
};

struct DihedralList {
  int type;
  int members[4];
};

struct OOPList {
  int type;
  int members[4];
};

struct AngleAngleList {
  int type;
  int members[4];
};

/* Internal coodinate Types Lists */


struct AtomTypeList
{
  char potential[5];
  float mass;
  double params[2];
  int no_connect;
};

struct BondTypeList {
  int types[2];
  double params[4];
};

struct AngleTypeList {
  int types[3];
  double params[4];
  double bondangle_cross_term[4];
  double bondbond_cross_term[3];
};

struct DihedralTypeList {
  int types[4];
  double params[6];
  double endbonddihedral_cross_term[8];
  double midbonddihedral_cross_term[4];
  double angledihedral_cross_term[8];
  double angleangledihedral_cross_term[3];
  double bond13_cross_term[3];
};

struct OOPTypeList {
  int types[4];
  double params[3];
  double angleangle_params[6];
};

struct AngleAngleTypeList {
  int types[4];
  double params[6];
};

/* ---------------------------------------------- */

struct Atom {
  int   molecule;      /* molecule id */
  int   no;            /* atom id */
  char  name[10];      /* atom name */
  double x[3];         /* position vector */
  char  potential[5];  /* atom potential type */
  char  element[2];    /* atom element */
  float q;             /* charge */
  char  residue_string[16]; /* residue string */
  int  no_connect;	/* number of connections to atom */
  char connections[MAX_CONNECTIONS][MAX_STRING];  /* long form, connection name*/
  double bond_order[6];
  int conn_no[6];	/* Atom number to which atom is connected */
  int type;
};


_EX  char   rootname[20];
_EX  char   path[20];
_EX  double pbc[9];
_EX  int    periodic   _ARG( 1 ); /* 0= nonperiodic 1= 3-D periodic */
// Added triclinic flag for non-orthogonal boxes Oct 5, 2010 SLTM
_EX  int TriclinicFlag; // 1 for non-orthoganal boxes, 0 for orthogonal boxes
_EX  int    forcefield _ARG( 0 ); /* 0= ClassI      1= ClassII */
_EX  int    pflag;
_EX  int    *no_atoms;
_EX  int    no_molecules;
_EX  int    replicate[3];
_EX  int    total_no_atoms;
_EX  int    total_no_bonds;
_EX  int    total_no_angles;
_EX  int    total_no_dihedrals;
_EX  int    total_no_angle_angles;
_EX  int    total_no_oops;
_EX  int    no_atom_types;
_EX  int    no_bond_types;
_EX  int    no_angle_types;
_EX  int    no_dihedral_types;
_EX  int    no_oop_types;
_EX  int    no_angleangle_types;
_EX  char   FrcFileName[MAX_LINE_LENGTH];
_EX  FILE   *CarF;
_EX  FILE   *FrcF;
_EX  FILE   *PrmF;
_EX  FILE   *MdfF;
_EX  FILE   *RptF;
_EX  struct Atom *atoms;
_EX  struct MoleculeList *molecule;
_EX  struct BondList *bonds;
_EX  struct AngleList *angles;
_EX  struct DihedralList *dihedrals;
_EX  struct OOPList *oops;
_EX  struct AngleAngleList *angleangles;
_EX  struct AtomTypeList *atomtypes;
_EX  struct BondTypeList *bondtypes;
_EX  struct AngleTypeList *angletypes;
_EX  struct DihedralTypeList *dihedraltypes;
_EX  struct OOPTypeList *ooptypes;
_EX  struct AngleAngleTypeList *angleangletypes;
#undef _EX
#undef _ARG

