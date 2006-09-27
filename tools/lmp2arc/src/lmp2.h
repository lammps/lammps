# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <stddef.h>
# include <math.h>

#define MAX_LINE_LENGTH   200
#define MAX_POS_FILES     50

#ifdef MAIN
#define _EX
#else
#define _EX extern  
#endif

struct Sys {
  int periodic;
  int ntypes;
  int natoms;
  int no_molecules;
  int nbonds;
  float celldim[3];
  float *masses;
  struct Mol *molinfo;
  struct Atom *atoms;
  int *bondindex;
};
  
struct Mol
{
  int start;
  int end;
};

struct Atom		/* atom information in .car file */
{
  int   molecule;	/* molecule id */
  float q;              /* charge */
  double xyz[3];         /* position vector */
  char  potential[5];	/* atom potential type */
  char  element[2];     /* atom element */
  char res_name[8];	/* residue name */
  char res_num[8];      /* residue numer */
  char  name[10];       /* atom name */
};


struct Boundary
{
  double low[3];
  double hi[3];
  double size[3];
};

struct NewAtomCoordinates
{
  int type;
  double fract[3];
  int truef[3];
};

struct Colors
{
  float rgb[3];
};

_EX    int    trueflag;           /* 0=no_true_flags; 1=true_flags */
_EX    int    move_molecules;     /* 0=don't move; 1=move */
_EX    int    nskip;              /* number of steps to skip in pos file */
_EX    int    npico;              /* number of steps per picosecond */


