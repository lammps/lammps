/*
*
*  msi2lmp.exe  V3.8
*
*   v3.6 KLA - Changes to output to either lammps 2001 (F90 version) or to
*              lammps 2005 (C++ version)
*
*   v3.4 JEC - a number of minor changes due to way newline and EOF are generated
*              on Materials Studio generated .car and .mdf files as well as odd
*              behavior out of newer Linux IO libraries. ReadMdfFile was restructured
*              in the process.
*
*   v3.1 JEC - changed IO interface to standard in/out, forcefield file
*              location can be indicated by environmental variable; added
*              printing options, consistency checks and forcefield
*              parameter versions sensitivity (highest one used)
*
*   v3.0 JEC - program substantially rewritten to reduce execution time
*              and be 98 % dynamic in memory use (still fixed limits on
*              number of parameter types for different internal coordinate
*              sets)
*
*   v2.0 MDP - got internal coordinate information from mdf file and
*              forcefield parameters from frc file thus eliminating
*              need for Discover
*
*   V1.0 SL  - original version. Used .car file and internal coordinate
*              information from Discover to produce LAMMPS data file.
*
*  This program uses the .car and .mdf files from MSI/Biosyms's INSIGHT
*  program to produce a LAMMPS data file.
*
*  The program is started by supplying information at the command prompt
* according to the usage described below.
*
*  USAGE: msi2lmp3 ROOTNAME {-print #} {-class #} {-frc FRC_FILE} {-ignore} {-nocenter}
*
*  -- msi2lmp3 is the name of the executable
*  -- ROOTNAME is the base name of the .car and .mdf files
*  -- all opther flags are optional and can be abbreviated (e.g. -p instead of -print)
*
*  -- -print
*        # is the print level:  0  - silent except for errors
*                               1  - minimal (default)
*                               2  - more verbose
*                               3  - even more verbose
*  -- -class
*        # is the class of forcefield to use (I  = Class I e.g., CVFF)
*                                            (II = Class II e.g., CFFx )
*     default is -class I
*
*  -- -ignore   - tells msi2lmp to ignore warnings and errors and keep going
*
*  -- -nocenter - tells msi2lmp to not center the box around the (geometrical)
*                 center of the atoms, but around the origin
*
*  -- -frc      - specifies name of the forcefield file (e.g., cff91)
*
*     If the name includes a hard wired directory (i.e., if the name
*     starts with . or /), then the name is used alone. Otherwise,
*     the program looks for the forcefield file in $BIOSYM_LIBRARY.
*     If $BIOSYM_LIBRARY is not set, then the current directory is
*     used.
*
*     If the file name does not include a dot after the first
*     character, then .frc is appended to the name.
*
*     For example,  -frc cvff (assumes cvff.frc is in $BIOSYM_LIBRARY
*                              or .)
*
*                   -frc cff/cff91 (assumes cff91.frc is in
*                                   $BIOSYM_LIBRARY/cff or ./cff)
*
*                   -frc /usr/local/biosym/forcefields/cff95 (absolute
*                                                             location)
*
*     By default, the program uses $BIOSYM_LIBRARY/cvff.frc
*
*  -- output is written to a file called ROOTNAME.data
*
*
****************************************************************
*
* Msi2lmp3
*
* This is the third version of a program that generates a LAMMPS
* data file based on the information in a MSI car file (atom
* coordinates) and mdf file (molecular topology). A key part of
* the program looks up forcefield parameters from an MSI frc file.
*
* The first version was written by Steve Lustig at Dupont, but
* required using Discover to derive internal coordinates and
* forcefield parameters
*
* The second version was written by Michael Peachey while an
* in intern in the Cray Chemistry Applications Group managed
* by John Carpenter. This version derived internal coordinates
* from the mdf file and looked up parameters in the frc file
* thus eliminating the need for Discover.
*
* The third version was written by John Carpenter to optimize
* the performance of the program for large molecular systems
* (the original  code for deriving atom numbers was quadratic in time)
* and to make the program fully dynamic. The second version used
* fixed dimension arrays for the internal coordinates.
*
* John Carpenter can be contacted by sending email to
* jec374@earthlink.net
*
* November 2000
*/

#include "msi2lmp.h"

#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/* global variables */

char  *rootname;
double pbc[6];
double box[3][3];
int    periodic = 1;
int    TriclinicFlag = 0;
int    forcefield = 0;
int    centerflag = 1;

int    pflag;
int    iflag;
int   *no_atoms;
int    no_molecules;
int    replicate[3];
int    total_no_atoms;
int    total_no_bonds;
int    total_no_angles;
int    total_no_dihedrals;
int    total_no_angle_angles;
int    total_no_oops;
int    no_atom_types;
int    no_bond_types;
int    no_angle_types;
int    no_dihedral_types;
int    no_oop_types;
int    no_angleangle_types;
char   *FrcFileName;
FILE   *CarF;
FILE   *FrcF;
FILE   *PrmF;
FILE   *MdfF;
FILE   *RptF;

struct Atom *atoms;
struct MoleculeList *molecule;
struct BondList *bonds;
struct AngleList *angles;
struct DihedralList *dihedrals;
struct OOPList *oops;
struct AngleAngleList *angleangles;
struct AtomTypeList *atomtypes;
struct BondTypeList *bondtypes;
struct AngleTypeList *angletypes;
struct DihedralTypeList *dihedraltypes;
struct OOPTypeList *ooptypes;
struct AngleAngleTypeList *angleangletypes;


static int check_arg(char **arg, const char *flag, int num, int argc)
{
  if (num >= argc) {
    printf("Missing argument to \"%s\" flag\n",flag);
    return 1;
  }
  if (arg[num][0] == '-') {
    printf("Incorrect argument to \"%s\" flag: %s\n",flag,arg[num]);
    return 1;
  }
  return 0;
}

int main (int argc, char *argv[])
{
  int n,i,found_sep;
  const char *frc_dir_name = NULL;
  const char *frc_file_name = NULL;

  pflag = 1;
  iflag = 0;
  forcefield = 1;

  frc_dir_name = getenv("BIOSYM_LIBRARY");

  if (argc < 2) {
    printf("usage: %s <rootname> [-class <I|1|II|2>] [-frc <path to frc file>] [-print #] [-ignore] [-nocenter]\n",argv[0]);
    return 1;
  } else { /* rootname was supplied as first argument, copy to rootname */
    int len = strlen(argv[1]) + 1;
    rootname = (char *)malloc(len);
    strcpy(rootname,argv[1]);
  }

  n = 2;
  while (n < argc) {
    if (strncmp(argv[n],"-c",2) == 0) {
      n++;
      if (check_arg(argv,"-class",n,argc))
        return 2;
      if ((strcmp(argv[n],"I") == 0) || (strcmp(argv[n],"1") == 0)) {
        forcefield = 1;
      } else if ((strcmp(argv[n],"II") == 0) || (strcmp(argv[n],"2") == 0)) {
        forcefield = 2;
      } else {
        printf("Unrecognized Forcefield class: %s\n",argv[n]);
        return 3;
      }
    } else if (strncmp(argv[n],"-f",2) == 0) {
      n++;
      if (check_arg(argv,"-frc",n,argc))
        return 4;
      frc_file_name = argv[n];
    } else if (strncmp(argv[n],"-i",2) == 0 ) {
      iflag = 1;
    } else if (strncmp(argv[n],"-n",2) == 0 ) {
      centerflag = 0;
    } else if (strncmp(argv[n],"-p",2) == 0) {
      n++;
      if (check_arg(argv,"-print",n,argc))
        return 5;
      pflag = atoi(argv[n]);
    } else {
      printf("Unrecognized option: %s\n",argv[n]);
      return 6;
    }
    n++;
  }

  /* set defaults, if nothing else was given */
  if (frc_dir_name == NULL)
#if (_WIN32)
    frc_dir_name = "..\\biosym_frc_files";
#else
  frc_dir_name = "../biosym_frc_files";
#endif
  if (frc_file_name == NULL)
    frc_file_name = "cvff.frc";

  found_sep=0;
#ifdef _WIN32
  if (isalpha(frc_file_name[0]) && (frc_file_name[1] == ':'))
    found_sep=1; /* windows drive letter => full path. */
#endif

  n = strlen(frc_file_name);
  for (i=0; i < n; ++i) {
#ifdef _WIN32
    if ((frc_file_name[i] == '/') || (frc_file_name[i] == '\\'))
      found_sep=1+i;
#else
    if (frc_file_name[i] == '/')
      found_sep=1+i;
#endif
  }

  /* full pathname given */
  if (found_sep) {
    i = 0;
    /* need to append extension? */
    if ((n < 5) || (strcmp(frc_file_name+n-4,".frc") !=0))
      i=1;

    FrcFileName = (char *)malloc(n+1+i*4);
    strcpy(FrcFileName,frc_file_name);
    if (i) strcat(FrcFileName,".frc");
  } else {
    i = 0;
    /* need to append extension? */
    if ((n < 5) || (strcmp(frc_file_name+n-4,".frc") !=0))
      i=1;

    FrcFileName = (char *)malloc(n+2+i*4+strlen(frc_dir_name));
    strcpy(FrcFileName,frc_dir_name);
#ifdef _WIN32
    strcat(FrcFileName,"\\");
#else
    strcat(FrcFileName,"/");
#endif
    strcat(FrcFileName,frc_file_name);
    if (i) strcat(FrcFileName,".frc");
  }


  if (pflag > 0) {
    printf("\nRunning msi2lmp.....\n\n");
    printf(" Forcefield file name: %s\n",FrcFileName);
    printf(" Forcefield class: %d\n\n",forcefield);
  }

  if (((forcefield == 1) && (strstr(FrcFileName,"cff") != NULL)) ||
      ((forcefield == 2) &&
       ! ((strstr(FrcFileName,"cvff") == NULL)
          || (strstr(FrcFileName,"clayff") == NULL)
          || (strstr(FrcFileName,"compass") == NULL)))) {
    fprintf(stderr," WARNING - forcefield name and class appear to\n");
    fprintf(stderr,"           be inconsistent - Errors may result\n\n");
  }

  /* Read in .car file */
  ReadCarFile();

  /*Read in .mdf file */

  ReadMdfFile();

  /* Define bonds, angles, etc...*/

  if (pflag > 0)
    printf("\n Building internal coordinate lists \n");
  MakeLists();

  /* Read .frc file into memory */

  if (pflag > 0)
    printf("\n Reading forcefield file \n");
  ReadFrcFile();

  /* Get forcefield parameters */

  if (pflag > 0)
    printf("\n Get force field parameters for this system\n");
  GetParameters(forcefield);

  /* Do internal check of internal coordinate lists */
  if (pflag > 0)
    printf("\n Check parameters for internal consistency\n");
  CheckLists();

  /* Write out the final data */
  WriteDataFile(rootname,forcefield);

  free(rootname);
  if (pflag > 0)
    printf("\nNormal program termination\n");
  return 0;
}

