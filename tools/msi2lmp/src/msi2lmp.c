/*
*
*  msi2lmp.exe
*
*   v3.9.9 AK- Teach msi2lmp to not generate dihedrals with identical 1-4 atoms
*
*   v3.9.8 AK- Improved whitespace handling in parsing topology and force
*              field files to avoid bogus warnings about type name truncation
*
*   v3.9.7 AK- Add check to enforce that Class1/OPLS-AA use A-B parameter
*              conventions in force field file and Class2 us r-eps conventions
*
*   v3.9.6 AK- Refactoring of MDF file parser with more consistent
*              handling of compile time constants MAX_NAME and MAX_STRING
*
*   v3.9.5 AK- Add TopoTools style force field parameter type hints
*
*   v3.9.4 AK- Make force field style hints optional with a flag
*
*   v3.9.3 AK- Bugfix for triclinic cells.
*
*   v3.9.2 AK- Support for writing out force field style hints
*
*   v3.9.1 AK- Bugfix for Class2. Free allocated memory. Print version number.
*
*   v3.9 AK  - Rudimentary support for OPLS-AA
*
*   v3.8 AK  - Some refactoring and cleanup of global variables
*            - Bugfixes for argument parsing and improper definitions
*            - improved handling of box dimensions and image flags
*            - port to compiling on windows using MinGW
*            - more consistent print level handling
*            - more consistent handling of missing parameters
*            - Added a regression test script with examples.
*
*   V3.7 STM - Added support for triclinic cells
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
*  USAGE: msi2lmp3 ROOTNAME {-print #} {-class #} {-frc FRC_FILE} {-ignore} {-nocenter} {-oldstyle}
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
*        # is the class of forcefield to use (I  or 1 = Class I e.g., CVFF, clayff)
*                                            (II or 2 = Class II e.g., CFFx, COMPASS)
*                                            (O  or 0 = OPLS-AA)
*     default is -class I
*
*  -- -ignore   - tells msi2lmp to ignore warnings and errors and keep going
*
*  -- -nocenter - tells msi2lmp to not center the box around the (geometrical)
*                 center of the atoms, but around the origin
*
*  -- -oldstyle - tells msi2lmp to write out a data file without style hints
*                 (to be compatible with older LAMMPS versions)
*
*  -- -shift    - tells msi2lmp to shift the entire system (box and coordinates)
*                 by a vector (default: 0.0 0.0 0.0)
*
*  -- -frc      - specifies name of the forcefield file (e.g., cff91)
*
*     If the name includes a hard wired directory (i.e., if the name
*     starts with . or /), then the name is used alone. Otherwise,
*     the program looks for the forcefield file in $MSI2LMP_LIBRARY.
*     If $MSI2LMP_LIBRARY is not set, then the current directory is
*     used.
*
*     If the file name does not include a dot after the first
*     character, then .frc is appended to the name.
*
*     For example,  -frc cvff (assumes cvff.frc is in $MSI2LMP_LIBRARY
*                              or .)
*
*                   -frc cff/cff91 (assumes cff91.frc is in
*                                   $MSI2LMP_LIBRARY/cff or ./cff)
*
*                   -frc /usr/local/forcefields/cff95 (absolute
*                                                             location)
*
*     By default, the program uses $MSI2LMP_LIBRARY/cvff.frc
*
*  -- output is written to a file called ROOTNAME.data
*
*
****************************************************************
*
* msi2lmp
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
* November 2000
*/

#include "msi2lmp.h"

#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <ctype.h>
#endif

/* global variables */

char  *rootname;
double pbc[6];
double box[3][3];
double shift[3];
int    periodic = 1;
int    TriclinicFlag = 0;
int    forcefield = 0;
int    centerflag = 1;
int    hintflag = 1;
int    ljtypeflag = 0;

int    pflag;
int    iflag;
int   *no_atoms;
int    no_molecules;
int    replicate[3];
int    total_no_atoms = 0;
int    total_no_bonds = 0;
int    total_no_angles = 0;
int    total_no_dihedrals = 0;
int    total_no_angle_angles = 0;
int    total_no_oops = 0;
int    no_atom_types = 0;
int    no_bond_types = 0;
int    no_angle_types = 0;
int    no_dihedral_types = 0;
int    no_oop_types = 0;
int    no_angleangle_types = 0;
char   *FrcFileName = NULL;
FILE   *CarF = NULL;
FILE   *FrcF = NULL;
FILE   *PrmF = NULL;
FILE   *MdfF = NULL;
FILE   *RptF = NULL;

struct Atom *atoms = NULL;
struct MoleculeList *molecule = NULL;
struct BondList *bonds = NULL;
struct AngleList *angles = NULL;
struct DihedralList *dihedrals = NULL;
struct OOPList *oops = NULL;
struct AngleAngleList *angleangles = NULL;
struct AtomTypeList *atomtypes = NULL;
struct BondTypeList *bondtypes = NULL;
struct AngleTypeList *angletypes = NULL;
struct DihedralTypeList *dihedraltypes = NULL;
struct OOPTypeList *ooptypes = NULL;
struct AngleAngleTypeList *angleangletypes = NULL;

void condexit(int val)
{
    if (iflag == 0) exit(val);
}

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
  forcefield = FF_TYPE_CLASS1 | FF_TYPE_COMMON;
  shift[0] = shift[1] = shift[2] = 0.0;

  frc_dir_name = getenv("MSI2LMP_LIBRARY");

  if (argc < 2) {
    printf("usage: %s <rootname> [-class <I|1|II|2>] [-frc <path to frc file>] [-print #] [-ignore] [-nocenter] [-oldstyle]\n",argv[0]);
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
        forcefield = FF_TYPE_CLASS1 | FF_TYPE_COMMON;
      } else if ((strcmp(argv[n],"II") == 0) || (strcmp(argv[n],"2") == 0)) {
        forcefield = FF_TYPE_CLASS2 | FF_TYPE_COMMON;
      } else if ((strcmp(argv[n],"O") == 0) || (strcmp(argv[n],"0") == 0)) {
        forcefield = FF_TYPE_OPLSAA | FF_TYPE_COMMON;
      } else {
        printf("Unrecognized Forcefield class: %s\n",argv[n]);
        return 3;
      }
    } else if (strncmp(argv[n],"-f",2) == 0) {
      n++;
      if (check_arg(argv,"-frc",n,argc))
        return 4;
      frc_file_name = argv[n];
    } else if (strncmp(argv[n],"-s",2) == 0) {
      if (n+3 > argc) {
        printf("Missing argument(s) to \"-shift\" flag\n");
        return 1;
      }
      shift[0] = atof(argv[++n]);
      shift[1] = atof(argv[++n]);
      shift[2] = atof(argv[++n]);
    } else if (strncmp(argv[n],"-i",2) == 0 ) {
      iflag = 1;
    } else if (strncmp(argv[n],"-n",2) == 0 ) {
      centerflag = 0;
    } else if (strncmp(argv[n],"-o",2) == 0 ) {
      hintflag = 0;
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
  if (frc_dir_name == NULL) {
#if (_WIN32)
    frc_dir_name = "..\\frc_files";
#else
    frc_dir_name = "../frc_files";
#endif
  }

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
    puts("\nRunning msi2lmp " MSI2LMP_VERSION "\n");
    if (forcefield & FF_TYPE_CLASS1) puts(" Forcefield: Class I");
    if (forcefield & FF_TYPE_CLASS2) puts(" Forcefield: Class II");
    if (forcefield & FF_TYPE_OPLSAA) puts(" Forcefield: OPLS-AA");
    printf(" Forcefield file name: %s\n",FrcFileName);
    if (centerflag) puts(" Output is recentered around geometrical center");
    if (hintflag) puts(" Output contains style flag hints");
    else puts(" Style flag hints disabled");
    printf(" System translated by: %g %g %g\n",shift[0],shift[1],shift[2]);
  }

  n = 0;
  if (forcefield & FF_TYPE_CLASS1) {
    if (strstr(FrcFileName,"cvff") != NULL) ++n;
    if (strstr(FrcFileName,"clayff") != NULL) ++n;
  } else if (forcefield & FF_TYPE_OPLSAA) {
    if (strstr(FrcFileName,"oplsaa") != NULL) ++n;
  } else if (forcefield & FF_TYPE_CLASS2) {
    if (strstr(FrcFileName,"pcff") != NULL) ++n;
    if (strstr(FrcFileName,"cff91") != NULL) ++n;
    if (strstr(FrcFileName,"compass") != NULL) ++n;
  }

  if (n == 0) {
    if (iflag > 0) fputs(" WARNING",stderr);
    else           fputs(" Error  ",stderr);

    fputs("- forcefield name and class appear to be inconsistent\n\n",stderr);
    if (iflag == 0) return 7;
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
  GetParameters();

  /* Do internal check of internal coordinate lists */
  if (pflag > 0)
    printf("\n Check parameters for internal consistency\n");
  CheckLists();

  /* Write out the final data */
  WriteDataFile(rootname);

  /* free up memory to detect possible memory corruption */
  free(rootname);
  free(FrcFileName);
  ClearFrcData();

  for (n=0; n < no_molecules; n++) {
    free(molecule[n].residue);
  }

  free(no_atoms);
  free(molecule);
  free(atoms);
  free(atomtypes);
  if (bonds) free(bonds);
  if (bondtypes) free(bondtypes);
  if (angles) free(angles);
  if (angletypes) free(angletypes);
  if (dihedrals) free(dihedrals);
  if (dihedraltypes) free(dihedraltypes);
  if (oops) free(oops);
  if (ooptypes) free(ooptypes);
  if (angleangles) free(angleangles);
  if (angleangletypes) free(angleangletypes);

  if (pflag > 0)
    printf("\nNormal program termination\n");
  return 0;
}
