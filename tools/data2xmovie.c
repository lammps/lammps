/* data2xmovie tool

   read LAMMPS data file as input
   write a snapshot in XMOVIE format

   Syntax:  data2xmovie [options] < infile > outfile

     Options:

     -style atom_style

          use the LAMMPS atom style that corresponds to this file
	    e.g. atomic or bond or angle or full or eam or granular
          will be used for reading and writing files
          if not specified, atom_style = full

      -unmap

          unmap all input atom positions using input image flags
          image flags must be specified in infile
	  default is to leave atoms mapped to periodic box

     -inbox xlo xhi ylo yhi zlo zhi

          use these values for the output bounding box of the system
	    useful if are unmapping atoms into a larger domain
	  if not specified use the box bounds read in from infile
*/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct box {
  double xlo,xhi,xsize;
  double ylo,yhi,ysize;
  double zlo,zhi,zsize;
};

#define MAXLINE 1000

#define ATOMIC 1
#define BOND   2
#define ANGLE  3
#define FULL   4
#define EAM    5
#define GRANULAR 6

main(int argc, char *argv[])

{
  char line[1000];            /* strings for reading/parsing input file */
  char linetoken[1000];
  char *token;             
                             /* numbers of various quantities as read in */

  int natoms,nbonds,nangles,ndihedrals,nimpropers;
  int ntypes,nbondtypes,nangletypes,ndihedtypes,nimprotypes;

  int nmolecules;            /* total # of mols specified by mol tags */

  int atomstyle;             /* atomstyle for read/write of Atoms */
  int unmap;                 /* 0 = no unmapping, 1 = unmap from input box */
  int outbox;                /* flag for whether out box is explicitly set */

  struct box in,out;             /* input/output bounding box */

  char *gettoken(char *, int);    /* function defs */
  void skipline(int); 

  /* default input values */

  atomstyle = 0;
  unmap = 0;
  outbox = 0;
  
  /* read input options from command line
     should do more error checking for missing args */

  int i = 1;
  while (i < argc) {

    if (!strcmp(argv[i],"-style")) {
      if (strcmp(argv[i+1],"atomic") == 0) atomstyle = ATOMIC;
      else if (strcmp(argv[i+1],"bond") == 0) atomstyle = BOND;
      else if (strcmp(argv[i+1],"angle") == 0) atomstyle = ANGLE;
      else if (strcmp(argv[i+1],"full") == 0) atomstyle = FULL;
      else if (strcmp(argv[i+1],"eam") == 0) atomstyle = EAM;
      else if (strcmp(argv[i+1],"granular") == 0) atomstyle = GRANULAR;
      else {
	fprintf(stderr,"Error with command-line arg style\n");
	exit(1);
      }
      i += 2;
    } else if (!strcmp(argv[i],"-unmap")) {
      unmap = 1;
      i += 1;
    } else if (!strcmp(argv[i],"-inbox")) {
      outbox = 1;
      sscanf(argv[i+1],"%lg",&out.xlo);
      sscanf(argv[i+2],"%lg",&out.xhi);
      sscanf(argv[i+3],"%lg",&out.ylo);
      sscanf(argv[i+4],"%lg",&out.yhi);
      sscanf(argv[i+5],"%lg",&out.zlo);
      sscanf(argv[i+6],"%lg",&out.zhi);
      i += 7;
    } else {
      fprintf(stderr,"Syntax error: data2xmovie [options] < infile > outfile\n");
      exit(1);
    }
  }
  
  /* error checks */

  if (atomstyle == 0) {
    fprintf(stderr,"ERROR: must use -style to set atom style\n");
    exit(1);
  }

  /* defaults */

  natoms = nbonds = nangles = ndihedrals = nimpropers = 0;
  ntypes = nbondtypes = nangletypes = ndihedtypes = nimprotypes = 0;

  /* read header */

  fgets(line,MAXLINE,stdin);

  while (1) {
    fgets(line,MAXLINE,stdin);

    if (strspn(line," \t\n\r") == strlen(line)) continue;
    else if (strstr(line,"atoms")) sscanf(line,"%d",&natoms);
    else if (strstr(line,"bonds")) sscanf(line,"%d",&nbonds);
    else if (strstr(line,"angles")) sscanf(line,"%d",&nangles);
    else if (strstr(line,"dihedrals")) sscanf(line,"%d",&ndihedrals);
    else if (strstr(line,"impropers")) sscanf(line,"%d",&nimpropers);
    else if (strstr(line,"atom types")) sscanf(line,"%d",&ntypes);
    else if (strstr(line,"bond types")) sscanf(line,"%d",&nbondtypes);
    else if (strstr(line,"angle types")) sscanf(line,"%d",&nangletypes);
    else if (strstr(line,"dihedral types")) sscanf(line,"%d",&ndihedtypes);
    else if (strstr(line,"improper types")) sscanf(line,"%d",&nimprotypes);
    else if (strstr(line,"xlo xhi")) sscanf(line,"%lg %lg",&in.xlo,&in.xhi);
    else if (strstr(line,"ylo yhi")) sscanf(line,"%lg %lg",&in.ylo,&in.yhi);
    else if (strstr(line,"zlo zhi")) sscanf(line,"%lg %lg",&in.zlo,&in.zhi);
    else break;
  }

  /* compute input box size */

  in.xsize = in.xhi - in.xlo;
  in.ysize = in.yhi - in.ylo;
  in.zsize = in.zhi - in.zlo;

  /* write XMOVIE header */

  printf("ITEM: TIMESTEP\n");
  printf("%d\n",0);
  printf("ITEM: NUMBER OF ATOMS\n");
  printf("%d\n",natoms);
  printf("ITEM: BOX BOUNDS\n");
  if (outbox) {
    printf("%g %g\n",out.xlo,out.xhi);
    printf("%g %g\n",out.ylo,out.yhi);
    printf("%g %g\n",out.zlo,out.zhi);
  } else {
    printf("%g %g\n",in.xlo,in.xhi);
    printf("%g %g\n",in.ylo,in.yhi);
    printf("%g %g\n",in.zlo,in.zhi);
  }

  /* read identifier strings one by one in free-form part of data file */

  token = gettoken(line,1);

  while (token) {

    /* read atoms */

    if (!strcmp(token,"Atoms")) {

      printf("ITEM: ATOMS\n");

      int tag,type,molecule,imagex,imagey,imagez;
      double x,y,z,radius,density,q;

      for (i = 0; i < natoms; i++) {
	fgets(line,MAXLINE,stdin);
	if (unmap) {
	  if (atomstyle == ATOMIC) 
	    sscanf(line,"%d %d %lg %lg %lg %d %d %d",
		   &tag,&type,&x,&y,&z,
		   &imagex,&imagey,&imagez);
	  else if (atomstyle == BOND) 
	    sscanf(line,"%d %d %d %lg %lg %lg %d %d %d",
		   &tag,&molecule,&type,&x,&y,&z,
		   &imagex,&imagey,&imagez);
	  else if (atomstyle == ANGLE) 
	    sscanf(line,"%d %d %d %lg %lg %lg %d %d %d",
		   &tag,&molecule,&type,&x,&y,&z,
		   &imagex,&imagey,&imagez);
	  else if (atomstyle == FULL) 
	    sscanf(line,"%d %d %d %lg %lg %lg %lg %d %d %d",
		   &tag,&molecule,&type,&q,&x,&y,&z,
		   &imagex,&imagey,&imagez);
	  else if (atomstyle == EAM) 
	    sscanf(line,"%d %d %lg %lg %lg %d %d %d",
		   &tag,&type,&x,&y,&z,
		   &imagex,&imagey,&imagez);
	  else if (atomstyle == GRANULAR) 
	    sscanf(line,"%d %d %lg %lg %lg %lg %lg %d %d %d",
		   &tag,&type,&radius,&density,&x,&y,&z,
		   &imagex,&imagey,&imagez);
	} else {
	  if (atomstyle == ATOMIC) 
	    sscanf(line,"%d %d %lg %lg %lg",
		   &tag,&type,&x,&y,&z);
	  else if (atomstyle == BOND) 
	    sscanf(line,"%d %d %d %lg %lg %lg",
		   &tag,&molecule,&type,&x,&y,&z);
	  else if (atomstyle == ANGLE) 
	    sscanf(line,"%d %d %d %lg %lg %lg",
		   &tag,&molecule,&type,&x,&y,&z);
	  else if (atomstyle == FULL) 
	    sscanf(line,"%d %d %d %lg %lg %lg %lg",
		   &tag,&molecule,&type,&q,&x,&y,&z);
	  else if (atomstyle == EAM) 
	    sscanf(line,"%d %d %lg %lg %lg",
		   &tag,&type,&x,&y,&z);
	  else if (atomstyle == GRANULAR) 
	    sscanf(line,"%d %d %lg %lg %lg %lg %lg",
		   &tag,&type,&radius,&density,&x,&y,&z);
	  imagez = imagey = imagex = 0;
	}

	/* unmap atom position if requested */

	if (unmap) {
	  x = x + imagex*in.xsize;
	  y = y + imagey*in.ysize;
	  z = z + imagez*in.zsize;
	}

	printf("%d %d %g %g %g\n",tag,type,x,y,z);
      }
    }
      
    /* read bonds and replicate */

    else if (!strcmp(token,"Bonds")) {

      printf("ITEM: BONDS\n");

      int n,btype,bond1,bond2;

      for (i = 0; i < nbonds; i++) {
	fgets(line,MAXLINE,stdin);
	sscanf(line,"%d %d %d %d",&n,&btype,&bond1,&bond2);
	printf("%d %d %d\n",btype,bond1,bond2);
      }
    }
    
	/* non-replicated sections - just skip lines */

    else if (!strcmp(token,"Velocities"))
      skipline(natoms);
    else if (!strcmp(token,"Angles"))
      skipline(nangles);
    else if (!strcmp(token,"Dihedrals"))
      skipline(ndihedrals);
    else if (!strcmp(token,"Impropers"))
      skipline(nimpropers);
    else if (!strcmp(token,"Masses"))
      skipline(ntypes);
    else if (!strcmp(token,"Dipoles"))
      skipline(ntypes);
    else if (!strcmp(token,"Pair Coeffs"))
      skipline(ntypes);
    else if (!strcmp(token,"Bond Coeffs"))
      skipline(nbondtypes);
    else if (!strcmp(token,"Angle Coeffs"))
      skipline(nangletypes);
    else if (!strcmp(token,"Dihedral Coeffs"))
      skipline(ndihedtypes);
    else if (!strcmp(token,"Improper Coeffs"))
      skipline(nimprotypes);
    else if (!strcmp(token,"BondBond Coeffs"))
      skipline(nangletypes);
    else if (!strcmp(token,"BondAngle Coeffs"))
      skipline(nangletypes);
    else if (!strcmp(token,"MiddleBondTorsion Coeffs"))
      skipline(ndihedtypes);
    else if (!strcmp(token,"EndBondTorsion Coeffs"))
      skipline(ndihedtypes);
    else if (!strcmp(token,"AngleTorsion Coeffs"))
      skipline(ndihedtypes);
    else if (!strcmp(token,"AngleAngleTorsion Coeffs"))
      skipline(ndihedtypes);
    else if (!strcmp(token,"BondBond13 Coeffs"))
      skipline(ndihedtypes);
    else if (!strcmp(token,"AngleAngle Coeffs"))
      skipline(nimprotypes);
    else {
      fprintf(stderr,
	      "Error in input data file - unknown identifier %s\n",token);
      exit(1);
    }

    token = gettoken(line,0);
  }
  
}

/* ------------------------------------------------------------------- */

/* return a LAMMPS keyword from data file
   if first is 1, non-blank line with token is in line
   else read until find a non-blank line
   keyword is all text on line w/out leading & trailing white space
   read one additional line after non-blank line (assumed blank)
   return ptr to keyword
   return NULL if end of file */

char *gettoken(char *line, int first)
{
  char buffer[MAXLINE];

  /* read upto non-blank line plus 1 following line
     eof is set to 1 if any read hits end-of-file */

  int eof = 0;

  if (!first) if (fgets(line,MAXLINE,stdin) == NULL) eof = 1;
  while (eof == 0 && strspn(line," \t\n\r") == strlen(line))
    if (fgets(line,MAXLINE,stdin) == NULL) eof = 1;
  if (fgets(buffer,MAXLINE,stdin) == NULL) eof = 1;

  /* if eof, return NULL */

  if (eof) return NULL;

  /* bracket non-whitespace portion of line */

  int start = strspn(line," \t\n\r");
  int stop = strlen(line) - 1;
  while (line[stop] == ' ' || line[stop] == '\t' 
	 || line[stop] == '\n' || line[stop] == '\r') stop--;
  line[stop+1] = '\0';

  /* return ptr to keyword */

  return &line[start];
}

/* read n lines and ignore */

void skipline(int n)
{
  char line[1000];

  while (n) {
    fgets(line,MAXLINE,stdin);
    n--;
  }
}
