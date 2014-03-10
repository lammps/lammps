/*
*  This function opens the .car file and extracts coordinate information
*  into the atoms Atom structure
*/

#include "msi2lmp.h"

#include <stdlib.h>
#include <string.h>
#include <math.h>


/* ----------------------------------------------------------------------
   set global box params
   assumes boxlo/hi and triclinic tilts are already set
------------------------------------------------------------------------- */

void set_box(double box[3][3], double *h, double *h_inv)
{
  h[0] = box[1][0] - box[0][0];
  h[1] = box[1][1] - box[0][1];
  h[2] = box[1][2] - box[0][2];

  h_inv[0] = 1.0/h[0];
  h_inv[1] = 1.0/h[1];
  h_inv[2] = 1.0/h[2];

  h[3] = box[2][2];
  h[4] = box[2][1];
  h[5] = box[2][0];
  h_inv[3] = -h[3] / (h[1]*h[2]);
  h_inv[4] = (h[3]*h[5] - h[1]*h[4]) / (h[0]*h[1]*h[2]);
  h_inv[5] = -h[5] / (h[0]*h[1]);
}


/* ----------------------------------------------------------------------
   convert triclinic 0-1 lamda coords to box coords for one atom
   x = H lamda + x0;
   lamda and x can point to same 3-vector
------------------------------------------------------------------------- */

void lamda2x(double *lamda, double *x, double *h, double *boxlo)
{
  x[0] = h[0]*lamda[0] + h[5]*lamda[1] + h[4]*lamda[2] + boxlo[0];
  x[1] = h[1]*lamda[1] + h[3]*lamda[2] + boxlo[1];
  x[2] = h[2]*lamda[2] + boxlo[2];
}

/* ----------------------------------------------------------------------
   convert box coords to triclinic 0-1 lamda coords for one atom
   lamda = H^-1 (x - x0)
   x and lamda can point to same 3-vector
------------------------------------------------------------------------- */

void x2lamda(double *x, double *lamda, double *h_inv, double *boxlo)
{
  double delta[3];
  delta[0] = x[0] - boxlo[0];
  delta[1] = x[1] - boxlo[1];
  delta[2] = x[2] - boxlo[2];

  lamda[0] = h_inv[0]*delta[0] + h_inv[5]*delta[1] + h_inv[4]*delta[2];
  lamda[1] = h_inv[1]*delta[1] + h_inv[3]*delta[2];
  lamda[2] = h_inv[2]*delta[2];
}


void ReadCarFile(void)
{
  char line[MAX_LINE_LENGTH];  /* Stores lines as they are read in */
  int k,m,n;                   /* counters */
  int skip;                    /* lines to skip at beginning of file */
  double lowest, highest;      /* temp coordinate finding variables */
  double total_q;
  double sq_c;
  double cos_alpha;  /* Added by SLTM Sept 13, 2010 */
  double cos_gamma;
  double sin_gamma;
  double cos_beta;
  double sin_beta;
  double A, B, C;
  double center[3];
  double hmat[6];
  double hinv[6];
  double lamda[3];

  /* Open .car file for reading */

  sprintf(line,"%s.car",rootname);
  if (pflag > 0) printf(" Reading car file: %s\n",line);
  if( (CarF = fopen(line,"r")) == NULL ) {
    printf("Cannot open %s\n",line);
    exit(33);
  }

  /* Determine Number of molecules & atoms */

  rewind(CarF);
  no_molecules = -1; /* Set to -1 because counter will be incremented an
                        extra time at the end of the file */

  fgets(line,MAX_LINE_LENGTH,CarF); /* Read header line */

  /* Check for periodicity, if present, read cell constants */

  if( strncmp(fgets(line,MAX_LINE_LENGTH,CarF),"PBC=ON",6) == 0) {
    periodic = 1;
    skip = 5; /* Data starts 5 lines from beginning of file */
    fgets(line,MAX_LINE_LENGTH,CarF); /* Comment line */
    fgets(line,MAX_LINE_LENGTH,CarF); /* Date stamp */
    fscanf(CarF,"%*s %lf %lf %lf %lf %lf %lf %*s",
           &pbc[0],&pbc[1],&pbc[2],&pbc[3],&pbc[4],&pbc[5]);

    /* Added triclinic flag for non-orthogonal boxes Oct 5, 2010 SLTM */
    if(pbc[3] != 90.0 || pbc[4] != 90.0 || pbc[5] != 90.0) {
      TriclinicFlag = 1;
    } else TriclinicFlag = 0;
  } else {
    periodic = 0;
    skip = 4;
    if (pflag > 1) {
      printf("   %s is not a periodic system\n", rootname);
      printf("   Assigning cell parameters based on coordinates\n");
    }
    fgets(line,MAX_LINE_LENGTH, CarF); /* Comment line */
    fgets(line,MAX_LINE_LENGTH, CarF); /* Date Stamp */
  }

  /* First pass through file -- Count molecules */

  while(fgets(line,MAX_LINE_LENGTH,CarF) != NULL )
    if( strncmp(line,"end",3) == 0 )
      no_molecules++;

  /* Allocate space to keep track of the number of atoms within a molecule */

  no_atoms = (int *) calloc(no_molecules,sizeof(int));
  if ( no_atoms == NULL ) {
    printf("Could not allocate memory for no_atoms\n");
    exit(32);
  }

  /* Second pass through file -- Count atoms */

  rewind(CarF);
  for(n=0; n < skip; n++)               /* Skip beginning lines */
    fgets(line,MAX_LINE_LENGTH,CarF);

  for(n=0; n < no_molecules; n++)
    while( strncmp(fgets(line,MAX_LINE_LENGTH,CarF),"end",3) )
      no_atoms[n]++;

  for( total_no_atoms=0, n=0; n < no_molecules; n++ )
    total_no_atoms += no_atoms[n];

  molecule = (struct MoleculeList *) calloc(no_molecules,
                                            sizeof(struct MoleculeList));
  if (molecule == NULL) {
    printf("Unable to allocate memory for molecule structure\n");
    exit(32);
  }
  molecule[0].start = 0;
  molecule[0].end = no_atoms[0];
  for (n=1; n < no_molecules; n++) {
    molecule[n].start = molecule[n-1].end;
    molecule[n].end = molecule[n].start + no_atoms[n];
  }

  /* Allocate space for atoms Atom structures */

  atoms = (struct Atom *) calloc(total_no_atoms,sizeof(struct Atom));
  if( atoms == NULL ) {
    printf("Could not allocate memory for AtomList\n");
    exit(32);
  }

  /* Third pass through file -- Read+Parse Car File */
  center[0] = center[1] = center[2] = 0.0;
  rewind(CarF);
  for(n=0; n < skip; n++)
    fgets(line,MAX_LINE_LENGTH,CarF);

  for(m=0; m < no_molecules; m++) {
    for(k=molecule[m].start; k <
          molecule[m].end; k++) {

      atoms[k].molecule = m;
      atoms[k].no = k;

      fscanf(CarF,"%s %lf %lf %lf %*s %d %s %s %lf",
             atoms[k].name,
             &(atoms[k].x[0]),
             &(atoms[k].x[1]),
             &(atoms[k].x[2]),
             &(atoms[k].molecule),
             atoms[k].potential,
             atoms[k].element,
             &(atoms[k].q));

      atoms[k].x[0] += shift[0];
      atoms[k].x[1] += shift[1];
      atoms[k].x[2] += shift[2];

      if (centerflag) {
        center[0] += atoms[k].x[0];
        center[1] += atoms[k].x[1];
        center[2] += atoms[k].x[2];
      }
    }
    fgets(line,MAX_LINE_LENGTH,CarF);
    fgets(line,MAX_LINE_LENGTH,CarF);

  } /* End m (molecule) loop */

  center[0] /= (double) total_no_atoms;
  center[1] /= (double) total_no_atoms;
  center[2] /= (double) total_no_atoms;

  for (total_q=0.0,k=0; k < total_no_atoms; k++)
    total_q += atoms[k].q;

  if (pflag > 1) {
    printf("   There are %d atoms in %d molecules in this file\n",
           total_no_atoms,no_molecules);
    printf("   The total charge in the system is %7.3f.\n",total_q);
  }

  /* Search coordinates to find lowest and highest for x, y, and z */

  if (periodic == 0) {
    /* Added if/else statment STLM Oct 5 2010 */
    if (TriclinicFlag == 0) {
      /* no need to re-center the box, if we use min/max values */
      center[0] = center[1] = center[2] = 0.0;
      for ( k = 0; k < 3; k++) {
        lowest  = atoms[0].x[k];
        highest = atoms[0].x[k];

        for ( m = 1; m < total_no_atoms; m++) {
          if (atoms[m].x[k] < lowest)  lowest = atoms[m].x[k];
          if (atoms[m].x[k] > highest) highest = atoms[m].x[k];
        }
        box[0][k] = lowest  - 0.5;
        box[1][k] = highest + 0.5;
        box[2][k] = 0.0;
      }
    } else {
      printf("This tool only works for periodic systems with triclinic boxes");
      exit(32);
    }

  } else {

      if (pflag > 2)
          printf(" pbc[0] %f pbc[1] %f pbc[2] %f\n", pbc[0] ,pbc[1] ,pbc[2]);
    if (TriclinicFlag == 0) {
      for (k=0; k < 3; k++) {
        box[0][k] = -0.5*pbc[k] + center[k] + shift[k];
        box[1][k] =  0.5*pbc[k] + center[k] + shift[k];
        box[2][k] =  0.0;
      }
    } else {
      sq_c = pbc[2]*pbc[2];
      cos_alpha = cos(pbc[3]*PI_180);
      cos_gamma = cos(pbc[5]*PI_180);
      sin_gamma = sin(pbc[5]*PI_180);
      cos_beta =  cos(pbc[4]*PI_180);
      sin_beta =  sin(pbc[4]*PI_180);
      if (pflag > 2) {
        printf(" pbc[3] %f pbc[4] %f pbc[5] %f\n", pbc[3] ,pbc[4] ,pbc[5]);
        printf(" cos_alpha %f cos_beta %f cos_gamma %f\n", cos_alpha ,cos_beta ,cos_gamma);
      }
      A = pbc[0];
      B = pbc[1];
      C = pbc[2];


      box[0][0] = -0.5*A + center[0] + shift[0];
      box[1][0] =  0.5*A + center[0] + shift[0];
      box[0][1] = -0.5*B*sin_gamma + center[1] + shift[1];
      box[1][1] =  0.5*B*sin_gamma + center[1] + shift[1];
      box[0][2] = -0.5*sqrt(sq_c * sin_beta*sin_beta - C*(cos_alpha-cos_gamma*cos_beta)/sin_gamma) + center[2] + shift[2];
      box[1][2] =  0.5*sqrt(sq_c * sin_beta*sin_beta - C*(cos_alpha-cos_gamma*cos_beta)/sin_gamma) + center[2] + shift[2];
      box[2][0] =  B * cos_gamma; /* This is xy SLTM */
      box[2][1] =  C * cos_beta;  /* This is xz SLTM */
      box[2][2] =  C*(cos_alpha-cos_gamma*cos_beta)/sin_gamma; /* This is yz SLTM */
    }
  }

  /* compute image flags */

  set_box(box,hmat,hinv);

  n = 0;
  for (m = 0; m < total_no_atoms; m++) {
    double tmp;
    int w=0;

    x2lamda(atoms[m].x,lamda,hinv,box[0]);
    for (k = 0; k < 3; ++k) {
      tmp = floor(lamda[k]);
      atoms[m].image[k] = tmp;
      lamda[k] -= tmp;
      if (tmp != 0.0) ++w;
    }
    lamda2x(lamda, atoms[m].x,hmat,box[0]);
    if (w > 0) ++n;
  }

  /* warn if atoms are outside the box */
  if (n > 0) {
    if (periodic) {
      if (pflag > 1)
        printf("   %d of %d atoms with nonzero image flags\n\n",n,total_no_atoms);
    } else {
      if (iflag == 0 || (pflag > 1))
        printf("   %d of %d atoms outside the box\n\n",n,total_no_atoms);

      condexit(32);
    }
  }

  /* Close .car file */

  if (fclose(CarF) !=0) {
    printf("Error closing %s.car\n", rootname);
    exit(31);
  }
}
/* End ReadCarFile() */
