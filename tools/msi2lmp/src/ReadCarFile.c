/*
*  This function opens the .car file and extracts coordinate information
*  into the atoms Atom structure
*/

#include "Msi2LMP2.h"

void ReadCarFile(void)
{
   char line[MAX_LINE_LENGTH];  /* Stores lines as they are read in */
   int k,m,n;			/* counters */	
   int skip;			/* lines to skip at beginning of file */
   double  lowest, highest;	/* temp coordinate finding variables */
   double total_q;

/* Open .car file for reading */

   sprintf(line,"%s.car",rootname);
   if (pflag > 0) fprintf(stderr," Reading car file: %s\n",line);
   if( (CarF = fopen(line,"r")) == NULL ) {
      fprintf(stderr,"Cannot open %s\n",line);
      exit(2);
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
      if(pbc[3] != 90.0 || pbc[4] != 90.0 || pbc[5] != 90.0) {
         fprintf(stderr,"The system is not rectangular- LAMMPS can't handle it!!");
         exit(2);
      }
   }
   else {
      periodic = 0;
      skip = 4;
      if (pflag > 1) {
	fprintf(stderr,"   %s is not a periodic system\n", rootname);
	fprintf(stderr,"   Assigning cell parameters based on coordinates\n");
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
     fprintf(stderr,"Could not allocate memory for no_atoms\n");
     exit(2);
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
     fprintf(stderr,"Unable to allocate memory for molecule structure\n");
     exit(2);
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
     fprintf(stderr,"Could not allocate memory for AtomList\n");
     exit(2);
   }
   
/* Third pass through file -- Read+Parse Car File */

   rewind(CarF);
   for(n=0; n < skip; n++)
     fgets(line,MAX_LINE_LENGTH,CarF);

   for(m=0; m < no_molecules; m++) {
     for(k=molecule[m].start; k <
	   molecule[m].end; k++) {

       atoms[k].molecule = m;
       atoms[k].no = k;

       fscanf(CarF,"%s %lf %lf %lf %*s %*s %s %s %f",
	      atoms[k].name,
	      &(atoms[k].x[0]),
	      &(atoms[k].x[1]),
	      &(atoms[k].x[2]),
	      atoms[k].potential,
	      atoms[k].element,
	      &(atoms[k].q));
     }
     fgets(line,MAX_LINE_LENGTH,CarF);
     fgets(line,MAX_LINE_LENGTH,CarF);

   } /* End m (molecule) loop */

   for (total_q=0.0,k=0; k < total_no_atoms; k++)
     total_q += atoms[k].q;

   if (pflag > 1) {
     fprintf(stderr,"   There are %d atoms in %d molecules in this file\n",
	     total_no_atoms,no_molecules);
     fprintf(stderr,"   The total charge in the system is %7.3f.\n\n",total_q);
   }

/* Search coordinates to find lowest and highest for x, y, and z */

   if (periodic == 0) {
     for ( k = 0; k < 3; k++) {
       lowest  = atoms[0].x[k];
       highest = atoms[0].x[k];

       for ( m = 1; m < total_no_atoms; m++) {
         if (atoms[m].x[k] < lowest)  lowest = atoms[m].x[k];
         if (atoms[m].x[k] > highest) highest = atoms[m].x[k];
       }
       pbc[k] = lowest;
       pbc[k+3] = highest;
     }
   }
   else {
     for (k=0; k < 3; k++) {
       pbc[k+3] = pbc[k];
       pbc[k] = 0.0;
     }
   }

/* Close .car file */
   
   if (fclose(CarF) !=0) {
      fprintf(stderr,"Error closing %s.car\n", rootname);
      exit(1);
   }
}
/* End ReadCarFile() */
