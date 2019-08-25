/****************************** ReadCarFile.c ******************************
*
*  This function opens the .car file and extracts coordinate information
*  into the atoms Atom structure
*
*  THIS FUNCTION HAS BEEN MODIFIED FROM Msi2LMP2 AND USES THE ATOM
*  STRUCTURE MEMBERS DIFFERENTLY
*/

#include "lmp2.h"

void ReadCarFile(FILE *CarF,struct Sys *sysinfo)
{
   char line[MAX_LINE_LENGTH];  /* Stores lines as they are read in */
   int  k,m,n;			/* counters */	
   int skip;			/* lines to skip at beginning of file */


/* Determine Number of molecules & atoms */

   rewind(CarF);
   sysinfo->no_molecules = -1; /* Set to -1 because counter will be 
			 incremented an extra time at the end of the file */

   fgets(line,MAX_LINE_LENGTH,CarF); /* Read header line */

/* Check for periodicity, set periodic and skip */

   if( strncmp(fgets(line,MAX_LINE_LENGTH,CarF),"PBC=ON",6) == 0)
   {
      sysinfo->periodic = 1;
      skip = 5; /* Data starts 5 lines from beginning of file */
   }
   else
   {
      sysinfo->periodic = 0;
      skip = 4;
   }

/* First pass through file -- Count molecules */

   while( fgets(line,MAX_LINE_LENGTH,CarF) != NULL )
      if( strncmp(line,"end",3) == 0 ) 
         sysinfo->no_molecules++;

/* Allocate space to keep track of the number of atoms within a molecule */

   sysinfo->molinfo = 
     (struct Mol *) calloc(sysinfo->no_molecules,sizeof(struct Mol));

   if ( sysinfo->molinfo == NULL )
   {
      puts("Could not allocate memory for molinfo");
      exit(2);
   }

/* Second pass through file -- Count atoms */

   rewind(CarF);
   for(n=0; n < skip; n++)               /* Skip beginning lines */
      fgets(line,MAX_LINE_LENGTH,CarF);

   sysinfo->molinfo[0].start = 0;
   sysinfo->molinfo[0].end   = 0;
   for(n=0; n < sysinfo->no_molecules; n++)
   {
      if (n != 0) {
         sysinfo->molinfo[n].start = sysinfo->molinfo[n-1].end;
	 sysinfo->molinfo[n].end   = sysinfo->molinfo[n-1].end;
      }
      while( strncmp(fgets(line,MAX_LINE_LENGTH,CarF),"end",3) ) 
         sysinfo->molinfo[n].end++;
   }
   n = sysinfo->no_molecules-1;
   sysinfo->natoms = sysinfo->molinfo[n].end; 

/* Allocate space for atoms Atom structures */

   sysinfo->atoms = (struct Atom *) calloc(sysinfo->natoms,
                                    sizeof(struct Atom));
   if( sysinfo->atoms == NULL )
   {
      puts("Could not allocate memory for AtomList");
      exit(2);
   }

/* Third pass through file -- Read+Parse Car File */

   rewind(CarF);

   for(n=0; n < skip; n++)
      fgets(line,MAX_LINE_LENGTH,CarF);

   for(m=0; m < sysinfo->no_molecules; m++)
   { /* m loops through molecules */
      for(n=sysinfo->molinfo[m].start; n < sysinfo->molinfo[m].end; n++)
      { /* n loops through atoms within a molecule */

         sysinfo->atoms[n].molecule = m;

         fscanf(CarF,"%s %lf %lf %lf %s %s %s %s %f",
                sysinfo->atoms[n].name,
                &(sysinfo->atoms[n].xyz[0]),
                &(sysinfo->atoms[n].xyz[1]),
                &(sysinfo->atoms[n].xyz[2]),
		sysinfo->atoms[n].res_name,
                sysinfo->atoms[n].res_num,
                sysinfo->atoms[n].potential,
                sysinfo->atoms[n].element,
                &(sysinfo->atoms[n].q));
      } /* End n (atoms) loop */

   fgets(line,MAX_LINE_LENGTH,CarF);
   fgets(line,MAX_LINE_LENGTH,CarF);

   } /* End m (molecule) loop */

} /* End ReadCarFile() */

