/*
*   This routine reads the data from a .frc forcefield file and stores it in
*   dynamically allocated memory. This allows for fast searches of the
*   file.
*
*/

#define FF_MAIN

#include "Forcefield.h"
#include "Msi2LMP2.h"

void ReadFrcFile(void)
{
   extern void InitializeItems(void);
   extern void SearchAndFill(struct FrcFieldItem *item);

  /* Open Forcefield File */
   if ( (FrcF = fopen(FrcFileName,"r")) == NULL ) {
      fprintf(stderr,"Cannot open %s\n", FrcFileName);
      exit(2);
   }
   InitializeItems(); /* sets keywords, number of members and number of
			 parameters for each structure */
  /* allocate memory to and search and fill each structure */


   SearchAndFill(&ff_atomtypes);
   SearchAndFill(&equivalence);
   SearchAndFill(&ff_vdw);
   SearchAndFill(&ff_bond);
   SearchAndFill(&ff_ang);
   SearchAndFill(&ff_tor);
   SearchAndFill(&ff_oop);
   
   if (forcefield != 1) {  /* Skip cross terms for class I */
     SearchAndFill(&ff_bonbon);
     SearchAndFill(&ff_bonang);
     SearchAndFill(&ff_angtor);
     SearchAndFill(&ff_angangtor);
     SearchAndFill(&ff_endbontor);
     SearchAndFill(&ff_midbontor);
     SearchAndFill(&ff_bonbon13);
     SearchAndFill(&ff_angang);
   }
   if (pflag > 1) { 

     fprintf(stderr,"\n Item %s has %d entries\n",
	     ff_atomtypes.keyword,ff_atomtypes.entries);
     fprintf(stderr," Item %s has %d entries\n",
	     equivalence.keyword,equivalence.entries);
     fprintf(stderr," Item %s has %d entries\n",
	     ff_vdw.keyword,ff_vdw.entries);
     fprintf(stderr," Item %s has %d entries\n",
	     ff_bond.keyword,ff_bond.entries);
     fprintf(stderr," Item %s has %d entries\n",
	     ff_ang.keyword,ff_ang.entries);
     if (forcefield > 1) {
       fprintf(stderr," Item %s has %d entries\n",
	       ff_bonbon.keyword,ff_bonbon.entries);
       fprintf(stderr," Item %s has %d entries\n",
	       ff_bonang.keyword,ff_bonang.entries);
     }
     fprintf(stderr," Item %s has %d entries\n",
	     ff_tor.keyword,ff_tor.entries);
     if (forcefield > 1) {
       fprintf(stderr," Item %s has %d entries\n",
	       ff_angtor.keyword,ff_angtor.entries);
       fprintf(stderr," Item %s has %d entries\n",
	       ff_angangtor.keyword,ff_angangtor.entries);
       fprintf(stderr," Item %s has %d entries\n",
	       ff_endbontor.keyword,ff_endbontor.entries);
       fprintf(stderr," Item %s has %d entries\n",
	       ff_midbontor.keyword,ff_midbontor.entries);
       fprintf(stderr," Item %s has %d entries\n",
	       ff_bonbon13.keyword,ff_bonbon13.entries);
     }
     fprintf(stderr," Item %s has %d entries\n",
	     ff_oop.keyword,ff_oop.entries);
     if (forcefield > 1) {
       fprintf(stderr," Item %s has %d entries\n",
	       ff_angang.keyword,ff_angang.entries);
     }
     fprintf(stderr,"\n");
   }
   fclose(FrcF);
}

