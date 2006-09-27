#include "lmp2.h"

void unwrap_molecules(struct NewAtomCoordinates *coord,struct Sys *sysinfo)
{
  int iflag[5],nflag[5],icase;
  register int i,iax,imol;
  int min_true,max_true;

  if (trueflag) {

    /* Use trueflags to "unwrap" molecules */

    for (imol=0; imol < sysinfo->no_molecules; imol++) {
      for (iax=0; iax < 3; iax++) {

	min_true =  1000;
	max_true = -1000;
	for (i=sysinfo->molinfo[imol].start; i < sysinfo->molinfo[imol].end; 
	     i++) {
	  if (coord[i].truef[iax] < min_true) min_true = coord[i].truef[iax];
	  if (coord[i].truef[iax] > max_true) max_true = coord[i].truef[iax];
	}

	if ((min_true > 0) || (max_true < 0)) {
	  for (i=sysinfo->molinfo[imol].start; i < sysinfo->molinfo[imol].end; 
	       i++) coord[i].truef[iax] -= min_true;
	}

      } /* end loop over iax */
    } /* end loop over imol */
  }
  else {

    /* Use coordinates to "unwrap" molecules */

    for (imol=0; imol < sysinfo->no_molecules; imol++) {

      for (iax=0; iax < 3; iax++) {

	for (i=0; i<5; i++) { iflag[i] = 0; nflag[i] = 0;}

	for (i=sysinfo->molinfo[imol].start; i < sysinfo->molinfo[imol].end; 
	     i++) {
	  if (coord[i].fract[iax] <= 0.20) 
	    nflag[0]++;
	  else if ((coord[i].fract[iax] > 0.20) && (coord[i].fract[iax] <= 0.40))
	    nflag[1]++;
	  else if ((coord[i].fract[iax] > 0.40) && (coord[i].fract[iax] <= 0.60))
	    nflag[2]++;
	  else if ((coord[i].fract[iax] > 0.60) && (coord[i].fract[iax] <= 0.80))
	    nflag[3]++;
	  else if ((coord[i].fract[iax] > 0.80) && (coord[i].fract[iax] <= 1.00))
	    nflag[4]++;
	}
      
	for (i=0; i<5; i++) { if (nflag[i] > 0) iflag[i] = 1; }

	icase = 10000*iflag[0] + 1000*iflag[1] + 100*iflag[2] + 10*iflag[3] + iflag[4];

	switch (icase) {
	 case 10001:
	   if (nflag[0] > nflag[4]) {
	     for (i=sysinfo->molinfo[imol].start; i < sysinfo->molinfo[imol].end; 
		  i++) if (coord[i].fract[iax] > 0.80) coord[i].fract[iax] -= 1.0;
	   }
	   else {
	     for (i=sysinfo->molinfo[imol].start; i < sysinfo->molinfo[imol].end; 
		  i++) if (coord[i].fract[iax] <= 0.20) coord[i].fract[iax] += 1.0;
	   }
	   break;
         case 11001:
	   for (i=sysinfo->molinfo[imol].start; i < sysinfo->molinfo[imol].end; 
		i++) if (coord[i].fract[iax] > 0.80) coord[i].fract[iax] -= 1.0;
	   break;
         case 10011:
	   for (i=sysinfo->molinfo[imol].start; i < sysinfo->molinfo[imol].end; 
		i++) if (coord[i].fract[iax] < 0.20) coord[i].fract[iax] += 1.0;
	   break;
         case 11101:
	   for (i=sysinfo->molinfo[imol].start; i < sysinfo->molinfo[imol].end; 
		i++) if (coord[i].fract[iax] > 0.80) coord[i].fract[iax] -= 1.0;
	   break;
         case 10111:
	   for (i=sysinfo->molinfo[imol].start; i < sysinfo->molinfo[imol].end; 
		i++) if (coord[i].fract[iax] < 0.20) coord[i].fract[iax] += 1.0;
	   break;
         case 11011:
	   if ((nflag[0]+nflag[1]) > (nflag[3]+nflag[4])) {
	     for (i=sysinfo->molinfo[imol].start; i < sysinfo->molinfo[imol].end; 
		  i++) if (coord[i].fract[iax] > 0.60) coord[i].fract[iax] -= 1.0;
	   }
	   else {
	     for (i=sysinfo->molinfo[imol].start; i < sysinfo->molinfo[imol].end; 
		  i++) if (coord[i].fract[iax] <= 0.40) coord[i].fract[iax] += 1.0;
	   }
	   break;
         default:
	   break;
	} /* end switch */

      } /* end loop over imol*/
    } /* end loop over iax */

  } /* end if on trueflag */
}
