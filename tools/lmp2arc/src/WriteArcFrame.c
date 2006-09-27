#include <time.h>
#include "lmp2.h"

void WriteArcFrame(FILE *ArcFile,int frame,int timestep,
                     struct Sys *sysinfo)
{
  time_t tim;
  char *datestr;
  const char MSIformat[] = "!BIOSYM archive 3";
  const char PBC_ON[] = "PBC=ON";
  const char PBC_OFF[] = "PBC=OFF";

  const float perp = 90.0;
  register int i,m;
  float number_of_psecs;

  /* Function prototypes */

  void dot(int);

  /* Begin execution */

  number_of_psecs = (float) timestep/(float) npico;

  tim = time(NULL);
  datestr = asctime(localtime(&tim));

  if (frame == 0) 
  {
    fprintf(ArcFile,"%s\n",MSIformat);
    if (sysinfo->periodic)
      fprintf(ArcFile,"%s\n",PBC_ON);
    else
      fprintf(ArcFile,"%s\n",PBC_OFF);
  }

  fprintf(ArcFile,"Frame %d  Number of Picoseconds = %10.4f\n",frame,
                 number_of_psecs);

  fprintf(ArcFile,"!DATE %s",datestr);

  if (sysinfo->periodic) 
   fprintf(ArcFile,"PBC%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f  (P1)\n", 
			      sysinfo->celldim[0],sysinfo->celldim[1],
			      sysinfo->celldim[2],perp,perp,perp);

  for (m=0; m<sysinfo->no_molecules; m++) 
  {
    for (i=sysinfo->molinfo[m].start; i<sysinfo->molinfo[m].end; i++)
    {
      fprintf(ArcFile,"%-5s %14.9f %14.9f %14.9f %-4s %-7s%-7s %-2s% 6.4f\n",
              sysinfo->atoms[i].name,
	      sysinfo->atoms[i].xyz[0],sysinfo->atoms[i].xyz[1],
	      sysinfo->atoms[i].xyz[2],sysinfo->atoms[i].res_name,
	      sysinfo->atoms[i].res_num,sysinfo->atoms[i].potential,
	      sysinfo->atoms[i].element,sysinfo->atoms[i].q);
    }
    fprintf(ArcFile, "end\n");
   }
   fprintf( ArcFile, "end\n");

   fflush(ArcFile);
}
