/****************************** ProcessDumpFile.c ***************************
*
*   Created 8/98 by John Carpenter
*   Based on routines ReadLMPDump and CreateArcFile written by
*   Mike Peachey during summer 1997
*
*   This function reads the output from LAMMPS dumpfile
*   and then writes it back out to an ACCELRYS arcfile using
*   info from ReadCarfile.
*
*/

#include "lmp2.h"

void ProcessPosFile01(int num_posfiles,char *posnames[],struct Sys *sysinfo,FILE *ArcFile)
{
  struct Boundary cell;
  struct NewAtomCoordinates *coord;

  FILE *PosFile;
  
  char line[MAX_LINE_LENGTH];

  int i,iax,itag,itype;
  int nxx,nyy,nzz,matoms;
  int frame, mskip,timestep;
  int ifile;

  float xx,yy,zz;

  /* Function prototype declarations */

  extern void unwrap_molecules(struct NewAtomCoordinates *,struct Sys *);
  extern void WriteArcFrame(FILE *,int,int,struct Sys *);

  /* Begin execution */

  coord = calloc(sysinfo->natoms, sizeof(struct NewAtomCoordinates));
  if ( coord == NULL )
  {
      fprintf(stderr,"Memory Allocation Problem(coord), exiting program");
      exit(3);
  }

  /* process data */

  fprintf(stderr,"\n Processing Timesteps\n\n");

  frame = 0;
  mskip = 0;

  for (ifile=0; ifile < num_posfiles; ifile++) {

  fprintf(stderr,"\n Opening file \n");

    if ( (PosFile = fopen(posnames[ifile],"r")) == NULL ) {
      fprintf(stderr,"Cannot open %s\n",posnames[ifile]);
      exit(2);
    }


  fprintf(stderr,"\n File open \n");

    while (fgets(line,MAX_LINE_LENGTH, PosFile) != NULL) {
      if (strstr(line,"NUMBER OF ATOMS")) {   
	fgets(line,MAX_LINE_LENGTH, PosFile);
	matoms = atoi(line);
	if (matoms != sysinfo->natoms) {
	  fprintf(stderr,"Number of atoms in car and dump files do not match/n");
	  exit(1);
	}
      }
      else if (strstr(line,"BOX BOUNDS")) {

	for (iax=0; iax < 3; iax++) {
	  fgets(line, MAX_LINE_LENGTH, PosFile);
	  cell.low[iax] = atof(strtok(line, " "));
	  cell.hi[iax]  = atof(strtok(NULL, " \0\n"));
	  cell.size[iax] = cell.hi[iax] - cell.low[iax];
	  sysinfo->celldim[iax] = cell.size[iax];
	}
      }
      else if (strstr(line,"TIMESTEP")) {

	fscanf(PosFile,"%d",&timestep);
      }
      else if (strstr(line,"ATOMS")) {
       
	if (trueflag) {
	  for (i=0; i<matoms; i++) {
	    fgets(line, MAX_LINE_LENGTH, PosFile);
	    sscanf(line,"%d %d %f %f %f %d %d %d",&itag,&itype,
		   &xx,&yy,&zz,&nxx,&nyy,&nzz);

	    coord[itag-1].type = itype;
	    coord[itag-1].fract[0]   = xx;
	    coord[itag-1].fract[1]   = yy;
	    coord[itag-1].fract[2]   = zz;
	    coord[itag-1].truef[0]   = nxx;
	    coord[itag-1].truef[1]   = nyy;
	    coord[itag-1].truef[2]   = nzz;
	  }

	  if (move_molecules) unwrap_molecules(coord,sysinfo);

	  for (i=0; i<matoms; i++) {
	    for (iax=0; iax<3; iax++)
	      coord[i].fract[iax] += coord[i].truef[iax];
	  }
	}
	else {
	 
	  for (i=0; i<matoms; i++) {
	    fgets(line, MAX_LINE_LENGTH, PosFile);
	    sscanf(line,"%d %d %f %f %f",&itag,&itype,&xx,&yy,&zz);

	    coord[itag-1].type = itype;
	    coord[itag-1].fract[0]   = xx;
	    coord[itag-1].fract[1]   = yy;
	    coord[itag-1].fract[2]   = zz;
	  }

	  if (move_molecules) unwrap_molecules(coord,sysinfo);

	} /* end if on trueflag */

	for (i=0; i<matoms; i++) {
	  sysinfo->atoms[i].xyz[0] = cell.size[0]*coord[i].fract[0];
	  sysinfo->atoms[i].xyz[1] = cell.size[1]*coord[i].fract[1];
	  sysinfo->atoms[i].xyz[2] = cell.size[2]*coord[i].fract[2];
	}

	/* Ready to write arc file entry */

	if ((mskip >= nskip) || (frame == 0)) {

	  WriteArcFrame(ArcFile,frame,timestep,sysinfo);
	  frame++;
	  mskip = 0;

	  if (frame%20 == 0) fprintf(stderr," %d",frame);
	  if (frame%400 == 0) fprintf(stderr,"\n");

	}
	else
	  mskip++;

      } /* end if on ITEM */

    } /* end while over LINES*/
    
    fclose(PosFile);
   
  } /* end for loop over POS FILES */
  
  fprintf(stderr,"\n\n %d frames were written to the ArcFile\n",frame);
}

	           
