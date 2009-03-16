/* -----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------ */

// Convert one or more LAMMPS binary dump files to ASCII text files
//
// this serial code must be compiled on a platform that can read the binary
//   dump files since binary formats are not compatible across all platforms
//
// Syntax: binary2txt file1 file2 ...
// Creates:           file1.txt file2.txt ...

#include "stdio.h"
#include "stdlib.h"
#include "string.h"

main(int narg, char **arg)
{
  int i,j,k,m,n;
  int ntimestep,natoms,size_one,nchunk,triclinic;
  double xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz;
  int maxbuf = 0;
  double *buf = NULL;

  if (narg == 1) {
    printf("Syntax: binary2txt file1 file2 ...\n");
    exit(1);
  }

  // loop over files

  for (int iarg = 1; iarg < narg; iarg++) {
    printf("%s:",arg[iarg]);
    fflush(stdout);
    FILE *fp = fopen(arg[iarg],"rb");
    if (!fp) {
      printf("ERROR: Could not open %s\n",arg[iarg]);
      exit(1);
    }

    n = strlen(arg[iarg]) + 1 + 4;
    char *filetxt = new char[n];
    strcpy(filetxt,arg[iarg]);
    strcat(filetxt,".txt");
    FILE *fptxt = fopen(filetxt,"w");
    delete [] filetxt;

    // loop over snapshots in file

    while (1) {

      fread(&ntimestep,sizeof(int),1,fp);

      // detect end-of-file

      if (feof(fp)) {
	fclose(fp);
	fclose(fptxt);
	break;
      }

      fread(&natoms,sizeof(int),1,fp);
      fread(&triclinic,sizeof(int),1,fp);
      fread(&xlo,sizeof(double),1,fp);
      fread(&xhi,sizeof(double),1,fp);
      fread(&ylo,sizeof(double),1,fp);
      fread(&yhi,sizeof(double),1,fp);
      fread(&zlo,sizeof(double),1,fp);
      fread(&zhi,sizeof(double),1,fp);
      if (triclinic) {
	fread(&xy,sizeof(double),1,fp);
	fread(&xz,sizeof(double),1,fp);
	fread(&yz,sizeof(double),1,fp);
      }
      fread(&size_one,sizeof(int),1,fp);
      fread(&nchunk,sizeof(int),1,fp);
      
      fprintf(fptxt,"ITEM: TIMESTEP\n");
      fprintf(fptxt,"%d\n",ntimestep);
      fprintf(fptxt,"ITEM: NUMBER OF ATOMS\n");
      fprintf(fptxt,"%d\n",natoms);
      if (!triclinic) {
	fprintf(fptxt,"ITEM: BOX BOUNDS\n");
	fprintf(fptxt,"%g %g\n",xlo,xhi);
	fprintf(fptxt,"%g %g\n",ylo,yhi);
	fprintf(fptxt,"%g %g\n",zlo,zhi);
      } else {
	fprintf(fptxt,"ITEM: BOX BOUNDS xy xz yz\n");
	fprintf(fptxt,"%g %g %g\n",xlo,xhi,xy);
	fprintf(fptxt,"%g %g %g\n",ylo,yhi,xz);
	fprintf(fptxt,"%g %g %g\n",zlo,zhi,yz);
      }
      fprintf(fptxt,"ITEM: ATOMS\n");

      // loop over processor chunks in file

      for (i = 0; i < nchunk; i++) {
	fread(&n,sizeof(int),1,fp);

	// extend buffer to fit chunk size
	
	if (n > maxbuf) {
	  if (buf) delete [] buf;
	  buf = new double[n];
	  maxbuf = n;
	}

	// read chunk and write as size_one values per line

	fread(buf,sizeof(double),n,fp);
	n /= size_one;
	m = 0;
	for (j = 0; j < n; j++) {
	  for (k = 0; k < size_one; k++) fprintf(fptxt,"%g ",buf[m++]);
	  fprintf(fptxt,"\n");
	}
      }

      printf(" %d",ntimestep);
      fflush(stdout);
    }
    printf("\n");
  }

  if (buf) delete [] buf;
}
