#include "stdlib.h"
#include "string.h"
#include "lammps_data_write.h"
#include "memory.h"
#include "error.h"

#define DELTA 4;

enum{INT,DOUBLE,DOUBLE2};

/* ---------------------------------------------------------------------- */

LAMMPSDataWrite::LAMMPSDataWrite(MPI_Comm caller_comm) : Send2One(caller_comm)
{
  outfile = NULL;

  nheader = maxheader = 0;
  format = NULL;
  headtype = NULL;
  ihead = NULL;
  dhead = NULL;
  ddhead = NULL;

  nper = maxper = 0;
  atomtype = NULL;
  ivec = NULL;
  dvec = NULL;
  stride = NULL;
}

/* ---------------------------------------------------------------------- */

LAMMPSDataWrite::~LAMMPSDataWrite()
{
  delete [] outfile;

  for (int i = 0; i < nheader; i++) delete [] format[i];
  memory->sfree(format);
  memory->sfree(headtype);
  memory->sfree(ihead);
  memory->sfree(dhead);
  memory->destroy_2d_double_array(ddhead);

  memory->sfree(atomtype);
  memory->sfree(ivec);
  memory->sfree(dvec);
  memory->sfree(stride);
}

/* ---------------------------------------------------------------------- */

void LAMMPSDataWrite::pre()
{
  if (me == 0) {
    fp = fopen(outfile,"w");
    if (fp == NULL)
      error->one("Could not open data_write output file");
  }

  if (me == 0) {
    fprintf(fp,"%s","LAMMPS data file\n\n");
    for (int i = 0; i < nheader; i++)
      if (headtype[i] == INT) fprintf(fp,format[i],ihead[i]);
      else if (headtype[i] == DOUBLE) fprintf(fp,format[i],dhead[i]);
      else if (headtype[i] == DOUBLE2) fprintf(fp,format[i],
					       ddhead[i][0],ddhead[i][1]);
    fprintf(fp,"\nAtoms\n\n");
  }
}

/* ---------------------------------------------------------------------- */

int LAMMPSDataWrite::size()
{
  return nper*nlocal*sizeof(double);
}

/* ---------------------------------------------------------------------- */

void LAMMPSDataWrite::pack(char *cbuf)
{
  int i,j;

  double *dbuf = (double *) cbuf;

  int m = 0;
  for (i = 0; i < nlocal; i++)
    for (j = 0; j < nper; j++) {
      if (i == 0) {
	if (atomtype[j] == 0) 
	  printf("AT %d %d %p %d\n",
		 atomtype[j],stride[j],ivec[j],ivec[j][0]);
	else
	  printf("AT %d %d %p %g\n",
		 atomtype[j],stride[j],dvec[j],dvec[j][0]);
      }
      if (atomtype[j] == INT) dbuf[m++] = ivec[j][i*stride[j]];
      else if (atomtype[j] == DOUBLE) dbuf[m++] = dvec[j][i*stride[j]];
    }
}

/* ---------------------------------------------------------------------- */

void LAMMPSDataWrite::process(int nbuf, char *cbuf)
{
  int i,j;

  double *dbuf = (double *) cbuf;
  int n = nbuf/nper/sizeof(double);

  int m = 0;
  for (i = 0; i < n; i++) {
    for (j = 0; j < nper; j++) {
      double dvalue = dbuf[m++];
      if (atomtype[j] == INT) fprintf(fp,"%d ",static_cast<int> (dvalue));
      else if (atomtype[j] == DOUBLE) fprintf(fp,"%g ",dvalue);
    }
    fprintf(fp,"\n");
  }
}

/* ---------------------------------------------------------------------- */

void LAMMPSDataWrite::post()
{
  if (me == 0) fclose(fp);
}

/* ---------------------------------------------------------------------- */

void LAMMPSDataWrite::file(char *fname)
{
  int n = strlen(fname) + 1;
  outfile = new char[n];
  strcpy(outfile,fname);
}

/* ---------------------------------------------------------------------- */

void LAMMPSDataWrite::header(char *str, int ivalue)
{
  if (nheader == maxheader) grow_header();
  int n = strlen(str) + 2;
  format[nheader] = new char[n];
  strcpy(format[nheader],str);
  format[nheader][n-2] = '\n';
  format[nheader][n-1] = '\0';
  headtype[nheader] = INT;
  ihead[nheader] = ivalue;
  nheader++;
}

/* ---------------------------------------------------------------------- */

void LAMMPSDataWrite::header(char *str, double dvalue)
{
  if (nheader == maxheader) grow_header();
  int n = strlen(str) + 2;
  format[nheader] = new char[n];
  strcpy(format[nheader],str);
  format[nheader][n-2] = '\n';
  format[nheader][n-1] = '\0';
  headtype[nheader] = DOUBLE;
  dhead[nheader] = dvalue;
  nheader++;
}

/* ---------------------------------------------------------------------- */

void LAMMPSDataWrite::header(char *str, double dvalue1, double dvalue2)
{
  if (nheader == maxheader) grow_header();
  int n = strlen(str) + 2;
  format[nheader] = new char[n];
  strcpy(format[nheader],str);
  format[nheader][n-2] = '\n';
  format[nheader][n-1] = '\0';
  headtype[nheader] = DOUBLE2;
  ddhead[nheader][0] = dvalue1;
  ddhead[nheader][1] = dvalue2;
  nheader++;
}

/* ---------------------------------------------------------------------- */

void LAMMPSDataWrite::atoms(int n)
{
  nlocal = n;
}

/* ---------------------------------------------------------------------- */

void LAMMPSDataWrite::atoms(int *vec)
{
  if (nper == maxper) grow_peratom();
  atomtype[nper] = INT;
  ivec[nper] = vec;
  stride[nper] = 1;
  nper++;
}

/* ---------------------------------------------------------------------- */

void LAMMPSDataWrite::atoms(double *vec)
{
  if (nper == maxper) grow_peratom();
  atomtype[nper] = DOUBLE;
  dvec[nper] = vec;
  stride[nper] = 1;
  nper++;
}

/* ---------------------------------------------------------------------- */

void LAMMPSDataWrite::atoms(int n, double **vec)
{
  if (nper+n >= maxper) grow_peratom();
  for (int i = 0; i < n; i++) {
    atomtype[nper] = DOUBLE;
    dvec[nper] = &vec[0][i];
    stride[nper] = n;
    nper++;
  }
}

/* ---------------------------------------------------------------------- */

void LAMMPSDataWrite::grow_header()
{
  int n = maxheader + DELTA;
  format = (char **) memory->srealloc(format,n*sizeof(char *),"ldw:format");
  headtype = (int *) memory->srealloc(headtype,n*sizeof(int),"ldw:headtype");
  ihead = (int *) memory->srealloc(ihead,n*sizeof(int),"ldw:ihead");
  dhead = (double *) memory->srealloc(dhead,n*sizeof(double),"ldw:dhead");
  ddhead = memory->grow_2d_double_array(ddhead,n,2,"ldw:ddhead");
  maxheader = n;
}

/* ---------------------------------------------------------------------- */

void LAMMPSDataWrite::grow_peratom()
{
  int n = maxper + DELTA;
  atomtype = (int *) memory->srealloc(atomtype,n*sizeof(int *),"ldw:atomtype");
  ivec = (int **) memory->srealloc(ivec,n*sizeof(int *),"ldw:ihead");
  dvec = (double **) memory->srealloc(dvec,n*sizeof(double *),"ldw:dhead");
  stride = (int *) memory->srealloc(stride,n*sizeof(int *),"ldw:stride");
  maxper = n;
}
