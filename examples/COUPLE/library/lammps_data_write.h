#ifndef LAMMPS_DATA_WRITE_H
#define LAMMPS_DATA_WRITE_H

#include "send2one.h"
#include "stdio.h"

class LAMMPSDataWrite : public Send2One {
 public:
  LAMMPSDataWrite(MPI_Comm);
  ~LAMMPSDataWrite();

  void pre();
  int size();
  void pack(char *);
  void process(int, char *);
  void post();

  void file(char *);
  void header(char *, int);
  void header(char *, double);
  void header(char *, double, double);
  void atoms(int);
  void atoms(int *);
  void atoms(double *);
  void atoms(int, double **);

 private:
  char *outfile;
  int nlocal;
  FILE *fp;

  int nheader,maxheader;
  char **format;
  int *headtype,*ihead;
  double *dhead;
  double **ddhead;

  int nper,maxper;
  int *atomtype;
  int **ivec;
  double **dvec;
  int *stride;

  void grow_header();
  void grow_peratom();
};

#endif
