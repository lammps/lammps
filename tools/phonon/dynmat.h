#ifndef DYNMAT_H
#define DYNMAT_H

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "memory.h"
#include "interpolate.h"

extern "C"{
#include "f2c.h"
#include "clapack.h"
}

using namespace std;

class DynMat {
public:

  DynMat(int, char**);
  ~DynMat();

  int nx, ny, nz, nucell;
  int sysdim, fftdim;
  double eml2f;
  char *funit;

  void getDMq(double *);
  void getDMq(double *, double *);
  void writeDMq(double *);
  void writeDMq(double *, const double, FILE *fp);
  int geteigen(double *, int);
  void reset_interp_method();

  doublecomplex **DM_q;

  int flag_latinfo;
  double Tmeasure, basevec[9], ibasevec[9];
  double **basis;
  int *attyp;

private:

  int flag_skip, flag_reset_gamma;
  Interpolate *interpolate;
  
  Memory *memory;
  int npt, fftdim2;

  int nasr;
  void EnforceASR();

  char *binfile, *dmfile;
  double boltz, q[3];
  double *M_inv_sqrt;

  doublecomplex **DM_all;

  void car2dir(int); // to convert basis from cartisian coordinate into factional.
  void real2rec();
  void GaussJordan(int, double *);

  void help();
  void ShowVersion();
};
#endif
