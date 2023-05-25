#ifndef DYNMAT_H
#define DYNMAT_H

#include "zheevd.h"

#include <cstdio>

class DynMat {
public:

  DynMat(int, char**);
  ~DynMat();

  int nx, ny, nz, nucell;
  int sysdim, fftdim;
  double eml2f, eml2fc, symprec;
  char *funit;

  void getDMq(double *);
  void getDMq(double *, double *);
  void writeDMq(double *);
  void writeDMq(double *, const double, FILE *fp);
  int geteigen(double *, int);
  void reset_interp_method();

  doublecomplex **DM_q;

  int flag_latinfo;
  int npt, fftdim2;
  double Tmeasure, basevec[9], ibasevec[9];
  double *M_inv_sqrt;
  double **basis;
  int *attyp;

  class UserInput *input;

private:

  int flag_skip, flag_reset_gamma;
  class Interpolate *interpolate;
  class Memory *memory;

  void EnforceASR();

  char *binfile, *dmfile;
  double boltz;

  doublecomplex **DM_all;

  void car2dir();      // to convert basis from cartisian coordinate into factional.
  void real2rec();
  void GaussJordan(int, double *);

  void help();
  void ShowInfo();
  void ShowVersion();
  void Define_Conversion_Factor();
};
#endif
