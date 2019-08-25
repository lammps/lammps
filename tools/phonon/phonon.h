#ifndef PHONON_H
#define PHONON_H

#include "stdio.h"
#include "stdlib.h"
#include <complex>
#include "dynmat.h"
#include "memory.h"

using namespace std;

class Phonon{
public:
  Phonon(DynMat *);
  ~Phonon();

  DynMat *dynmat;

private:
  int nq, ndim, sysdim;
  double **qpts, *wt;
  double **eigs;

  int ndos, nlocal, *locals;
  double *dos, fmin, fmax, df, rdf;
  double ***ldos;

  Memory *memory;

  void QMesh();
  void ComputeAll();

  void pdos();
  void pdisp();
  void therm();

  void ldos_egv();
  void ldos_rsgf();
  void local_therm();

  void dmanyq();
  void vfanyq();
  void DMdisp();
  void vecanyq();

  void ShowCell();

  void smooth(double *, const int);
  void writeDOS();
  void writeLDOS();
  void Normalize();

  int count_words(const char *);

#ifdef UseSPG
  int num_atom, *attyp;
  double latvec[3][3], **atpos;
#endif
};

#endif
