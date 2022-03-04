#ifndef UseSPG
#define KPATH_H
#endif
#ifndef KPATH_H
#define KPATH_H

#include "qnodes.h"
#include "dynmat.h"
#include "memory.h"

class kPath{
public:

  kPath(DynMat *, QNodes *);
  ~kPath();

  void kpath();
  void show_path();
  void show_info();

private:

  Memory *memory;

  DynMat *dynmat;
  QNodes *q;
  char symbol[11];
  int spgnum, sysdim, fftdim, num_atom, *attyp;
  double latvec[3][3], **atpos;

};
#endif
