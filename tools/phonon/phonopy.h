#ifndef FFTW3
#define PHONOPY_H
#endif
#ifndef PHONOPY_H
#define PHONOPY_H

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "memory.h"
#include "qnodes.h"
#include "dynmat.h"
#include "global.h"

class Phonopy {
public:
   Phonopy(DynMat *);
   ~Phonopy();

private:
   Memory *memory;
   char str[MAXLINE];
   int npt, fftdim2;       // local variables
   int nx, ny, nz, nucell; // local variables
   int sysdim, fftdim;     // local variables
   double *mass;

   doublecomplex **FC_all;

   DynMat *dm;
   void write(int);
   void get_my_FC();
   void phonopy();
   int count_words(const char *line);
};
#endif
