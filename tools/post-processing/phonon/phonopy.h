#ifndef FFTW3
#define PHONOPY_H
#endif
#ifndef PHONOPY_H
#define PHONOPY_H

#include "zheevd.h"

class Phonopy {
public:
   Phonopy(class DynMat *);
   ~Phonopy();

private:
   class UserInput *input;
   class Memory *memory;
   int npt, fftdim2;       // local variables
   int nx, ny, nz, nucell; // local variables
   int sysdim, fftdim;     // local variables
   double *mass;

   doublecomplex **FC_all;

   class DynMat *dm;
   void write(int);
   void get_my_FC();
   void phonopy();
   int count_words(const char *line);
};
#endif
