#ifndef ONE2MANY_H
#define ONE2MANY_H

#include "mpi.h"

#include <map>

class One2Many {
 public:
  One2Many(MPI_Comm);
  ~One2Many();

  void setup(int, int, int *);
  void scatter(double *, int, double *);

 protected:
  int me,nprocs;
  MPI_Comm comm;
  class Memory *memory;
  std::map<int,int> *hash;
  int nsrc;
};

#endif
