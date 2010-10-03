#ifndef MANY2ONE_H
#define MANY2ONE_H

#include "mpi.h"

class Many2One {
 public:
  Many2One(MPI_Comm);
  ~Many2One();

  void setup(int, int *, int);
  void gather(double *, int, double *);

 protected:
  int me,nprocs;
  MPI_Comm comm;
  class Memory *memory;

  int nsrc,nall;
  int *counts,*multicounts;
  int *displs,*multidispls;
  int *idall;
};

#endif
