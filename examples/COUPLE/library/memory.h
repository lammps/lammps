#ifndef MEMORY_H
#define MEMORY_H

#include "mpi.h"

class Memory {
 public:
  Memory(MPI_Comm);
  ~Memory();

  void *smalloc(int n, const char *);
  void sfree(void *);
  void *srealloc(void *, int n, const char *name);

  double **create_2d_double_array(int, int, const char *);
  double **grow_2d_double_array(double **, int, int, const char *);
  void destroy_2d_double_array(double **);

 private:
  class Error *error;
};

#endif
