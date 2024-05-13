#ifndef MEMORY_LIB_H
#define MEMORY_LIB_H

#include <mpi.h>

class MemoryLib {
 public:
  MemoryLib(MPI_Comm);
  ~MemoryLib();

  void *smalloc(int n, const char *);
  void sfree(void *);
  void *srealloc(void *, int n, const char *name);

  double **create_2d_double_array(int, int, const char *);
  double **grow_2d_double_array(double **, int, int, const char *);
  void destroy_2d_double_array(double **);

 private:
  class ErrorLib *error;
};

#endif
