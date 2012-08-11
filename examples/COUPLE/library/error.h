#ifndef ERROR_H
#define ERROR_H

#include "mpi.h"

class Error {
 public:
  Error(MPI_Comm);

  void all(const char *);
  void one(const char *);
  void warning(const char *);

 private:
  MPI_Comm comm;
  int me;
};

#endif
