#ifndef ERRORLIB_H
#define ERRORLIB_H

#include <mpi.h>

class ErrorLib {
 public:
  ErrorLib(MPI_Comm);

  void all(const char *);
  void one(const char *);
  void warning(const char *);

 private:
  MPI_Comm comm;
  int me;
};

#endif
