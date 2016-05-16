#ifndef SEND2ONE_H
#define SEND2ONE_H

#include <mpi.h>

class Send2One {
 public:
  Send2One(MPI_Comm);
  virtual ~Send2One();

  void execute();

 protected:
  int me,nprocs;
  MPI_Comm comm;
  class Memory *memory;
  class Error *error;

  int maxbuf;
  char *buf;

  virtual void pre() = 0;
  virtual int size() = 0;
  virtual void pack(char *) = 0;
  virtual void process(int, char *) = 0;
  virtual void post() = 0;
};

#endif
