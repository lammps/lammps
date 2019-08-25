#ifndef MANY2MANY_H
#define MANY2MANY_H

#include <mpi.h>

class Many2Many {
 public:
  Many2Many(MPI_Comm);
  ~Many2Many();

  void setup(int, int *, int, int *);
  void exchange(int *, int *);
  void exchange(double *, double *);

 protected:
  int me,nprocs;
  MPI_Comm comm;
  class Memory *memory;
  class Error *error;

  int nown;                       // # of IDs common to src and dest
  int nsrc_off,ndest_off;         // # of off-processor IDs

  int *src_own,*dest_own;         // indices of the owned IDs
  int *src_off,*dest_off;         // indices of the off-proc IDs

  int *src_iwork,*dest_iwork;     // work arrays for comm of ints
  double *src_dwork,*dest_dwork;  // work arrays for comm of doubles

  class Irregular *irregular;     // irregular comm from src->dest

  struct Datum1 {
    int id;            // src or dest global ID
    int proc;          // owning proc
    int index;         // local index on owning proc
  };

  struct Datum2 {
    int slocal;        // local index of src ID on sending proc
    int dlocal;        // local index of dest ID on receiving proc
    int dproc;         // receiving proc
  };

  void deallocate();
};

#endif
