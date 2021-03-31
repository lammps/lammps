/* ----------------------------------------------------------------------
   CSlib - Client/server library for code coupling
   http://cslib.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright 2018 National Technology & Engineering Solutions of
   Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
   NTESS, the U.S. Government retains certain rights in this software.
   This software is distributed under the modified Berkeley Software
   Distribution (BSD) License.

   See the README file in the top-level CSlib directory.
------------------------------------------------------------------------- */

#ifndef CSLIB_H
#define CSLIB_H

#include <stdint.h>

#if defined(LAMMPS_BIGBIG)
#error CSlib is not compatible with -DLAMMPS_BIGBIG
#endif

namespace CSLIB_NS {

class CSlib {
 public:
  int nsend,nrecv;

  CSlib(int, const char *, const void *, const void *);
  ~CSlib();

  void send(int, int);

  void pack_int(int, int);
  void pack_int64(int, int64_t);
  void pack_float(int, float);
  void pack_double(int, double);
  void pack_string(int, char *);
  void pack(int, int, int, void *);
  void pack_parallel(int, int, int, int *, int, void *);

  int recv(int &, int *&, int *&, int *&);

  int unpack_int(int);
  int64_t unpack_int64(int);
  float unpack_float(int);
  double unpack_double(int);
  char *unpack_string(int);
  void *unpack(int);
  void unpack(int, void *);
  void unpack_parallel(int, int, int *, int, void *);

  int extract(int);
  
 private:
  uint64_t myworld;    // really MPI_Comm, but avoids use of mpi.h in this file
                       // so apps can include this file w/ no MPI on system
  int me,nprocs;
  int client,server;
  int nfield,maxfield;
  int msgID,fieldcount;
  int nheader,maxheader;
  int nbuf,maxbuf;
  int maxglobal,maxfieldbytes;
  int *fieldID,*fieldtype,*fieldlen,*fieldoffset;
  int *header;
  int *recvcounts,*displs;    // nprocs size for Allgathers
  int *allids;                // nglobal size for pack_parallel()
  char *buf;                  // maxbuf size for msg with all fields
  char *fielddata;            // maxfieldbytes size for one global field
  const char *pad;

  class Msg *msg;

  void send_message();
  void onefield(int, int, int &, int &);
  int find_field(int, int);
  void allocate_fields();
  void deallocate_fields();
  int64_t roundup(int64_t, int);
  void *smalloc(int);
  void *srealloc(void *, int);
  void sfree(void *);
  void error_all(const char *);
  void error_one(const char *);
};

}

#endif
