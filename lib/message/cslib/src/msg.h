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

#ifndef MSG_H
#define MSG_H

#ifdef MPI_YES
#include <mpi.h>
#else
#include <mpi_dummy.h>
#endif

namespace CSLIB_NS {

class Msg {
 public:
  int nsend,nrecv;
  MPI_Comm world;

  Msg(int, const void *, MPI_Comm);
  Msg(int, const void *);
  virtual ~Msg() {}
  virtual void send(int, int *, int, char *) = 0;
  virtual void recv(int &, int *&, int &, char *&) = 0;

 protected:
  int me,nprocs;
  int client,server;

  int nfield;
  int *fieldID,*fieldtype,*fieldlen;
  int lengths[2];

  void init(int);
  void allocate(int, int &, int *&, int, int &, char *&);
  void *smalloc(int);
  void sfree(void *);
  void error_all(const char *);
  void error_one(const char *);
};


}

#endif
