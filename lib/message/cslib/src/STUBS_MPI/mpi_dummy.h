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

// MPI constants and dummy functions

#ifndef MPI_DUMMY_H
#define MPI_DUMMY_H

#include <stdlib.h>
#include <stdint.h>
#include <string.h>

namespace CSLIB_NS {

typedef int MPI_Comm;
typedef int MPI_Fint;
typedef int MPI_Datatype;
typedef int MPI_Status;
typedef int MPI_Op;
typedef int MPI_Info;

#define MPI_COMM_WORLD 0
#define MPI_MAX_PORT_NAME 0
#define MPI_INFO_NULL 0
#define MPI_INT 1
#define MPI_LONG_LONG 2
#define MPI_FLOAT 3
#define MPI_DOUBLE 4
#define MPI_CHAR 5
#define MPI_SUM 0

static void MPI_Init(int *, char ***) {}
static MPI_Comm MPI_Comm_f2c(MPI_Comm world) {return world;}
static void MPI_Comm_rank(MPI_Comm, int *) {}
static void MPI_Comm_size(MPI_Comm, int *) {}

static void MPI_Open_port(MPI_Info, char *) {}
static void MPI_Close_port(const char *) {}
static void MPI_Comm_accept(const char *, MPI_Info, int,
                            MPI_Comm, MPI_Comm *) {}
static void MPI_Comm_connect(const char *, MPI_Info, int,
                             MPI_Comm, MPI_Comm *) {}

static void MPI_Comm_split(MPI_Comm, int, int, MPI_Comm *) {}
static void MPI_Comm_free(MPI_Comm *) {}

static void MPI_Send(const void *, int, MPI_Datatype, int, int, MPI_Comm) {}
static void MPI_Recv(void *, int, MPI_Datatype, int, int,
                     MPI_Comm, MPI_Status *) {}

static void MPI_Allreduce(const void *in, void *out, int, MPI_Datatype type,
                          MPI_Op op, MPI_Comm)
{
  if (type == MPI_INT) *((int *) out) = *((int *) in);
}
static void MPI_Scan(const void *in, void *out, int, MPI_Datatype intype,
                     MPI_Op op,MPI_Comm)
{
  if (intype == MPI_INT) *((int *) out) = *((int *) in);
}

static void MPI_Bcast(void *, int, MPI_Datatype, int, MPI_Comm) {}
static void MPI_Allgather(const void *in, int incount, MPI_Datatype intype,
                          void *out, int, MPI_Datatype, MPI_Comm)
{
  // assuming incount = 1
  if (intype == MPI_INT) *((int *) out) = *((int *) in);
}
static void MPI_Allgatherv(const void *in, int incount, MPI_Datatype intype,
                           void *out, const int *, const int *,
                           MPI_Datatype, MPI_Comm)
{
  if (intype == MPI_INT) memcpy(out,in,incount*sizeof(int));
  else if (intype == MPI_LONG_LONG) memcpy(out,in,incount*sizeof(int64_t));
  else if (intype == MPI_FLOAT) memcpy(out,in,incount*sizeof(float));
  else if (intype == MPI_DOUBLE) memcpy(out,in,incount*sizeof(double));
  else if (intype == MPI_CHAR) memcpy(out,in,incount*sizeof(char));
}

static void MPI_Abort(MPI_Comm, int) {exit(1);}
static void MPI_Finalize() {}

}

#endif
