/* ----------------------------------------------------------------------
   CSlib - Client/server library for code coupling
   https://cslib.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright 2018 National Technology & Engineering Solutions of
   Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 with
   NTESS, the U.S. Government retains certain rights in this software.
   This software is distributed under the modified Berkeley Software
   Distribution (BSD) License.

   See the README file in the top-level CSlib directory.
------------------------------------------------------------------------- */

#ifdef MPI_YES
#include <mpi.h>
#else
#include <mpi_dummy.h>
#endif
#include <cstdio>
#include <cstring>
#include <cstdlib>

#include "msg.h"

using namespace CSLIB_NS;

/* ---------------------------------------------------------------------- */

Msg::Msg(int csflag, const void * /* ptr */, MPI_Comm cworld)
{
  world = cworld;
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  init(csflag);
}

/* ---------------------------------------------------------------------- */

Msg::Msg(int csflag, const void * /* ptr */)
{
  world = 0;
  me = 0;
  nprocs = 1;

  init(csflag);
}

/* ---------------------------------------------------------------------- */

void Msg::init(int csflag)
{
  client = server = 0;
  if (csflag == 0) client = 1;
  else if (csflag == 1) server = 1;

  nsend = nrecv = 0;
}

/* ---------------------------------------------------------------------- */

void Msg::allocate(int nheader, int &maxheader, int *&header,
                   int nbuf, int &maxbuf, char *&buf)
{
  if (nheader > maxheader) {
    sfree(header);
    maxheader = nheader;
    header = (int *) smalloc(maxheader*sizeof(int));
  }

  if (nbuf > maxbuf) {
    sfree(buf);
    maxbuf = nbuf;
    buf = (char *) smalloc(maxbuf*sizeof(char));
  }
}

/* ---------------------------------------------------------------------- */

void *Msg::smalloc(int nbytes)
{
  if (nbytes == 0) return nullptr;
  void *ptr = (void *) malloc(nbytes);
  if (ptr == nullptr) {
    char str[128];
    sprintf(str,"Failed to allocate %d bytes",nbytes);
  }
  return ptr;
}

/* ---------------------------------------------------------------------- */

void Msg::sfree(void *ptr)
{
  if (ptr == nullptr) return;
  free(ptr);
}

/* ---------------------------------------------------------------------- */

void Msg::error_all(const char *str)
{
  if (me == 0) printf("CSlib ERROR: %s\n",str);
  MPI_Abort(world,1);
}

/* ---------------------------------------------------------------------- */

void Msg::error_one(const char *str)
{
  printf("CSlib ERROR: %s\n",str);
  MPI_Abort(world,1);
}
