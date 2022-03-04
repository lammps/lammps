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
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <unistd.h>

#include "msg_mpi_one.h"

using namespace CSLIB_NS;

/* ---------------------------------------------------------------------- */

MsgMPIOne::MsgMPIOne(int csflag, const void *ptr, MPI_Comm cworld) :
  Msg(csflag, ptr, cworld)
{
  // NOTE: ideally would skip this call if mpi/two
  init(ptr);
}

/* ---------------------------------------------------------------------- */

void MsgMPIOne::init(const void *ptr)
{
  MPI_Comm *pbothcomm = (MPI_Comm *) ptr;
  bothcomm = *pbothcomm;

  if (client) {
    MPI_Comm_size(world,&nprocs);
    otherroot = nprocs;
  } else if (server) {
    otherroot = 0;
  }
}

/* ---------------------------------------------------------------------- */

void MsgMPIOne::send(int nheader, int *header, int nbuf, char *buf)
{
  lengths[0] = nheader;
  lengths[1] = nbuf;

  if (me == 0) {
    MPI_Send(lengths,2,MPI_INT,otherroot,0,bothcomm);
    MPI_Send(header,nheader,MPI_INT,otherroot,0,bothcomm);
    MPI_Send(buf,nbuf,MPI_CHAR,otherroot,0,bothcomm);
  }
}

/* ---------------------------------------------------------------------- */

void MsgMPIOne::recv(int &maxheader, int *&header, int &maxbuf, char *&buf)
{
  MPI_Status status;

  if (me == 0) MPI_Recv(lengths,2,MPI_INT,otherroot,0,bothcomm,&status);
  if (nprocs > 1) MPI_Bcast(lengths,2,MPI_INT,0,world);

  int nheader = lengths[0];
  int nbuf = lengths[1];
  allocate(nheader,maxheader,header,nbuf,maxbuf,buf);

  if (me == 0) MPI_Recv(header,nheader,MPI_INT,otherroot,0,bothcomm,&status);
  if (nprocs > 1) MPI_Bcast(header,nheader,MPI_INT,0,world);

  if (me == 0) MPI_Recv(buf,nbuf,MPI_CHAR,otherroot,0,bothcomm,&status);
  if (nprocs > 1) MPI_Bcast(buf,nbuf,MPI_CHAR,0,world);
}
