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

#include "msg_mpi_two.h"

using namespace CSLIB_NS;

/* ---------------------------------------------------------------------- */

MsgMPITwo::MsgMPITwo(int csflag, const void *ptr, MPI_Comm cworld) :
  MsgMPIOne(csflag, ptr, cworld)
{
  char *filename = (char *) ptr;
  init(filename);
}

/* ---------------------------------------------------------------------- */

MsgMPITwo::~MsgMPITwo()
{
  // free the inter comm that spans both client and server

  MPI_Comm_free(&bothcomm);
  MPI_Close_port(port);
}

/* ---------------------------------------------------------------------- */

void MsgMPITwo::init(char *filename)
{
  if (client) {
    if (me == 0) {
      FILE *fp = nullptr;
      while (!fp) {
        fp = fopen(filename,"r");
        if (!fp) sleep(1);
      }
      fgets(port,MPI_MAX_PORT_NAME,fp);
      //printf("Client port: %s\n",port);
      fclose(fp);
    }

    MPI_Bcast(port,MPI_MAX_PORT_NAME,MPI_CHAR,0,world);
    MPI_Comm_connect(port,MPI_INFO_NULL,0,world,&bothcomm);
    //if (me == 0) printf("CLIENT comm connect\n");
    if (me == 0) unlink(filename);

  } else if (server) {
    MPI_Open_port(MPI_INFO_NULL,port);

    if (me == 0) {
      //printf("Server name: %s\n",port);
      FILE *fp = fopen(filename,"w");
      fprintf(fp,"%s",port);
      fclose(fp);
    }

    MPI_Comm_accept(port,MPI_INFO_NULL,0,world,&bothcomm);
    //if (me == 0) printf("SERVER comm accept\n");
  }

  otherroot = 0;
}
