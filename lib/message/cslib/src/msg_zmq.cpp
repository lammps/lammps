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

#include <mpi.h>
#include <zmq.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

#include "msg_zmq.h"

using namespace CSLIB_NS;

/* ---------------------------------------------------------------------- */

MsgZMQ::MsgZMQ(int csflag, const void *ptr, MPI_Comm cworld) :
  Msg(csflag, ptr, cworld)
{
  char *port = (char *) ptr;
  init(port);
}

MsgZMQ::MsgZMQ(int csflag, const void *ptr) : Msg(csflag, ptr)
{
  char *port = (char *) ptr;
  init(port);
}

/* ---------------------------------------------------------------------- */

MsgZMQ::~MsgZMQ()
{
  if (me == 0) {
    zmq_close(socket);
    zmq_ctx_destroy(context);
  }
}

/* ---------------------------------------------------------------------- */

void MsgZMQ::init(char *port)
{
#ifdef ZMQ_NO
  error_all("constructor(): Library not built with ZMQ support");
#endif

  if (me == 0) {
    int n = strlen(port) + 8;
    char *socket_name = new char[n];
    strcpy(socket_name,"tcp://");
    strcat(socket_name,port);
  
    if (client) {
      context = zmq_ctx_new();
      socket = zmq_socket(context,ZMQ_REQ);
      zmq_connect(socket,socket_name);
    } else if (server) {
      context = zmq_ctx_new();
      socket = zmq_socket(context,ZMQ_REP);
      int rc = zmq_bind(socket,socket_name);
      if (rc) error_one("constructor(): Server could not make socket connection");
    }

    delete [] socket_name;
  }
}

/* ----------------------------------------------------------------------
   client/server sockets (REQ/REP) must follow this protocol:
     client sends request (REQ) which server receives
     server sends response (REP) which client receives
     every exchange is of this form, server cannot initiate a send
   thus each ZMQ send below has a following ZMQ recv, except last one
     if client calls send(), it will next call recv()
     if server calls send(), it will next call recv() from its wait loop
     in either case, recv() issues a ZMQ recv to match last ZMQ send here
------------------------------------------------------------------------- */

void MsgZMQ::send(int nheader, int *header, int nbuf, char *buf)
{
  lengths[0] = nheader;
  lengths[1] = nbuf;

  if (me == 0) {
    zmq_send(socket,lengths,2*sizeof(int),0);
    zmq_recv(socket,NULL,0,0);
  }

  if (me == 0) {
    zmq_send(socket,header,nheader*sizeof(int),0);
    zmq_recv(socket,NULL,0,0);
  }

  if (me == 0) zmq_send(socket,buf,nbuf,0);
}

/* ----------------------------------------------------------------------
   client/server sockets (REQ/REP) must follow this protocol:
     client sends request (REQ) which server receives
     server sends response (REP) which client receives
     every exchange is of this form, server cannot initiate a send
   thus each ZMQ recv below has a following ZMQ send, except last one
     if client calls recv(), it will next call send() to ping server again,
     if server calls recv(), it will next call send() to respond to client
     in either case, send() issues a ZMQ send to match last ZMQ recv here
------------------------------------------------------------------------- */

void MsgZMQ::recv(int &maxheader, int *&header, int &maxbuf, char *&buf)
{
  if (me == 0) {
    zmq_recv(socket,lengths,2*sizeof(int),0);
    zmq_send(socket,NULL,0,0);
  }
  if (nprocs > 1) MPI_Bcast(lengths,2,MPI_INT,0,world);

  int nheader = lengths[0];
  int nbuf = lengths[1];
  allocate(nheader,maxheader,header,nbuf,maxbuf,buf);

  if (me == 0) {
    zmq_recv(socket,header,nheader*sizeof(int),0);
    zmq_send(socket,NULL,0,0);
  }
  if (nprocs > 1) MPI_Bcast(header,nheader,MPI_INT,0,world);

  if (me == 0) zmq_recv(socket,buf,nbuf,0);
  if (nprocs > 1) MPI_Bcast(buf,nbuf,MPI_CHAR,0,world);
}
