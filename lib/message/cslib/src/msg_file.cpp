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
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>

#include "msg_file.h"

using namespace CSLIB_NS;

#define MAXLINE 256
#define SLEEP 0.1       // delay in CPU secs to check for message file

/* ---------------------------------------------------------------------- */

MsgFile::MsgFile(int csflag, const void *ptr, MPI_Comm cworld) : 
  Msg(csflag, ptr, cworld)
{
  char *filename = (char *) ptr;
  init(filename);
}

/* ---------------------------------------------------------------------- */

MsgFile::MsgFile(int csflag, const void *ptr) : Msg(csflag, ptr)
{
  char *filename = (char *) ptr;
  init(filename);
}

/* ---------------------------------------------------------------------- */

MsgFile::~MsgFile()
{
  delete [] fileroot;
}

/* ---------------------------------------------------------------------- */

void MsgFile::init(char *filename)
{
  int n = strlen(filename) + 1;
  fileroot = new char[n];
  strcpy(fileroot,filename);
}

/* ---------------------------------------------------------------------- */

void MsgFile::send(int nheader, int *header, int nbuf, char *buf)
{
  char filename[MAXLINE];
  
  lengths[0] = nheader;
  lengths[1] = nbuf;
  
  if (me == 0) {
    if (client) sprintf(filename,"%s.%s",fileroot,"client");
    else if (server) sprintf(filename,"%s.%s",fileroot,"server");
    
    fp = fopen(filename,"wb");
    if (!fp) error_one("send(): Could not open send message file");
    fwrite(lengths,sizeof(int),2,fp);
    fwrite(header,sizeof(int),nheader,fp);
    fwrite(buf,1,nbuf,fp);
    fclose(fp);
  }
  
  // create empty signal file

  if (me == 0) {
    if (client) sprintf(filename,"%s.%s",fileroot,"client.signal");
    else if (server) sprintf(filename,"%s.%s",fileroot,"server.signal");
    fp = fopen(filename,"w");
    fclose(fp);
  }
}

/* ---------------------------------------------------------------------- */

void MsgFile::recv(int &maxheader, int *&header, int &maxbuf, char *&buf)
{
  char filename[MAXLINE];

  // wait until signal file exists to open message file

  if (me == 0) {
    if (client) sprintf(filename,"%s.%s",fileroot,"server.signal");
    else if (server) sprintf(filename,"%s.%s",fileroot,"client.signal");

    int delay = (int) (1000000 * SLEEP);
    while (1) {
      fp = fopen(filename,"r");
      if (fp) break;
      usleep(delay);
    }
    fclose(fp);
  
    if (client) sprintf(filename,"%s.%s",fileroot,"server");
    else if (server) sprintf(filename,"%s.%s",fileroot,"client");
    fp = fopen(filename,"rb");
    if (!fp) error_one("recv(): Could not open recv message file");
  }

  // read and broadcast data
  
  if (me == 0) fread(lengths,sizeof(int),2,fp);
  if (nprocs > 1) MPI_Bcast(lengths,2,MPI_INT,0,world);

  int nheader = lengths[0];
  int nbuf = lengths[1];
  allocate(nheader,maxheader,header,nbuf,maxbuf,buf);
  
  if (me == 0) fread(header,sizeof(int),nheader,fp);
  if (nprocs > 1) MPI_Bcast(header,nheader,MPI_INT,0,world);

  if (me == 0) fread(buf,1,nbuf,fp);
  if (nprocs > 1) MPI_Bcast(buf,nbuf,MPI_CHAR,0,world);

  // delete both message and signal file

  if (me == 0) {
    fclose(fp);
    unlink(filename);
    if (client) sprintf(filename,"%s.%s",fileroot,"server.signal");
    else if (server) sprintf(filename,"%s.%s",fileroot,"client.signal");
    unlink(filename);
  }
}
