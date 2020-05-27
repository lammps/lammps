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

#ifndef MSG_FILE_H
#define MSG_FILE_H

#include <stdio.h>
#include "msg.h"

namespace CSLIB_NS {

class MsgFile : public Msg {
 public:
  MsgFile(int, const void *, MPI_Comm);
  MsgFile(int, const void *);
  ~MsgFile();
  void send(int, int *, int, char *);
  void recv(int &, int *&, int &, char *&);

 private:
  char *fileroot;
  FILE *fp;

  void init(char *);
};

}

#endif
