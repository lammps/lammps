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

#ifndef MSG_ZMQ_H
#define MSG_ZMQ_H

#include "msg.h"

namespace CSLIB_NS {

class MsgZMQ : public Msg {
 public:
  MsgZMQ(int, const void *, MPI_Comm);
  MsgZMQ(int, const void *);
  ~MsgZMQ();
  void send(int, int *, int, char *);
  void recv(int &, int *&, int &, char *&);

 private:
  void *context,*socket;

  void init(char *);
};

}

#endif
