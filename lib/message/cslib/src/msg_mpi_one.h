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

#ifndef MSG_MPI_ONE_H
#define MSG_MPI_ONE_H

#include "msg.h"

namespace CSLIB_NS {

class MsgMPIOne : public Msg {
 public:
  MsgMPIOne(int, const void *, MPI_Comm);
  virtual ~MsgMPIOne() {}
  void send(int, int *, int, char *);
  void recv(int &, int *&, int &, char *&);

 protected:
  MPI_Comm bothcomm;
  int otherroot;

  void init(const void *);
};

}

#endif
