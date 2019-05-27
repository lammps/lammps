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

#ifndef MSG_MPI_TWO_H
#define MSG_MPI_TWO_H

#include "msg_mpi_one.h"

namespace CSLIB_NS {

class MsgMPITwo : public MsgMPIOne {
 public:
  MsgMPITwo(int, const void *, MPI_Comm);
  ~MsgMPITwo();

 private:
  char port[MPI_MAX_PORT_NAME];

  void init(char *);
};

}

#endif
