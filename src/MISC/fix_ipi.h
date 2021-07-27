/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(ipi,FixIPI);
// clang-format on
#else

#ifndef LMP_FIX_IPI_H
#define LMP_FIX_IPI_H

#include "fix.h"

namespace LAMMPS_NS {

class FixIPI : public Fix {
 public:
  FixIPI(class LAMMPS *, int, char **);
  virtual ~FixIPI();
  int setmask();
  virtual void init();
  virtual void initial_integrate(int);
  virtual void final_integrate();

 protected:
  char *host;
  int port;
  int inet, master, hasdata;
  int ipisock, me, socketflag;
  double *buffer;
  long bsize;
  int kspace_flag;
  int reset_flag;

 private:
  class Irregular *irregular;
};

}    // namespace LAMMPS_NS

#endif
#endif
