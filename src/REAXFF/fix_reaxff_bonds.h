// clang-format off
/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(reaxff/bonds,FixReaxFFBonds);
// clang-format on
#else

#ifndef LMP_FIX_REAXC_BONDS_H
#define LMP_FIX_REAXC_BONDS_H

#include "fix.h"

namespace LAMMPS_NS {

class FixReaxFFBonds : public Fix {
 public:
  FixReaxFFBonds(class LAMMPS *, int, char **);
  ~FixReaxFFBonds() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void end_of_step() override;

 protected:
  int me, nprocs, nmax, ntypes, maxsize, compressed;
  int *numneigh;
  tagint **neighid;
  double **abo;
  FILE *fp;

  void allocate();
  void destroy();
  virtual void Output_ReaxFF_Bonds();
  int FindBond();
  void PassBuffer(double *, int &);
  void RecvBuffer(double *, int, int, int, int);
  int nint(const double &);
  double memory_usage() override;

  bigint nvalid, nextvalid();
  struct _reax_list *lists;
  class PairReaxFF *reaxff;
  class NeighList *list;
};
}    // namespace LAMMPS_NS

#endif
#endif
