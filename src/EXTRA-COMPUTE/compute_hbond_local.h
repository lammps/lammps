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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(hbond/local,ComputeHBondLocal);
// clang-format on
#else

#ifndef LMP_COMPUTE_HBOND_LOCAL_H
#define LMP_COMPUTE_HBOND_LOCAL_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeHBondLocal : public Compute {
 public:
  ComputeHBondLocal(class LAMMPS *, int, char **);
  ~ComputeHBondLocal() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  double compute_scalar() override;
  int compute_image(int *&, double **&) override;
  void compute_local() override;
  double memory_usage() override;

 private:
  double distcutoff, anglecutoff, distcutoffsq, ehbcutoff;
  int ncount;
  int singleflag, hdistflag;
  int hydrogenmask, donormask, acceptormask;    // group bitmask for groups of atoms
  int nmax;
  std::vector<int> vflag;    // value flags
  double **alocal;           // local array storage
  class NeighList *list;     // full neighbor list

  // arrays for dump image rendering

  int numobjs;
  int *imgobjs;
  double **imgparms;

  int compute_hbonds(int);
  void reallocate(int);
};
}    // namespace LAMMPS_NS
#endif
#endif
