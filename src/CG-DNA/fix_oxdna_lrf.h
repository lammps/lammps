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
FixStyle(OXDNA/LRF,FixOxdnaLRF);
// clang-format on
#else

#ifndef LMP_FIX_OXDNA_LRF_H
#define LMP_FIX_OXDNA_LRF_H

#include "fix.h"

namespace LAMMPS_NS {

class FixOxdnaLRF : public Fix {
 public:
  FixOxdnaLRF(class LAMMPS *, int, char **);
  ~FixOxdnaLRF() override;
  int setmask() override;
  void min_setup_pre_force(int);
  void min_pre_force(int) override;
  void setup_pre_force(int) override;
  void pre_force(int) override;

  double memory_usage() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  void set_arrays(int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

 private:
  double **nxyz;    // per-atom arrays for local unit vectors
  void compute_lrf();

  class AtomVecEllipsoid *avec;
};

}    // namespace LAMMPS_NS
#endif
#endif
