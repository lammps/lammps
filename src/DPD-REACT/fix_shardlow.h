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
FixStyle(shardlow,FixShardlow);
// clang-format on
#else

#ifndef LMP_FIX_SHARDLOW_H
#define LMP_FIX_SHARDLOW_H

#include "fix.h"
#include "random_external_state.h"

namespace LAMMPS_NS {

class FixShardlow : public Fix {
 public:
  class NeighList *list;    // The SSA specific neighbor list

  FixShardlow(class LAMMPS *, int, char **);
  ~FixShardlow() override;
  int setmask() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void setup(int) override;
  void initial_integrate(int) override;

  double memory_usage() override;

#ifdef DEBUG_SSA_PAIR_CT
  int counters[2][3];
  int hist[32];
#endif

 protected:
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

  class PairDPDfdt *pairDPD;
  class PairDPDfdtEnergy *pairDPDE;
  double (*v_t0)[3];
  int maxRNG;

 private:
  random_external_state::es_RNG_t *rand_state;
  double dtsqrt;    // = sqrt(update->dt);

  void ssa_update_dpd(int, int, int);     // Constant Temperature
  void ssa_update_dpde(int, int, int);    // Constant Energy
};

}    // namespace LAMMPS_NS

#endif
#endif
