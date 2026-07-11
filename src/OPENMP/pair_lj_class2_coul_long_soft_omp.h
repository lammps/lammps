/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(lj/class2/coul/long/soft/omp,PairLJClass2CoulLongSoftOMP);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_CLASS2_COUL_LONG_SOFT_OMP_H
#define LMP_PAIR_LJ_CLASS2_COUL_LONG_SOFT_OMP_H

#include "pair_lj_class2_coul_long_soft.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class PairLJClass2CoulLongSoftOMP : public PairLJClass2CoulLongSoft, public ThrOMP {

 public:
  PairLJClass2CoulLongSoftOMP(class LAMMPS *);

  void compute(int, int) override;
  double memory_usage() override;

 private:
  template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
  void eval(int ifrom, int ito, ThrData *const thr);
};

}    // namespace LAMMPS_NS

#endif
#endif
