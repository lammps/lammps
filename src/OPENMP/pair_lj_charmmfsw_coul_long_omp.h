/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(lj/charmmfsw/coul/long/omp,PairLJCharmmfswCoulLongOMP);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_CHARMMFSW_COUL_LONG_OMP_H
#define LMP_PAIR_LJ_CHARMMFSW_COUL_LONG_OMP_H

#include "pair_lj_charmmfsw_coul_long.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class PairLJCharmmfswCoulLongOMP : public PairLJCharmmfswCoulLong, public ThrOMP {

 public:
  PairLJCharmmfswCoulLongOMP(class LAMMPS *);

  void compute(int, int) override;
  double memory_usage() override;

 private:
  template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
  void eval(int ifrom, int ito, ThrData *const thr);
};

}    // namespace LAMMPS_NS

#endif
#endif
