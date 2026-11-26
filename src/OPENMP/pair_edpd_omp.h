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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(edpd/omp,PairEDPDOMP);
// clang-format on
#else

#ifndef LMP_PAIR_EDPD_OMP_H
#define LMP_PAIR_EDPD_OMP_H

#include "pair_edpd.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class PairEDPDOMP : public PairEDPD, public ThrOMP {

 public:
  PairEDPDOMP(class LAMMPS *);
  ~PairEDPDOMP() override;

  void settings(int, char **) override;
  void read_restart_settings(FILE *) override;
  void compute(int, int) override;
  double memory_usage() override;

 protected:
  class RanMars **random_thr, **randomT_thr;
  int nthreads;

  void destroy_random_thr();

 private:
  template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
  void eval(int ifrom, int ito, ThrData *const thr);
};

}    // namespace LAMMPS_NS

#endif
#endif
