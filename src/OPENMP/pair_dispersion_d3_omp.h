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
PairStyle(dispersion/d3/omp,PairDispersionD3OMP);
// clang-format on
#else

#ifndef LMP_PAIR_DISPERSION_D3_OMP_H
#define LMP_PAIR_DISPERSION_D3_OMP_H

#include "pair_dispersion_d3.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class PairDispersionD3OMP : public PairDispersionD3, public ThrOMP {

 public:
  PairDispersionD3OMP(class LAMMPS *);

  void compute(int, int) override;
  double memory_usage() override;

 private:
  template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
  void eval_first_phase(int iifrom, int iito, ThrData *const thr);

  template <int NEWTON_PAIR> void eval_coordination(int iifrom, int iito, ThrData *const thr);

  void calc_coordination_number() override;

  template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
  void eval_second_phase(int ifrom, int ito, ThrData *const thr);
};

}    // namespace LAMMPS_NS

#endif
#endif
