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

/* ----------------------------------------------------------------------
   Contributing author: Trung Dac Nguyen (ndactrung@gmail.com)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(body/rounded/polygon/omp,PairBodyRoundedPolygonOMP);
// clang-format on
#else

#ifndef LMP_PAIR_BODY_ROUNDED_POLYGON_OMP_H
#define LMP_PAIR_BODY_ROUNDED_POLYGON_OMP_H

#include "pair_body_rounded_polygon.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class PairBodyRoundedPolygonOMP : public PairBodyRoundedPolygon, public ThrOMP {

 public:
  PairBodyRoundedPolygonOMP(class LAMMPS *);
  ~PairBodyRoundedPolygonOMP() override;

  void compute(int, int) override;
  double memory_usage() override;

 private:
  double **fnc_thr;                     // per-thread fnc, see PairBodyRoundedPolygon
  int nmax_thr;                         // allocated number of atoms per thread in fnc_thr
  std::vector<Scratch> scratch_thr;    // per-thread scratch space

  template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
  void eval(int ifrom, int ito, ThrData *const thr, double **fnc_t, Scratch &s);
};

}    // namespace LAMMPS_NS

#endif
#endif
