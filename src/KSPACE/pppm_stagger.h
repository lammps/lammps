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

#ifdef KSPACE_CLASS
// clang-format off
KSpaceStyle(pppm/stagger,PPPMStagger);
// clang-format on
#else

#ifndef LMP_PPPM_STAGGER_H
#define LMP_PPPM_STAGGER_H

#include "pppm.h"

namespace LAMMPS_NS {

class PPPMStagger : public PPPM {
 public:
  PPPMStagger(class LAMMPS *);
  ~PPPMStagger() override;
  void init() override;
  void compute(int, int) override;
  int timing_1d(int, double &) override;
  int timing_3d(int, double &) override;

 protected:
  int nstagger;
  double stagger;
  double **gf_b2;

  double compute_qopt() override;
  double compute_qopt_ad();
  void compute_gf_denom() override;
  void compute_gf_ik() override;
  void compute_gf_ad() override;

  void particle_map() override;
  void make_rho() override;
  void fieldforce_ik() override;
  void fieldforce_ad() override;
  void fieldforce_peratom() override;

  inline double gf_denom2(const double &x, const double &y, const double &z) const
  {
    double sx, sy, sz;
    double x2 = x * x;
    double y2 = y * y;
    double z2 = z * z;
    double xl = x;
    double yl = y;
    double zl = z;
    sx = sy = sz = 0.0;
    for (int l = 0; l < order; l++) {
      sx += gf_b2[order][l] * xl;
      sy += gf_b2[order][l] * yl;
      sz += gf_b2[order][l] * zl;
      xl *= x2;
      yl *= y2;
      zl *= z2;
    }
    double s = sx * sy * sz;
    return s * s;
  };
};

}    // namespace LAMMPS_NS

#endif
#endif
