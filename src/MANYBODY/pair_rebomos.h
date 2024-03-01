/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(rebomos,PairREBOMoS);
// clang-format on
#else

#ifndef LMP_PAIR_REBOMOS_H
#define LMP_PAIR_REBOMOS_H

#include "math_const.h"
#include "pair.h"

#include <cmath>

namespace LAMMPS_NS {

class PairREBOMoS : public Pair {
 public:
  PairREBOMoS(class LAMMPS *);
  ~PairREBOMoS() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  double memory_usage() override;

 protected:
  double **lj1, **lj2, **lj3, **lj4;    // pre-computed LJ coeffs for M,S types
  double cut3rebo;                      // maximum distance for 3rd REBO neigh

  int maxlocal;             // size of numneigh, firstneigh arrays
  int pgsize;               // size of neighbor page
  int oneatom;              // max # of neighbors for one atom
  MyPage<int> *ipage;       // neighbor list pages
  int *REBO_numneigh;       // # of pair neighbors for each atom
  int **REBO_firstneigh;    // ptr to 1st neighbor of each atom

  double *nM, *nS;          // sum of weighting fns with REBO neighs

  double rcmin[2][2], rcmax[2][2], rcmaxsq[2][2], rcmaxp[2][2];
  double Q[2][2], alpha[2][2], A[2][2], BIJc[2][2], Beta[2][2];
  double b0[2], b1[2], b2[2], b3[2], b4[2], b5[2], b6[2];
  double bg0[2], bg1[2], bg2[2], bg3[2], bg4[2], bg5[2], bg6[2];
  double a0[2], a1[2], a2[2], a3[2];
  double rcLJmin[2][2], rcLJmax[2][2];
  double epsilon[2][2], sigma[2][2];

  void REBO_neigh();
  void FREBO(int);
  void FLJ(int);

  double bondorder(int, int, double *, double, double, double **);

  inline double gSpline(const double costh, const int typei, double &dgdc) const
  {
    const double b0i = b0[typei];
    const double b1i = b1[typei];
    const double b2i = b2[typei];
    const double b3i = b3[typei];
    const double b4i = b4[typei];
    const double b5i = b5[typei];
    const double b6i = b6[typei];
    double g = 0.0;

    if (costh >= -1.0 && costh < 0.5) {
      g = b6i * costh;
      double dg = 6.0 * b6i * costh;
      g += b5i;
      dg += 5.0 * b5i;
      g *= costh;
      dg *= costh;
      g += b4i;
      dg += 4.0 * b4i;
      g *= costh;
      dg *= costh;
      g += b3i;
      dg += 3.0 * b3i;
      g *= costh;
      dg *= costh;
      g += b2i;
      dg += 2.0 * b2i;
      g *= costh;
      dg *= costh;
      g += b1i;
      dg += b1i;
      g *= costh;
      g += b0i;
      dgdc = dg;

    } else if (costh >= 0.5 && costh <= 1.0) {
      double gcos = b6i * costh;
      double dgcos = 6.0 * b6i * costh;
      gcos += b5i;
      dgcos += 5.0 * b5i;
      gcos *= costh;
      dgcos *= costh;
      gcos += b4i;
      dgcos += 4.0 * b4i;
      gcos *= costh;
      dgcos *= costh;
      gcos += b3i;
      dgcos += 3.0 * b3i;
      gcos *= costh;
      dgcos *= costh;
      gcos += b2i;
      dgcos += 2.0 * b2i;
      gcos *= costh;
      dgcos *= costh;
      gcos += b1i;
      dgcos += b1i;
      gcos *= costh;
      gcos += b0i;

      const double bg0i = bg0[typei];
      const double bg1i = bg1[typei];
      const double bg2i = bg2[typei];
      const double bg3i = bg3[typei];
      const double bg4i = bg4[typei];
      const double bg5i = bg5[typei];
      const double bg6i = bg6[typei];
      double gamma = bg6i * costh;
      double dgamma = 6.0 * bg6i * costh;
      gamma += bg5i;
      dgamma += 5.0 * bg5i;
      gamma *= costh;
      dgamma *= costh;
      gamma += bg4i;
      dgamma += 4.0 * bg4i;
      gamma *= costh;
      dgamma *= costh;
      gamma += bg3i;
      dgamma += 3.0 * bg3i;
      gamma *= costh;
      dgamma *= costh;
      gamma += bg2i;
      dgamma += 2.0 * bg2i;
      gamma *= costh;
      dgamma *= costh;
      gamma += bg1i;
      dgamma += bg1i;
      gamma *= costh;
      gamma += bg0i;

      const double tmp = MathConst::MY_2PI * (costh - 0.5);
      const double psi = 0.5 * (1 - cos(tmp));
      const double dpsi = MathConst::MY_PI * sin(tmp);
      g = gcos + psi * (gamma - gcos);
      dgdc = dgcos + dpsi * (gamma - gcos) + psi * (dgamma - dgcos);
    } else {
      dgdc = 0.0;
    }
    return g;
  }

  /* ----------------------------------------------------------------------
    Pij calculation
 ------------------------------------------------------------------------- */

  inline double PijSpline(const double NM, const double NS, const int typei, double &dp) const
  {
    const double N = NM + NS;

    dp = -a0[typei] + a1[typei] * a2[typei] * exp(-a2[typei] * N);
    return -a0[typei] * (N - 1) - a1[typei] * exp(-a2[typei] * N) + a3[typei];
  }

  void read_file(char *);
  void allocate();

  // ----------------------------------------------------------------------
  // S'(t) and S(t) cutoff functions
  // added to header for inlining
  // ----------------------------------------------------------------------

  /* ----------------------------------------------------------------------
     cutoff function Sprime
     return cutoff and dX = derivative
     no side effects
  ------------------------------------------------------------------------- */

  inline double Sp(double Xij, double Xmin, double Xmax, double &dX) const
  {
    double cutoff;

    const double t = (Xij - Xmin) / (Xmax - Xmin);
    if (t <= 0.0) {
      cutoff = 1.0;
      dX = 0.0;
    } else if (t >= 1.0) {
      cutoff = 0.0;
      dX = 0.0;
    } else {
      cutoff = 0.5 * (1.0 + cos(t * MathConst::MY_PI));
      dX = (-0.5 * MathConst::MY_PI * sin(t * MathConst::MY_PI)) / (Xmax - Xmin);
    }
    return cutoff;
  };
};
}    // namespace LAMMPS_NS

#endif
#endif
