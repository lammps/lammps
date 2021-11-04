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

#ifndef LMP_PAIR_PERI_H
#define LMP_PAIR_PERI_H

#include "pair.h"
#include <cmath>

namespace LAMMPS_NS {

class PairPeri : public Pair {
 public:
  PairPeri(class LAMMPS *);
  virtual ~PairPeri();

  virtual int pack_forward_comm(int, int *, double *, int, int *);
  virtual void unpack_forward_comm(int, int, double *);

  virtual void init_style();
  virtual void settings(int, char **);
  static constexpr double NEAR_ZERO = 2.2204e-16;

  double influence_function(const double &xi_x, const double &xi_y, const double &xi_z) const
  {
    const double r = sqrt((xi_x * xi_x) + (xi_y * xi_y) + (xi_z * xi_z));
    return (fabs(r) < NEAR_ZERO) ? 1.0 / NEAR_ZERO : (1.0 / r);
  }
  void compute_dilatation(int, int);

  double memory_usage();
  virtual void *extract(const char *, int &);

 protected:
  class FixPeriNeigh *fix_peri_neigh;
  double **bulkmodulus, **shearmodulus, **m_lambdai, **m_taubi, **m_yieldstress;
  double **s00, **alpha, **cut, **kspring;
  double *s0_new, *theta, *elastic_energy;

  int nmax;

 protected:
  void allocate();
};

}    // namespace LAMMPS_NS

#endif
