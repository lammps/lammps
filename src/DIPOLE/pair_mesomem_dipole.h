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
PairStyle(mesomem/dipole,PairMesomemDipole);
// clang-format on
#else

#ifndef LMP_PAIR_MESOMEM_DIPOLE_H
#define LMP_PAIR_MESOMEM_DIPOLE_H

#include "pair.h"

namespace LAMMPS_NS {

class PairMesomemDipole : public Pair {
 public:
  PairMesomemDipole(class LAMMPS *);
  ~PairMesomemDipole() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  void write_data(FILE *) override;
  void write_data_all(FILE *) override;

 protected:
  double **cut, **sigma, **eps;
  double **ktilt, **ksplay;
  double **weight_rcut;
  double **zeta, **c0;

  // per type pair constants derived from the coefficients in init_one()
  double **gscale;         // pi/2 / (r_c - sigma)
  double **wc_inv;         // 1 / w_c
  double **wc_half2inv;    // 1 / (w_c/2)^2
  int **zpow;              // 2*zeta - 1 if it is a small non-negative integer, otherwise -1

  virtual void allocate();
  double mesomem_analytic(int, int, double, const double *, const double *, const double *,
                          double *, double *, double *) const;
};

}    // namespace LAMMPS_NS

#endif
#endif
