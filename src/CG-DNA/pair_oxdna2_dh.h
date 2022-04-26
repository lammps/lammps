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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(oxdna2/dh,PairOxdna2Dh);
// clang-format on
#else

#ifndef LMP_PAIR_OXDNA2_DH_H
#define LMP_PAIR_OXDNA2_DH_H

#include "pair.h"

namespace LAMMPS_NS {

class PairOxdna2Dh : public Pair {
 public:
  PairOxdna2Dh(class LAMMPS *);
  ~PairOxdna2Dh() override;
  virtual void compute_interaction_sites(double *, double *, double *, double *);
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_list(int, class NeighList *) override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  void write_data(FILE *) override;
  void write_data_all(FILE *) override;
  void *extract(const char *, int &) override;

 protected:
  double **qeff_dh_pf, **kappa_dh;
  double **b_dh, **cut_dh_ast, **cutsq_dh_ast, **cut_dh_c, **cutsq_dh_c;
  double **nx_xtrct, **ny_xtrct, **nz_xtrct;    // per-atom arrays for local unit vectors

  virtual void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
