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
PairStyle(oxdna/excv,PairOxdnaExcv);
// clang-format on
#else

#ifndef LMP_PAIR_OXDNA_EXCV_H
#define LMP_PAIR_OXDNA_EXCV_H

#include "pair.h"

namespace LAMMPS_NS {

class PairOxdnaExcv : public Pair {
 public:
  PairOxdnaExcv(class LAMMPS *);
  ~PairOxdnaExcv() override;
  virtual void compute_backbone_site(double *, double *, double *, double *) const;
  virtual void compute_base_site(int, double *, double *, double *, double *) const;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  void init_list(int, class NeighList *) override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  void *extract(const char *, int &) override;

 protected:
  // s=sugar-phosphate backbone site, b=base site, st=stacking site

  // excluded volume interaction
  // base step-dependent coefficients
  double **epsilon_bkbk, **sigma_bkbk, **cut_bkbk_ast, **cutsq_bkbk_ast;
  double **lj1_bkbk, **lj2_bkbk, **b_bkbk, **cut_bkbk_c, **cutsq_bkbk_c;
  double **epsilon_bkbs, **sigma_bkbs, **cut_bkbs_ast, **cutsq_bkbs_ast;
  double **lj1_bkbs, **lj2_bkbs, **b_bkbs, **cut_bkbs_c, **cutsq_bkbs_c;
  double **epsilon_bsbs, **sigma_bsbs, **cut_bsbs_ast, **cutsq_bsbs_ast;
  double **lj1_bsbs, **lj2_bsbs, **b_bsbs, **cut_bsbs_c, **cutsq_bsbs_c;
  // tetramer-dependent coefficients
  double ****sigma4_bsbs, ****cut4_bsbs_ast, ****cut4sq_bsbs_ast;
  double ****lj14_bsbs, ****lj24_bsbs, ****b4_bsbs, ****cut4_bsbs_c, ****cut4sq_bsbs_c;

  double **nxyz_xtrct;    // per-atom arrays for local unit vectors
  virtual void allocate();

  class FixOxdnaLRF *fix_lrf;    // ptr to oxdna/lrf fix
};

}    // namespace LAMMPS_NS

#endif
#endif
