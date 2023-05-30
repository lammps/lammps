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
   Contributing authors: Eugen Rozic (University College London)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(cosine/squared, PairCosineSquared);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_COS_SQ_H
#define LMP_PAIR_LJ_COS_SQ_H

#include "pair.h"

namespace LAMMPS_NS {

class PairCosineSquared : public Pair {
 public:
  PairCosineSquared(class LAMMPS *);
  ~PairCosineSquared() override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  // void init_style();
  double init_one(int, int) override;
  void modify_params(int, char **) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  void write_data(FILE *) override;
  void write_data_all(FILE *) override;
  void compute(int, int) override;
  double single(int, int, int, int, double, double, double, double &) override;
  // void *extract(const char *, int &);

  /* RESPA stuff not implemented...
  void compute_inner();
  void compute_middle();
  void compute_outer(int, int);
*/

 protected:
  double cut_global;
  double **epsilon, **sigma, **w, **cut;
  int **wcaflag;
  double **lj12_e, **lj6_e, **lj12_f, **lj6_f;

  virtual void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
