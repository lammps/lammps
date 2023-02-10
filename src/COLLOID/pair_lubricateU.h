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
PairStyle(lubricateU,PairLubricateU);
// clang-format on
#else

#ifndef LMP_PAIR_LUBRICATEU_H
#define LMP_PAIR_LUBRICATEU_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLubricateU : public Pair {
 public:
  PairLubricateU(class LAMMPS *);
  ~PairLubricateU() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  void init_style() override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

 protected:
  double cut_inner_global, cut_global;
  double mu;
  double rad;
  int flaglog;
  int flagdeform, flagwall;
  int flagVF, flagHI;
  double vol_P;
  class FixWall *wallfix;

  double gdot, Ef[3][3];
  double **cut_inner, **cut;
  void allocate();
  double R0, RT0, RS0;

  int nmax;
  double **fl, **Tl, **xl;

  int cgmax;
  double *bcg, *xcg, *rcg, *rcg1, *pcg, *RU;

  void compute_RE();
  virtual void compute_RE(double **);
  void compute_RU();
  virtual void compute_RU(double **);
  virtual void compute_Fh(double **);
  void stage_one();
  void intermediates(int, double **);
  void stage_two(double **);
  void copy_vec_uo(int, double *, double **, double **);
  void copy_uo_vec(int, double **, double **, double *);
  double dot_vec_vec(int, double *, double *);
};

}    // namespace LAMMPS_NS

#endif
#endif
