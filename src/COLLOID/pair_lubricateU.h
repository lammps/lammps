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
  virtual ~PairLubricateU();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  virtual void init_style();
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);

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

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

W: Cannot include log terms without 1/r terms; setting flagHI to 1.

Self-explanatory.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair lubricateU requires atom style sphere

Self-explanatory.

E: Pair lubricateU requires ghost atoms store velocity

Use the comm_modify vel yes command to enable this.

E: Pair lubricateU requires monodisperse particles

All particles must be the same finite size.

E: Cannot use multiple fix wall commands with pair lubricateU

Self-explanatory.

*/
