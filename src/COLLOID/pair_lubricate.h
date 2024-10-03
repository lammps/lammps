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
PairStyle(lubricate,PairLubricate);
// clang-format on
#else

#ifndef LMP_PAIR_LUBRICATE_H
#define LMP_PAIR_LUBRICATE_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLubricate : public Pair {
 public:
  PairLubricate(class LAMMPS *);
  ~PairLubricate() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  void init_style() override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  int pre_adapt(char *, int, int, int, int);
  void adapt(int, int, int, int, int, double);

  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

 protected:
  double mu, cut_inner_global, cut_global;
  double rad;
  int flaglog, flagfld, shearing;
  int flagdeform, flagwall;
  double vol_P;
  class FixWall *wallfix;
  int flagVF, flagHI;

  double Ef[3][3];
  double R0, RT0, RS0;
  double **cut_inner, **cut;

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
