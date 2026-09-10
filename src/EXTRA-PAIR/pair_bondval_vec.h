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
PairStyle(bondval/vec,PairBondValVec)
// clang-format on
#else

#ifndef LMP_PAIR_BONDVAL_VEC_H
#define LMP_PAIR_BONDVAL_VEC_H

#include "pair.h"

namespace LAMMPS_NS {

class PairBondValVec : public Pair {
 public:
  PairBondValVec(class LAMMPS *);
  ~PairBondValVec() override;

  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;

  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

 protected:
  double cut_global;
  double power_global;
  double **cut;
  double **r0, **alpha;
  double *fp;
  double **bvvsparam, **bvvv0;
  double *s0x, *s0y, *s0z;
  double *Dix, *Diy, *Diz;
  double *s0temp;
  double *fpx, *fpy, *fpz;
  int index, index2;
  int nmax;
  double **offset;
  double *cut_respa;
  virtual void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
