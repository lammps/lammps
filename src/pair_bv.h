/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(bv,PairBV)

#else

#ifndef LMP_PAIR_BV_H
#define LMP_PAIR_BV_H

#include "pair.h"

namespace LAMMPS_NS {

class PairBV : public Pair {
 public:
  PairBV(class LAMMPS *);
  virtual ~PairBV() override;

  virtual void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;

  virtual int pack_forward_comm(int, int *, double *, int, int *) override;
  virtual void unpack_forward_comm(int, int, double *) override;

  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;


 protected:
  double cut_global;
  int nmax;  
  double power_global;
  double **cut;
  double **r0,**alpha,**sparam,**v0;
  double *s0,*fp,*energy0;
  double **offset;
  double *cut_respa;

  virtual void allocate();
};

}

#endif
#endif

