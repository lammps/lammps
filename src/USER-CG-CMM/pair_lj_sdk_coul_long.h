/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(lj/sdk/coul/long,PairLJSDKCoulLong)
PairStyle(cg/cmm/coul/long,PairLJSDKCoulLong)

#else

#ifndef LMP_PAIR_LJ_SDK_COUL_LONG_H
#define LMP_PAIR_LJ_SDK_COUL_LONG_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLJSDKCoulLong : public Pair {
 public:
  PairLJSDKCoulLong(class LAMMPS *);
  virtual ~PairLJSDKCoulLong();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);
  virtual void write_restart_settings(FILE *);
  virtual void read_restart_settings(FILE *);
  virtual double single(int, int, int, int, double, double, double, double &);
  virtual void *extract(const char *, int &);
  virtual double memory_usage();

 protected:
  double **cut_lj,**cut_ljsq;
  double cut_coul,cut_coulsq;
  double **epsilon,**sigma;
  double **lj1,**lj2,**lj3,**lj4,**offset;
  int **lj_type;

  // cutoff and offset for minimum of LJ potential
  // to be used in SDK angle potential, which
  // uses only the repulsive part of the potential

  double **rminsq, **emin;

  double cut_lj_global;
  double g_ewald;

  double tabinnersq;
  double *rtable,*drtable,*ftable,*dftable,*ctable,*dctable;
  double *etable,*detable;
  int ncoulshiftbits,ncoulmask;

  void allocate();
  void init_tables();
  void free_tables();

 private:
  template <int EVFLAG, int EFLAG, int NEWTON_PAIR> void eval();

};

}

#endif
#endif
