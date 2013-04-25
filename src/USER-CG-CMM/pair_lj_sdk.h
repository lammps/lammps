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

PairStyle(lj/sdk,PairLJSDK)
PairStyle(cg/cmm,PairLJSDK)

#else

#ifndef LMP_PAIR_LJ_SDK_H
#define LMP_PAIR_LJ_SDK_H

#include "pair.h"

namespace LAMMPS_NS {
class LAMMPS;

class PairLJSDK : public Pair {
 public:
  PairLJSDK(LAMMPS *);
  virtual ~PairLJSDK();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);
  double single(int, int, int, int, double, double, double, double &);
  void *extract(const char *, int &);
  virtual double memory_usage();

 protected:
  int **lj_type; // type of lennard jones potential

  double **cut;
  double **epsilon,**sigma;
  double **lj1,**lj2,**lj3,**lj4,**offset;

  // cutoff and offset for minimum of LJ potential
  // to be used in SDK angle potential, which
  // uses only the repulsive part of the potential

  double **rminsq, **emin;

  double cut_global;

  void allocate();

 private:
  template <int EVFLAG, int EFLAG, int NEWTON_PAIR> void eval();

};

}

#endif
#endif
