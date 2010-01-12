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

PairStyle(peri/pmb,PairPeriPMB)

#else

#ifndef LMP_PAIR_PERI_PMB_H
#define LMP_PAIR_PERI_PMB_H

#include "pair.h"

namespace LAMMPS_NS {

class PairPeriPMB : public Pair {
 public:
  PairPeriPMB(class LAMMPS *);
  ~PairPeriPMB();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void init_style();
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *) {}
  void read_restart_settings(FILE *) {}
  double single(int, int, int, int, double, double, double, double &);
  double memory_usage();

 protected:
  int ifix_peri;
  double **kspring;
  double **s00, **alpha;
  double **cut;
 
  double *s0_new;
  int nmax;

  void allocate();
};

}

#endif
#endif
