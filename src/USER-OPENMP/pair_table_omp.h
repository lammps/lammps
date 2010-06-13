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

#ifdef PAIR_CLASS

PairStyle(table/omp,PairTableOMP)

#else

#ifndef LMP_PAIR_TABLE_OMP_H
#define LMP_PAIR_TABLE_OMP_H

#include "pair_omp.h"

namespace LAMMPS_NS {

class PairTableOMP : public PairOMP {
 public:
  PairTableOMP(class LAMMPS *);
  ~PairTableOMP();

  virtual void compute(int, int);
  virtual double single(int, int, int, int, double, double, double, double &);

  virtual void settings(int, char **);
  virtual void coeff(int, char **);

  virtual double init_one(int, int);

  virtual void write_restart(FILE *);
  virtual void read_restart(FILE *);
  virtual void write_restart_settings(FILE *);
  virtual void read_restart_settings(FILE *);
  virtual void *extract(char *);

 protected:
  template <int TABSTYLE, int EVFLAG, int EFLAG, int NEWTON_PAIR> void eval();

 protected:
  int tabstyle,tablength;
  struct Table {
    int ninput,rflag,fpflag,match,ntablebits;
    int nshiftbits,nmask;
    double rlo,rhi,fplo,fphi,cut;
    double *rfile,*efile,*ffile;
    double *e2file,*f2file;
    double innersq,delta,invdelta,deltasq6;
    double *rsq,*drsq,*e,*de,*f,*df,*e2,*f2;
  };
  int ntables;
  Table *tables;

  int **tabindex;

  void allocate();
  void read_table(Table *, char *, char *);
  void param_extract(Table *, char *);
  void bcast_table(Table *);
  void spline_table(Table *);
  void compute_table(Table *);
  void null_table(Table *);
  void free_table(Table *);
  void spline(double *, double *, int, double, double, double *);
  double splint(double *, double *, double *, int, double);
};

}

#endif
#endif
