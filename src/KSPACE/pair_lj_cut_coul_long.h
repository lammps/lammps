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

#ifndef PAIR_LJ_CUT_COUL_LONG_H
#define PAIR_LJ_CUT_COUL_LONG_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLJCutCoulLong : public Pair {
 public:
  PairLJCutCoulLong(class LAMMPS *);
  virtual ~PairLJCutCoulLong();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  void coeff(int, char **);
  virtual void init_style();
  void init_list(int, class NeighList *);
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  virtual void write_restart_settings(FILE *);
  virtual void read_restart_settings(FILE *);
  virtual double single(int, int, int, int, double, double, double, double &);

  void compute_inner();
  void compute_middle();
  void compute_outer(int, int);
  void *extract(char *);

 protected:
  double cut_lj_global;
  double **cut_lj,**cut_ljsq;
  double cut_coul,cut_coulsq;
  double **epsilon,**sigma;
  double **lj1,**lj2,**lj3,**lj4,**offset;
  double *cut_respa;
  double g_ewald;

  double tabinnersq;
  double *rtable,*drtable,*ftable,*dftable,*ctable,*dctable;
  double *etable,*detable,*ptable,*dptable,*vtable,*dvtable;
  int ncoulshiftbits,ncoulmask;

  void allocate();
  void init_tables();
  void free_tables();
};

}

#endif
