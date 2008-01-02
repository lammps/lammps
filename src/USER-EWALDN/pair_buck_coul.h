/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef PAIR_BUCK_COUL_H
#define PAIR_BUCK_COUL_H

#include "pair.h"

namespace LAMMPS_NS {

class PairBuckCoul : public Pair {
 public:
  double cut_coul;

  PairBuckCoul(class LAMMPS *);
  ~PairBuckCoul();
  virtual void compute(int, int);

  virtual void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  void init_list(int, class NeighList *);
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  double single(int, int, int, int, double, double, double, double &);
  void *extract(char *);

  void compute_inner();
  void compute_middle();
  void compute_outer(int, int);

 protected:
  double cut_buck_global;
  double **cut_buck, **cut_buck_read, **cut_bucksq;
  double cut_coulsq;
  double **buck_a_read, **buck_a, **buck_c_read, **buck_c;
  double **buck1, **buck2, **buck_rho_read, **buck_rho, **rhoinv, **offset;
  double *cut_respa;
  double g_ewald;
  int ewald_order, ewald_off;

  double tabinnersq;
  double *rtable, *drtable, *ftable, *dftable, *ctable, *dctable;
  double *etable, *detable, *ptable, *dptable, *vtable, *dvtable;
  int ncoulshiftbits, ncoulmask;

  void options(char **arg, int order);
  void allocate();
  void init_tables();
  void free_tables();
};

}

#endif
