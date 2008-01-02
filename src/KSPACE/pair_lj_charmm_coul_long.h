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

#ifndef PAIR_LJ_CHARMM_COUL_LONG_H
#define PAIR_LJ_CHARMM_COUL_LONG_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLJCharmmCoulLong : public Pair {
 public:
  PairLJCharmmCoulLong(class LAMMPS *);
  ~PairLJCharmmCoulLong();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  void init_list(int, class NeighList *);
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  double single(int, int, int, int, double, double, double, double &);

  void compute_inner();
  void compute_middle();
  void compute_outer(int, int);
  void *extract(char *);

 protected:
  int implicit;
  double cut_lj_inner,cut_lj;
  double cut_lj_innersq,cut_ljsq;
  double cut_coul,cut_coulsq;
  double cut_bothsq;
  double denom_lj;
  double **epsilon,**sigma,**eps14,**sigma14;
  double **lj1,**lj2,**lj3,**lj4,**offset;
  double **lj14_1,**lj14_2,**lj14_3,**lj14_4;
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
