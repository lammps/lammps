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

#ifndef PAIR_COUL_LONG_H
#define PAIR_COUL_LONG_H

#include "pair.h"

namespace LAMMPS_NS {

class PairCoulLong : public Pair {
 public:
  PairCoulLong(class LAMMPS *);
  ~PairCoulLong();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  double single(int, int, int, int, double, double, double, double &);
  void *extract(char *);

 private:
  double cut_coul,cut_coulsq;
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
