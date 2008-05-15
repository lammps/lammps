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

#ifndef PAIR_LJ_GROMACS_COUL_GROMACS_H
#define PAIR_LJ_GROMACS_COUL_GROMACS_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLJGromacsCoulGromacs : public Pair {
 public:
  PairLJGromacsCoulGromacs(class LAMMPS *);
  virtual ~PairLJGromacsCoulGromacs();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  virtual void write_restart_settings(FILE *);
  virtual void read_restart_settings(FILE *);
  virtual double single(int, int, int, int, double, double, double, double &);

 protected:
  double cut_lj_inner,cut_lj,cut_coul_inner,cut_coul;
  double cut_lj_innersq,cut_ljsq,cut_coul_innersq,cut_coulsq,cut_bothsq;
  double **epsilon,**sigma;
  double **lj1,**lj2,**lj3,**lj4;
  double **ljsw1,**ljsw2,**ljsw3,**ljsw4,**ljsw5;
  double coulsw1,coulsw2,coulsw3,coulsw4,coulsw5;

  void allocate();
};

}

#endif
