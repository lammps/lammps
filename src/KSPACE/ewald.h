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

#ifndef EWALD_H
#define EWALD_H

#include "kspace.h"

namespace LAMMPS_NS {

class Ewald : public KSpace {
 public:
  Ewald(class LAMMPS *, int, char **);
  ~Ewald();
  void init();
  void setup();
  void compute(int, int);
  double memory_usage();

 private:
  double PI;
  double precision;
  int kcount,kmax,kmax3d,kmax_created;
  double qqrd2e;
  double gsqmx,qsum,qsqsum,volume;
  int nmax;

  double unitk[3];
  int *kxvecs,*kyvecs,*kzvecs;
  double *ug;
  double **eg,**vg;
  double **ek;
  double *sfacrl,*sfacim,*sfacrl_all,*sfacim_all;
  double ***cs,***sn;

  void eik_dot_r();
  void coeffs();
  void allocate();
  void deallocate();
  void slabcorr(int);
};

}

#endif


