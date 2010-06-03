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

PairStyle(cg/cmm/omp,PairCGCMMOMP)

#else

#ifndef LMP_PAIR_CG_CMM_OMP_H
#define LMP_PAIR_CG_CMM_OMP_H

#include "pair_omp.h"
#include "cg_cmm_parms.h"

namespace LAMMPS_NS {

class PairCGCMMOMP : public PairOMP, public CGCMMParms {
    
 public:

  PairCGCMMOMP(class LAMMPS *);
  virtual ~PairCGCMMOMP();

  virtual void compute(int, int);
  virtual double single(int, int, int, int, double, double, double, double &);

  virtual void settings(int, char **);
  virtual void coeff(int, char **);
  virtual double init_one(int, int);

  virtual void write_restart(FILE *);
  virtual void read_restart(FILE *);
  virtual void write_restart_settings(FILE *);
  virtual void read_restart_settings(FILE *);

  virtual double memory_usage();

 protected:
  // coarse grain flags
  int **cg_type;
    
  // lennard jones parameters
  double cut_lj_global, **cut, **cut_lj, **cut_ljsq;
  double **epsilon, **sigma;
  double **lj1, **lj2, **lj3, **lj4, **offset;

  // coulomb parameters
  int allocated_coul; // 0/1 = whether coulomb arrays are allocated
  double cut_coul_global, cut_coulsq_global, kappa, g_ewald;
  double **cut_coul, **cut_coulsq;

  // tables
  double tabinnersq;
  double *rtable,*drtable,*ftable,*dftable,*ctable,*dctable;
  double *etable,*detable,*ptable,*dptable,*vtable,*dvtable;
  int ncoulshiftbits,ncoulmask;

  // r-RESPA parameters
  double *cut_respa;

protected:
  virtual void allocate();

  template < int EVFLAG, int EFLAG, int NEWTON_PAIR > void eval();
 
private:
  // disable default constructor
  PairCGCMMOMP();
};
}

#endif
#endif
