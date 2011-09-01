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

#ifdef FIX_CLASS

FixStyle(gcmc,FixGCMC)

#else

#ifndef LMP_FIX_GCMC_H
#define LMP_FIX_GCMC_H

#include "stdio.h"
#include "fix.h"

namespace LAMMPS_NS {

class FixGCMC : public Fix {
 public:
  FixGCMC(class LAMMPS *, int, char **);
  ~FixGCMC();
  int setmask();
  void init();
  void pre_exchange();
  void attempt_move();
  void attempt_deletion();
  void attempt_insertion();
  double energy(int, double *);
  double compute_vector(int);
  double memory_usage();
  void write_restart(FILE *);
  void restart(char *);

 private:
  int ntype,nevery,seed;
  int ncycles,nexchanges,nmcmoves;
  int ngas;           // # of gas molecules (or atoms) on all procs 
  int ngas_local;     // # of gas molecules (or atoms) on this proc 
  int ngas_before;    // # of gas molecules (or atoms) on procs < this proc
  int molflag;        // 0 = atomic, 1 = molecular system                                                
  double nmove_attempts;   
  double nmove_successes;  
  double ndel_attempts;    
  double ndel_successes;   
  double ninsert_attempts; 
  double ninsert_successes;
  
  int nmax;
  double reservoir_temperature;
  double chemical_potential;
  double displace;
  double beta,zz,sigma,volume;
  double xlo,xhi,ylo,yhi,zlo,zhi;
  double *sublo,*subhi;
  int *local_gas_list;                           
  double **cutsq;
  class Pair *pair;
 
  class RanPark *random_equal;
  class RanPark *random_unequal;

  void options(int, char **);
};

}

#endif
#endif
