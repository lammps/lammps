/* -*- c++ -*- -------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#ifndef LMP_THR_OMP_H
#define LMP_THR_OMP_H

#include "pointers.h"

namespace LAMMPS_NS {

class ThrOMP : protected Pointers {

 protected:
  const int thr_style;
  enum {PAIR=1, BOND, ANGLE, DIHEDRAL, IMPROPER, KSPACE, FIX, COMPUTE};

  double *eng_vdwl_thr;  // per thread accumulated vdw energy
  double *eng_coul_thr;  // per thread accumulated coulomb energies
  double *eng_bond_thr;  // per thread accumlated bonded energy

  double **virial_thr;   // per thread virial
  double **eatom_thr;    // per thread per atom energy
  double ***vatom_thr;   // per thread per atom virial

  int maxeatom_thr, maxvatom_thr;
  
 public:
  ThrOMP(class LAMMPS *, int);
  virtual ~ThrOMP();

 protected:
  // set loop range for, thread id, and force array offset for threaded runs.
  double **loop_setup_thr(double **f, int &ifrom, int &ito, int &tid,
                          const int inum, const int nall, const int nthreads);
  
#if 0
  // threading adapted versions of the ev_tally infrastructure.
  void ev_setup_thr(int, int);
  void ev_reduce_thr();
  void ev_tally_thr(int, int, int, int, double, double, double,
		    double, double, double, int);
  void ev_tally_xyz_thr(int, int, int, int, double, double,
			double, double, double, double, double, double, int);
  void ev_tally3_thr(int, int, int, double, double,
		     double *, double *, double *, double *, int);
  void ev_tally4_thr(int, int, int, int, double,
		     double *, double *, double *, double *, double *, double *, int);
  void ev_tally_list_thr(int, int *, double, double *, int);
  void v_tally2_thr(int, int, double, double *, int);
  void v_tally3_thr(int, int, int, double *, double *, double *, double *, int);
  void v_tally4_thr(int, int, int, int, double *, double *, double *,
                    double *, double *, double *, int);
#endif

  // reduce per thread forces into the first part of the force array
  void force_reduce_thr(double **fall, const int nall,
			const int nthreads, const int tid);
#if 0
  // reduce per thread density into the first part of the rho array
  void rho_reduce_thr(double *rho, const int nmax, const int nrange, 
		      const int nthreads, const int tid);
#endif
};

}
#endif
