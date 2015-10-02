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

#ifndef LMP_THR_DATA_H
#define LMP_THR_DATA_H

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "timer.h"

namespace LAMMPS_NS {

// per thread data accumulators
// there should be one instance
// of this class for each thread.
class ThrData {
  friend class FixOMP;
  friend class ThrOMP;

 public:
  ThrData(int tid, class Timer *t);
  ~ThrData() { delete _timer; _timer = NULL; };

  void check_tid(int);    // thread id consistency check
  int get_tid() const { return _tid; }; // our thread id.

  // inline wrapper, to make this more efficient
  // when per-thread timers are off
  void timer(enum Timer::ttype flag) { if (_timer) _stamp(flag); };
  double get_time(enum Timer::ttype flag);

  // erase accumulator contents and hook up force arrays
  void init_force(int, double **, double **, double *, double *, double *);

  // give access to per-thread offset arrays
  double **get_f() const { return _f; };
  double **get_torque() const { return _torque; };
  double *get_de() const { return _de; };
  double *get_drho() const { return _drho; };

  // setup and erase per atom arrays
  void init_adp(int, double *, double **, double **); // ADP (+ EAM)
  void init_cdeam(int, double *, double *, double *); // CDEAM (+ EAM)
  void init_eam(int, double *);                       // EAM
  void init_eim(int, double *, double *);             // EIM (+ EAM)

  void init_pppm(int, class Memory *);
  void init_pppm_disp(int, class Memory *);

  // access methods for arrays that we handle in this class
  double **get_lambda() const { return _lambda; };
  double **get_mu() const { return _mu; };
  double *get_D_values() const { return _D_values; };
  double *get_fp() const { return _fp; };
  double *get_rho() const { return _rho; };
  double *get_rhoB() const { return _rhoB; };
  void *get_rho1d() const { return _rho1d; };
  void *get_drho1d() const { return _drho1d; };
  void *get_rho1d_6() const { return _rho1d_6; };
  void *get_drho1d_6() const { return _drho1d_6; };

 private:
  double eng_vdwl;        // non-bonded non-coulomb energy
  double eng_coul;        // non-bonded coulomb energy
  double eng_bond;        // bond energy
  double eng_angle;       // angle energy
  double eng_dihed;       // dihedral energy
  double eng_imprp;       // improper energy
  double eng_kspce;       // kspace energy
  double virial_pair[6];  // virial contribution from non-bonded
  double virial_bond[6];  // virial contribution from bonds
  double virial_angle[6]; // virial contribution from angles
  double virial_dihed[6]; // virial contribution from dihedrals
  double virial_imprp[6]; // virial contribution from impropers
  double virial_kspce[6]; // virial contribution from kspace
  double *eatom_pair;
  double *eatom_bond;
  double *eatom_angle;
  double *eatom_dihed;
  double *eatom_imprp;
  double *eatom_kspce;
  double **vatom_pair;
  double **vatom_bond;
  double **vatom_angle;
  double **vatom_dihed;
  double **vatom_imprp;
  double **vatom_kspce;

  // per thread segments of various force or similar arrays

  // these are maintained by atom styles
  double **_f;
  double **_torque;
  double *_erforce;
  double *_de;
  double *_drho;

  // these are maintained by individual pair styles
  double **_mu, **_lambda;   // ADP (+ EAM)
  double *_rhoB, *_D_values; // CDEAM (+ EAM)
  double *_rho;              // EAM
  double *_fp;               // EIM (+ EAM)

  // this is for pppm/omp
  void *_rho1d;
  void *_drho1d;
  // this is for pppm/disp/omp
  void *_rho1d_6;
  void *_drho1d_6;
  // my thread id
  const int _tid;
  // timer info
  int _timer_active;
  class Timer *_timer;

 private:
  void _stamp(enum Timer::ttype flag);

 public:
  // compute global per thread virial contribution from global forces and positions
  void virial_fdotr_compute(double **, int, int, int);

  double memory_usage();

 // disabled default methods
 private:
  ThrData() : _tid(-1), _timer(NULL) {};
};

////////////////////////////////////////////////////////////////////////
//  helper functions operating on data replicated for thread support  //
////////////////////////////////////////////////////////////////////////
// generic per thread data reduction for continous arrays of nthreads*nmax size
void data_reduce_thr(double *, int, int, int, int);
}
#endif
