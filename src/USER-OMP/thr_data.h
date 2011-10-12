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

namespace LAMMPS_NS {

enum {THR_NONE=0,THR_PAIR=1,THR_BOND=1<<1,THR_ANGLE=1<<2,
      THR_DIHEDRAL=1<<3,THR_IMPROPER=1<<4,THR_KSPACE=1<<5};
    
// per thread data accumulators
class ThrData {
  friend class FixOMP;
  friend class ThrOMP;

 public:
  ThrData(int tid) : _tid(tid) {};
  ~ThrData() {};

  void clear();           // erase accumulator contents
  void check_tid(int);    // thread id consistency check

 public:
  double eng_vdwl;        // non-bonded non-coulomb energy
  double eng_coul;        // non-bonded coulomb energy
  double virial_pair[6];  // virial contribution from non-bonded
  double eng_bond;        // bond energy
  double virial_bond[6];  // virial contribution from bonds
  double eng_angle;       // angle energy
  double virial_angle[6]; // virial contribution from angles
  double eng_dihed;       // dihedral energy
  double virial_dihed[6]; // virial contribution from dihedrals
  double eng_imprp;       // improper energy
  double virial_imprp[6]; // virial contribution from impropers
  double eng_kspce;       // kspace energy
  double virial_kspce[6]; // virial contribution from kspace

  double *eatom;          // per atom energy array segment for this thread
  double **vatom;         // per atom virial array segment for this thread

 private:
  int _tid;               // my thread id

 public:
  double memory_usage();

 // disabled default methods
 private:
  ThrData() {};
};

////////////////////////////////////////////////////////////////////////
//  helper functions operating on data replicated for thread support  //
////////////////////////////////////////////////////////////////////////
// compute global per thread virial contribution from per-thread force
void virial_fdotr_compute_thr(double * const, double *, double *, int, int, int);
// generic per thread data reduction for continous arrays of nthreads*nmax size
void data_reduce_thr(double *, int, int, int, int);
/* ---------------------------------------------------------------------- */

// set loop range thread id, and force array offset for threaded runs.
static double **loop_setup_thr(double **f, int &ifrom, int &ito, int &tid,
			       int inum, int nall, int nthreads)
{
#if defined(_OPENMP)
  tid = omp_get_thread_num();

  // each thread works on a fixed chunk of atoms.
  const int idelta = 1 + inum/nthreads;
  ifrom = tid*idelta;
  ito   = ((ifrom + idelta) > inum) ? inum : ifrom + idelta;

  return f + nall*tid;
#else
  tid = 0;
  ifrom = 0;
  ito = inum;
  return f;
#endif
}

}
#endif
