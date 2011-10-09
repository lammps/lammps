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

namespace LAMMPS_NS {

// per thread data accumulators
class ThrData {
  friend class FixOMP;
  friend class ThrOMP;

 public:
  ThrData(int tid);
  ~ThrData();

  void clear(int);             // erase accumulator contents
  void grow_arrays(int);       // grow per atom arrays
  void set_accflags(int flags) { _accflags = flags; }; // flag which accumulators to prepare
  void signal_reduce(int flag) { _redflags |= flag; }; // signal which reductions are needed

 protected:
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

  double *eatom_pair;     // per atom total energy from non-bonded
  double *vatom_pair;     // per atom virial from non-bonded
  double *eatom_bond;     // per atom total energy from bonds
  double *vatom_bond;     // per atom virial from bonds
  double *eatom_angle;    // per atom total energy from angles
  double *vatom_angle;    // per atom virial from angles
  double *eatom_dihed;    // per atom total energy from dihedrals
  double *vatom_dihed;    // per atom virial from dihedrals
  double *eatom_imprp;    // per atom total energy from impropers
  double *vatom_imprp;    // per atom virial from impropers
  double *eatom_kspce;    // per atom total energy from kspace
  double *vatom_kspce;    // per atom virial from kspace

 private:
  int _maxeatom;          // size of eatom array
  int _maxvatom;          // size of vatom array
  int _tid;               // my thread id
  int _accflags;          // flags indicating which accumulators to provide
  int _redflags;          // flags indicating which property to reduce

 public:
  enum {THR_NONE=0,THR_ENERGY=1<<0,THR_VIRIAL=1<<1,THR_EATOM=1<<2,THR_VATOM=1<<3,
	THR_PAIR=1<<4,THR_BOND=1<<5,THR_ANGLE=1<<6,THR_DIHEDRAL=1<<7,
	THR_IMPROPER=1<<8,THR_KSPACE=1<<9,THR_VFDOTR=1<<10};
    
 public:
  double memory_usage();

 // disabled default methods
 private:
  ThrData() {};
};

// reduce per thread data into the first part of the data
// array that is used for the non-threaded parts and reset
// the temporary storage to 0.0. this routine depends on
// multi-dimensional arrays like force stored in this order
// x1,y1,z1,x2,y2,z2,...
// we need to post a barrier to wait until all threads are done
// with writing to the array .
static void data_reduce_thr(double *dall, int nall, int nthreads, int ndim, int tid)
{
#if defined(_OPENMP)
  // NOOP in non-threaded execution.
  if (nthreads == 1) return;
#pragma omp barrier
  {
    const int nvals = ndim*nall;
    const int idelta = nvals/nthreads + 1;
    const int ifrom = tid*idelta;
    const int ito   = ((ifrom + idelta) > nvals) ? nvals : (ifrom + idelta);

    for (int m = ifrom; m < ito; ++m) {
      for (int n = 1; n < nthreads; ++n) {
	dall[m] += dall[n*nvals + m];
	dall[n*nvals + m] = 0.0;
      }
    }
  }
#else
  // NOOP in non-threaded execution.
  return;
#endif
}

}
#endif
