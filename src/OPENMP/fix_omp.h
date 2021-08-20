/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(OMP,FixOMP);
// clang-format on
#else

#ifndef LMP_FIX_OMP_H
#define LMP_FIX_OMP_H

#include "fix.h"

namespace LAMMPS_NS {

class ThrData;

class FixOMP : public Fix {
  friend class ThrOMP;
  friend class RespaOMP;

 public:
  FixOMP(class LAMMPS *, int, char **);
  virtual ~FixOMP();
  virtual int setmask();
  virtual void init();
  virtual void setup(int);
  virtual void min_setup(int flag) { setup(flag); }
  virtual void pre_force(int);

  virtual void setup_pre_force(int vflag) { pre_force(vflag); }
  virtual void min_setup_pre_force(int vflag) { pre_force(vflag); }
  virtual void min_pre_force(int vflag) { pre_force(vflag); }
  virtual void setup_pre_force_respa(int vflag, int) { pre_force(vflag); }
  virtual void pre_force_respa(int vflag, int, int) { pre_force(vflag); }

  virtual double memory_usage();

 protected:
  ThrData **thr;
  void *last_omp_style;      // pointer to the style that needs
                             // to do the general force reduction
  void *last_pair_hybrid;    // pointer to the pair style that needs
                             // to call virial_fdot_compute()
  // signal that an /omp style did the force reduction. needed by respa/omp
  void did_reduce() { _reduced = true; }

 public:
  ThrData *get_thr(int tid) { return thr[tid]; }
  int get_nthr() const { return _nthr; }

  bool get_neighbor() const { return _neighbor; }
  bool get_mixed() const { return _mixed; }
  bool get_reduced() const { return _reduced; }

 private:
  int _nthr;                    // number of currently active ThrData objects
  bool _neighbor;               // en/disable threads for neighbor list construction
  bool _mixed;                  // whether to prefer mixed precision compute kernels
  bool _reduced;                // whether forces have been reduced for this step
  bool _pair_compute_flag;      // whether pair_compute is called
  bool _kspace_compute_flag;    // whether kspace_compute is called

  void set_neighbor_omp();
};

}    // namespace LAMMPS_NS

#endif
#endif
