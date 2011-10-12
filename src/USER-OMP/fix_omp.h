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

#ifdef FIX_CLASS

FixStyle(OMP,FixOMP)

#else

#ifndef LMP_FIX_OMP_H
#define LMP_FIX_OMP_H

#include "fix.h"


namespace LAMMPS_NS {

class ThrData;

class FixOMP : public Fix {
  friend class ThrOMP;

 public:
  FixOMP(class LAMMPS *, int, char **);
  virtual ~FixOMP();
  virtual int setmask();
  virtual void grow_arrays(int);

  virtual void setup_pre_force(int);
  virtual void setup(int);
  virtual void pre_force(int);
  virtual void post_force(int);

  virtual double memory_usage();

  ThrData *get_thr(int tid)  { return thr[tid]; };

 protected:
  ThrData **thr;
  int _ndata;    // number of per thread data sets
  int _nall;     // cached number of local + ghost atoms
  int _nmax;     // cached number of max atoms per MPI task
  int _redflags; // flags for reduction
  
  void add_reduce(int flags) {_redflags |= flags; };

#if 0
  virtual void initial_integrate(int);
  virtual void initial_integrate(int,int,int);
  virtual void setup_pre_force(int);
  virtual void setup_post_force(int);
  virtual void pre_force_respa(int);
  virtual void post_force_respa(int);
  virtual void min_pre_force(int);
  virtual void min_post_force(int);


 public:
  bool get_newton() const {return _newton;};

 private:
  bool _newton;  // en/disable newton's 3rd law for local atoms.
#endif

};

}

#endif
#endif
