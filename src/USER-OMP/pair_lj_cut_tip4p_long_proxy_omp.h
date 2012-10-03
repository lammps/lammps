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

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(lj/cut/tip4p/long/proxy/omp,PairLJCutTIP4PLongProxyOMP)

#else

#ifndef LMP_PAIR_LJ_CUT_TIP4P_LONG_PROXY_OMP_H
#define LMP_PAIR_LJ_CUT_TIP4P_LONG_PROXY_OMP_H

#include "pair_lj_cut_tip4p_long.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class PairLJCutTIP4PLongProxyOMP : 
    public PairLJCutTIP4PLong, public ThrOMP {

 public:
  PairLJCutTIP4PLongProxyOMP(class LAMMPS *);
  virtual ~PairLJCutTIP4PLongProxyOMP() {};

  virtual void init_style();

  virtual void compute(int, int);
  virtual double memory_usage();

 private:
  template <int, int, int, int>
  void eval(int ifrom, int ito, ThrData * const thr);
  void compute_newsite_thr(const double *, const double *,
                           const double *, double *) const;

  class PPPMTIP4PProxy *kspace;
  int nproxy;
};

}

#endif
#endif
