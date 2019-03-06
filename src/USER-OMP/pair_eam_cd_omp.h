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

PairStyle(eam/cd/omp,PairEAMCD_OneSiteOMP)
PairStyle(eam/cd/old/omp,PairEAMCD_TwoSiteOMP)

#else

#ifndef LMP_PAIR_EAM_CD_OMP_H
#define LMP_PAIR_EAM_CD_OMP_H

#include "pair_eam_cd.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class PairEAMCDOMP : public PairEAMCD, public ThrOMP {

 public:
  PairEAMCDOMP(class LAMMPS *, int);

  virtual void compute(int, int);
  virtual double memory_usage();

 private:
  template <int EVFLAG, int EFLAG, int NEWTON_PAIR, int CDEAMVERSION>
  void eval(int iifrom, int iito, ThrData * const thr);
};

  /// The one-site concentration formulation of CD-EAM.
  class PairEAMCD_OneSiteOMP : public PairEAMCDOMP
  {
  public:
    /// Constructor.
    PairEAMCD_OneSiteOMP(class LAMMPS* lmp) : PairEAM(lmp), PairEAMCDOMP(lmp, 1) {}
  };

  /// The two-site concentration formulation of CD-EAM.
  class PairEAMCD_TwoSiteOMP : public PairEAMCDOMP
  {
  public:
    /// Constructor.
    PairEAMCD_TwoSiteOMP(class LAMMPS* lmp) : PairEAM(lmp), PairEAMCDOMP(lmp, 2) {}
  };

}

#endif
#endif
