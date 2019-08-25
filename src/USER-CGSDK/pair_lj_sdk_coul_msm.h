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

PairStyle(lj/sdk/coul/msm,PairLJSDKCoulMSM)

#else

#ifndef LMP_PAIR_LJ_SDK_COUL_MSM_H
#define LMP_PAIR_LJ_SDK_COUL_MSM_H

#include "pair_lj_sdk_coul_long.h"

namespace LAMMPS_NS {

class PairLJSDKCoulMSM : public PairLJSDKCoulLong {
 public:
  PairLJSDKCoulMSM(class LAMMPS *);
  virtual ~PairLJSDKCoulMSM() {};
  virtual void compute(int, int);
  virtual double single(int, int, int, int, double, double, double, double &);
  virtual void *extract(const char *, int &);

private:
  template <int EVFLAG, int EFLAG, int NEWTON_PAIR> void eval_msm();

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Must use 'kspace_modify pressure/scalar no' with Pair style

The kspace scalar pressure option is not (yet) compatible with at least one of
the defined Pair styles.

*/
