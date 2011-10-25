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

#if 0

#ifdef PAIR_CLASS

PairStyle(lj/sdk,PairLJSDK)
PairStyle(cg/cmm,PairLJSDK)

#else

#ifndef LMP_PAIR_LJ_SDK_H
#define LMP_PAIR_LJ_SDK_H

#include "lj_sdk_common.h"

namespace LAMMPS_NS {

  class PairLJSDK : PairLJSDKCommon {

    public:

    PairLJSDK(class LAMMPS *);
    virtual ~PairLJSDK();

    virtual void compute(int, int);

    void write_restart(FILE *);
    void read_restart(FILE *);

    virtual double single(int, int, int, int, double, double, double, double &);

    private:
    // disable default constructor
    PairLJSDK();
  };
}

#endif
#endif


#endif
