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

#ifdef PAIR_CLASS

PairStyle(cg/cmm/coul/msm,PairCGCMMCoulMSM)

#else

#ifndef LMP_PAIR_CG_CMM_COUL_MSM_H
#define LMP_PAIR_CG_CMM_COUL_MSM_H

#include "pair_cmm_common.h"

namespace LAMMPS_NS {

class PairCGCMMCoulMSM : public PairCMMCommon {
 public:
  PairCGCMMCoulMSM(class LAMMPS *);
  ~PairCGCMMCoulMSM();

  void compute(int, int);
  void settings(int, char **);
  void init_style();
  double init_one(int, int);

  void write_restart(FILE *);
  void read_restart(FILE *);

  double memory_usage();

  void *extract(char *str);

 protected:
  void allocate();
  double _ia, _ia2, _ia3;
  int _smooth;
};

}

#endif
#endif
