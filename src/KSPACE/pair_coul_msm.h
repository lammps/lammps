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

#ifdef PAIR_CLASS

PairStyle(coul/msm,PairCoulMSM)

#else

#ifndef LMP_PAIR_COUL_MSM_H
#define LMP_PAIR_COUL_MSM_H

#include "pair_coul_long.h"

namespace LAMMPS_NS {

class PairCoulMSM : public PairCoulLong {
 public:
  PairCoulMSM(class LAMMPS *);
  virtual ~PairCoulMSM(){};
  virtual void compute(int, int);
  virtual double single(int, int, int, int, double, double, double, double &);
  virtual void *extract(const char *, int &);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Must use 'kspace_modify pressure/scalar no' to obtain per-atom virial with kspace_style MSM

The kspace scalar pressure option cannot be used to obtain per-atom virial.

*/
