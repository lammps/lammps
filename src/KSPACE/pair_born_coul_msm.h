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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(born/coul/msm,PairBornCoulMSM);
// clang-format on
#else

#ifndef LMP_PAIR_BORN_COUL_MSM_H
#define LMP_PAIR_BORN_COUL_MSM_H

#include "pair_born_coul_long.h"

namespace LAMMPS_NS {

class PairBornCoulMSM : public PairBornCoulLong {
 public:
  PairBornCoulMSM(class LAMMPS *);
  virtual ~PairBornCoulMSM();
  virtual void compute(int, int);
  virtual double single(int, int, int, int, double, double, double, double &);
  virtual void *extract(const char *, int &);

 protected:
  int nmax;
  double **ftmp;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Must use 'kspace_modify pressure/scalar no' to obtain per-atom virial with kspace_style MSM

The kspace scalar pressure option cannot be used to obtain per-atom virial.

*/
