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

#ifdef IMPROPER_CLASS

ImproperStyle(ring/omp,ImproperRingOMP)

#else

#ifndef LMP_IMPROPER_RING_OMP_H
#define LMP_IMPROPER_RING_OMP_H

#include "improper_ring.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class ImproperRingOMP : public ImproperRing, public ThrOMP {

 public:
  ImproperRingOMP(class LAMMPS *lmp);
  virtual void compute(int, int);

 private:
  template <int EVFLAG, int EFLAG, int NEWTON_BOND>
  void eval(int ifrom, int ito, ThrData * const thr);
};

}

#endif
#endif
