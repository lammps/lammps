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

// NOTE: this file is *supposed* to be included multiple times

#ifdef LMP_USER_OMP

// true interface to USER-OMP

// provide a DomainOMP class with some overrides for Domain
#include "domain.h"

#ifndef LMP_DOMAIN_OMP_H
#define LMP_DOMAIN_OMP_H

namespace LAMMPS_NS {

class DomainOMP : public Domain {
 public:
  DomainOMP(class LAMMPS *lmp) : Domain(lmp) {}
  virtual ~DomainOMP() {}

  // multi-threaded versions
  virtual void pbc();
  virtual void lamda2x(int);
  virtual void lamda2x(double *lamda, double *x) {Domain::lamda2x(lamda,x);}
  virtual void x2lamda(int);
  virtual void x2lamda(double *x, double *lamda) {Domain::x2lamda(x,lamda);}
};
}

#endif /* LMP_DOMAIN_OMP_H */

#endif /* !LMP_USER_OMP */
