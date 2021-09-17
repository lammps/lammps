/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_MLIAP_DESCRIPTOR_SO3_H
#define LMP_MLIAP_DESCRIPTOR_SO3_H

#include "mliap_descriptor.h"

namespace LAMMPS_NS {

class MLIAPDescriptorSO3 : public MLIAPDescriptor {

 public:
  MLIAPDescriptorSO3(LAMMPS *, char *);
  virtual ~MLIAPDescriptorSO3();

  virtual void compute_descriptors(class MLIAPData *);
  virtual void compute_forces(class MLIAPData *);
  virtual void compute_force_gradients(class MLIAPData *){};
  virtual void compute_descriptor_gradients(class MLIAPData *){};
  virtual void init();
  virtual double memory_usage();

  double rcutfac;

 protected:
  class MLIAP_SO3 *so3ptr;
  void read_paramfile(char *);
  inline int equal(double *x, double *y);
  inline double dist2(double *x, double *y);

  int nmax, lmax;
  double alpha;

  int twojmax, switchflag, bzeroflag;
  int chemflag, bnormflag, wselfallflag;
  double rfac0, rmin0;
};
}    // namespace LAMMPS_NS

#endif
