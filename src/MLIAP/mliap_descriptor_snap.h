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

#ifndef LMP_MLIAP_DESCRIPTOR_SNAP_H
#define LMP_MLIAP_DESCRIPTOR_SNAP_H

#include "mliap_descriptor.h"

namespace LAMMPS_NS {

class MLIAPDescriptorSNAP : public MLIAPDescriptor  {
public:
  MLIAPDescriptorSNAP(LAMMPS*, char*, class PairMLIAP*);
  ~MLIAPDescriptorSNAP();
  virtual void forward(class NeighList*, double**);
  virtual void backward(class NeighList*, double**, int);
  virtual void init();
  virtual double get_cutoff(int, int);
  virtual double memory_usage();

  double rcutfac;                // declared public to workaround gcc 4.9
                                 // compiler bug, manifest in KOKKOS package
protected:
  class SNA* snaptr;
  void read_paramfile(char *);
  inline int equal(double* x,double* y);
  inline double dist2(double* x,double* y);

  double *radelem;              // element radii
  double *wjelem;               // elements weights
  int twojmax, switchflag, bzeroflag, bnormflag;
  int alloyflag, wselfallflag;
  double rfac0, rmin0;
};

}

#endif

