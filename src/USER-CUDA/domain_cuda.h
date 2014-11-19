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

#ifndef LMP_DOMAIN_CUDA_H
#define LMP_DOMAIN_CUDA_H

#include "pointers.h"
#include "domain.h"

namespace LAMMPS_NS {

class DomainCuda : public Domain {
 public:
  DomainCuda(class LAMMPS *);
  void init();
  void set_global_box();
  void set_lamda_box();
  void set_local_box();
  void reset_box();
  void pbc();

  virtual void lamda2x(int);
  virtual void x2lamda(int);

 protected:
  class Cuda *cuda;
};

}

#endif
