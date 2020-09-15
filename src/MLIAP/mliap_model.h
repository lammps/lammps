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

#ifndef LMP_MLIAP_MODEL_H
#define LMP_MLIAP_MODEL_H

#include "pointers.h"

namespace LAMMPS_NS {

class MLIAPModel : protected Pointers {
public:
  MLIAPModel(LAMMPS*, char*);
  ~MLIAPModel();
  void set_ndescriptors(int);
  void set_nelements(int);
  virtual int get_nparams()=0;
  virtual int get_gamma_nnz(class MLIAPData*)=0;
  virtual void compute_gradients(class MLIAPData*)=0;
  virtual void compute_gradgrads(class MLIAPData*)=0;
  virtual void compute_force_gradients(class MLIAPData*)=0;
  virtual void init();
  virtual double memory_usage();
  int nelements;                 // # of unique elements
  int nonlinearflag;             // 1 if gradient() requires escriptors
  int ndescriptors;              // number of descriptors
  int nparams;                   // number of parameters per element

protected:
  void read_coeffs(char *);
  double **coeffelem;            // element coefficients
};

}

#endif

