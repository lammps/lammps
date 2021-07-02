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

#ifndef LMP_MLIAP_MODEL_QUADRATIC_H
#define LMP_MLIAP_MODEL_QUADRATIC_H

#include "mliap_model.h"

namespace LAMMPS_NS {

class MLIAPModelQuadratic : public MLIAPModelSimple {
 public:
  MLIAPModelQuadratic(LAMMPS *, char * = nullptr);
  ~MLIAPModelQuadratic();
  virtual int get_nparams();
  virtual int get_gamma_nnz(class MLIAPData *);
  virtual void compute_gradients(class MLIAPData *);
  virtual void compute_gradgrads(class MLIAPData *);
  virtual void compute_force_gradients(class MLIAPData *);

 protected:
};

}    // namespace LAMMPS_NS

#endif
