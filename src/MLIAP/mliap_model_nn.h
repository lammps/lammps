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

#ifndef LMP_MLIAP_MODEL_NN_H
#define LMP_MLIAP_MODEL_NN_H

#include "mliap_model.h"
#include <math.h>

namespace LAMMPS_NS {

class MLIAPModelNN : public MLIAPModelSimple {
public:
  MLIAPModelNN(LAMMPS*, char* = nullptr);
  ~MLIAPModelNN();
  virtual int get_nparams();
  virtual int get_gamma_nnz(class MLIAPData*);
  virtual void compute_gradients(class MLIAPData*);
  virtual void compute_gradgrads(class MLIAPData*);
  virtual void compute_force_gradients(class MLIAPData*);

protected:
};

}

static inline double sigm(double x, double &deriv) {
    double expl = 1./(1.+exp(-x));
    deriv = expl*(1-expl);
    return expl;
}

static inline double tanh(double x, double &deriv) {
    double expl = 2./(1.+exp(-2.*x))-1;
    deriv = 1.-expl*expl;
    return expl;
}

static inline double relu(double x, double &deriv) {
    if (x > 0) {
        deriv = 1.;
        return x;
    } else {
        deriv = 0.;
        return 0;
    }
}

#endif

