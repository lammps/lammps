/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_MLIAP_MODEL_NN_H
#define LMP_MLIAP_MODEL_NN_H

#include "mliap_model.h"

#include <cmath>

namespace LAMMPS_NS {

class MLIAPModelNN : public MLIAPModel {
 public:
  MLIAPModelNN(LAMMPS *, char * = nullptr);
  ~MLIAPModelNN() override;
  int get_nparams() override;
  int get_gamma_nnz(class MLIAPData *) override;
  void compute_gradients(class MLIAPData *) override;
  void compute_gradgrads(class MLIAPData *) override;
  void compute_force_gradients(class MLIAPData *) override;
  double memory_usage() override;

  int nlayers;    // number of layers per element

 protected:
  int *activation;    // activation functions
  int *nnodes;        // number of nodes per layer
  double ***scale;    // element scale values
  void read_coeffs(char *) override;

  inline double sigm(double x, double &deriv)
  {
    double expl = 1. / (1. + exp(-x));
    deriv = expl * (1 - expl);
    return expl;
  }

  inline double tanh(double x, double &deriv)
  {
    double expl = 2. / (1. + exp(-2. * x)) - 1;
    deriv = 1. - expl * expl;
    return expl;
  }

  inline double relu(double x, double &deriv)
  {
    if (x > 0) {
      deriv = 1.;
      return x;
    } else {
      deriv = 0.;
      return 0;
    }
  }
};

}    // namespace LAMMPS_NS

#endif
