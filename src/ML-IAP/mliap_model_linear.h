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

#ifndef LMP_MLIAP_MODEL_LINEAR_H
#define LMP_MLIAP_MODEL_LINEAR_H

#include "mliap_model.h"

namespace LAMMPS_NS {

class MLIAPModelLinear : public MLIAPModelSimple {
 public:
  MLIAPModelLinear(LAMMPS *, char * = nullptr);

  int get_nparams() override;
  int get_gamma_nnz(class MLIAPData *) override;
  void compute_gradients(class MLIAPData *) override;
  void compute_gradgrads(class MLIAPData *) override;
  void compute_force_gradients(class MLIAPData *) override;

 protected:
};

}    // namespace LAMMPS_NS

#endif
