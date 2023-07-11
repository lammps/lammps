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

#ifndef LMP_MLIAP_MODEL_PYTHON_H
#define LMP_MLIAP_MODEL_PYTHON_H

#include "mliap_model.h"

namespace LAMMPS_NS {

class MLIAPModelPython : public MLIAPModel {
 public:
  MLIAPModelPython(LAMMPS *, char * = nullptr, bool is_child=false);
  ~MLIAPModelPython() override;
  int get_nparams() override;
  int get_gamma_nnz(class MLIAPData *) override;
  void compute_gradients(class MLIAPData *) override;
  void compute_gradgrads(class MLIAPData *) override;
  void compute_force_gradients(class MLIAPData *) override;
  double memory_usage() override;
  void connect_param_counts();    // If possible convert this to protected/private and
                                  // and figure out how to declare cython fn
                                  // load_from_python as a friend.
  int model_loaded;

 protected:
  void read_coeffs(char *) override;

 private:
};

}    // namespace LAMMPS_NS

#endif
