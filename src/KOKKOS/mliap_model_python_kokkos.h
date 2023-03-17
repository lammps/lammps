/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Matt Bettencourt (NVIDIA)
 ------------------------------------------------------------------------- */

#ifndef LMP_MLIAP_MODEL_PYTHON_KOKKOS_H
#define LMP_MLIAP_MODEL_PYTHON_KOKKOS_H

#include "mliap_model_python.h"
#include "mliap_model_kokkos.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template <class DeviceType>
class MLIAPModelPythonKokkos : public MLIAPModelPython, public MLIAPModelKokkos<DeviceType> {
 public:
  MLIAPModelPythonKokkos(LAMMPS *, char * = nullptr);
  ~MLIAPModelPythonKokkos();
  void read_coeffs(char *fname);

  void compute_gradients(class MLIAPData *) override;
  void compute_gradgrads(class MLIAPData *) override;
  void compute_force_gradients(class MLIAPData *) override;
  void connect_param_counts();
};

class  MLIAPModelPythonKokkosDevice: public MLIAPModelPythonKokkos<LMPDeviceType> {
};

}    // namespace LAMMPS_NS


#endif
