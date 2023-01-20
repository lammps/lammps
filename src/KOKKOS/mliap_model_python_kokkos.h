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
}    // namespace LAMMPS_NS




#include "mliap_data_kokkos.h"

namespace LAMMPS_NS {
class  MLIAPModelPythonKokkosDevice: public MLIAPModelPythonKokkos<LMPDeviceType> {
};

class MLIAPDataKokkosDevice {
public:

  MLIAPDataKokkosDevice(MLIAPDataKokkos<LMPDeviceType> &base) :
    ndescriptors(base.ndescriptors),
    nlistatoms(base.nlistatoms),
    ielems(base.k_ielems.d_view.data()),
    descriptors(base.k_descriptors.d_view.data()),
    betas(base.k_betas.d_view.data()),
    eatoms(base.k_eatoms.d_view.data()),
    energy(&base.energy),
#if defined(KOKKOS_ENABLE_CUDA)
    dev(1)
#else
    dev(0)
#endif
    {  }

  const int ndescriptors;
  const int nlistatoms;
  int *ielems;
  double *descriptors;
  double *betas;
  double *eatoms;
  double *energy;
  int dev;

#ifdef LMP_KOKKOS_GPU
  MLIAPDataKokkosDevice(MLIAPDataKokkos<LMPHostType> &base) : ndescriptors(-1),nlistatoms(-1)
  {
    // It cannot get here, but needed for compilation
  }
#endif
};
}

#endif
