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

/* ----------------------------------------------------------------------
   Contributing author: Matt Bettencourt (NVIDIA)
 ------------------------------------------------------------------------- */

#ifndef LMP_MLIAP_MODEL_KOKKOS_H
#define LMP_MLIAP_MODEL_KOKKOS_H

#include "kokkos_type.h"
#include "memory_kokkos.h"
#include "pointers.h"

namespace LAMMPS_NS {

template <class DeviceType> class MLIAPModelKokkos : protected Pointers {
 public:
  MLIAPModelKokkos(LAMMPS *lmp, MLIAPModel *model_in) : Pointers(lmp), model(model_in) {}
  virtual ~MLIAPModelKokkos()
  {
    memoryKK->destroy_kokkos(k_coeffelem);
    model->coeffelem = nullptr;
  }

  void set_k_coeffelem()
  {
    double **tmp = nullptr;
    memoryKK->destroy_kokkos(k_coeffelem);
    memoryKK->create_kokkos(k_coeffelem, tmp, model->nelements, model->nparams,
                            "MLIAPModelKokkos::coeffelem");
    for (int i = 0; i < model->nelements; ++i)
      for (int j = 0; j < model->nparams; ++j) tmp[i][j] = model->coeffelem[i][j];
    delete model->coeffelem;
    model->coeffelem = tmp;
    k_coeffelem.modify<LMPHostType>();
    k_coeffelem.sync<LMPDeviceType>();
  }

  MLIAPModel *model;
  DAT::tdual_float_2d k_coeffelem;
};
}    // namespace LAMMPS_NS
#endif
