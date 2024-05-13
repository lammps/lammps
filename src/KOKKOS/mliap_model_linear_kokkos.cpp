// clang-format off
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

#include "mliap_model_linear_kokkos.h"

#include "mliap_data_kokkos.h"

using namespace LAMMPS_NS;

template<class DeviceType>
MLIAPModelLinearKokkos<DeviceType>::MLIAPModelLinearKokkos(LAMMPS *lmp, char *args) :
  MLIAPModelLinear(lmp,args),
  MLIAPModelKokkos<DeviceType>(lmp, this)
{
  if (args) MLIAPModelKokkos<DeviceType>::set_k_coeffelem();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void MLIAPModelLinearKokkos<DeviceType>::compute_gradients(class MLIAPData *data)
{
  MLIAPDataKokkos<DeviceType> *k_data = (MLIAPDataKokkos<DeviceType>*)(data);

  // read but never changes
  auto d_coeffelem = this->k_coeffelem.template view<DeviceType>();

  // read
  auto d_ielems = k_data->k_ielems.template view<DeviceType>();
  auto d_descriptors = k_data->k_descriptors.template view<DeviceType>();

  // written
  auto d_betas = k_data->k_betas.template view<DeviceType>();
  auto d_eatoms = k_data->k_eatoms.template view<DeviceType>();

  const auto eflag = data->eflag;
  const int ndescriptors=data->ndescriptors;

  Kokkos::parallel_reduce(data->nlistatoms, KOKKOS_LAMBDA (int ii, double &update) {
    const int ielem = d_ielems(ii);

    for (int icoeff = 0; icoeff < ndescriptors; icoeff++)
      d_betas(ii,icoeff) = d_coeffelem(ielem,icoeff+1);

    // add in contributions to global and per-atom energy
    // this is optional and has no effect on force calculation
    if (eflag) {
      // energy of atom I
      double etmp = d_coeffelem(ielem,0);

      // E_i = beta.B_i
      for (int icoeff = 0; icoeff < ndescriptors; icoeff++)
        etmp += d_coeffelem(ielem,icoeff+1)*d_descriptors(ii, icoeff);
      update += etmp;
      d_eatoms(ii) = etmp;
    }
  }, data->energy);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void MLIAPModelLinearKokkos<DeviceType>::compute_gradgrads(class MLIAPData *data)
{
  MLIAPModelLinear::compute_gradgrads(data);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void MLIAPModelLinearKokkos<DeviceType>::compute_force_gradients(class MLIAPData *data)
{
  MLIAPModelLinear::compute_force_gradients(data);
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class MLIAPModelLinearKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class MLIAPModelLinearKokkos<LMPHostType>;
#endif
}
