// clang-format off
/* ----------------------------------------------------------------------
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

#include "mliap_descriptor_so3_kokkos.h"

#include "atom_kokkos.h"
#include "comm.h"
#include "error.h"
#include "memory.h"
#include "mliap_data_kokkos.h"
#include "mliap_so3_kokkos.h"
#include "pair_mliap.h"
#include "tokenizer.h"

#include <cstring>

using namespace LAMMPS_NS;

static constexpr int MAXLINE = 1024;
static constexpr int MAXWORD = 3;

/* ---------------------------------------------------------------------- */
template <class DeviceType>
MLIAPDescriptorSO3Kokkos<DeviceType>::MLIAPDescriptorSO3Kokkos(LAMMPS *lmp, char *paramfilename)
 // TODO: why take self as param, shouldn't be needed
  :  Pointers(lmp), MLIAPDescriptorSO3(lmp, paramfilename), MLIAPDescriptorKokkos<DeviceType>(lmp, this)
{
  // TODO: the MLIAP_SO3 object likely needs a kokkos-ified version
  so3ptr_kokkos = new MLIAP_SO3Kokkos<DeviceType>(lmp, rcutfac, lmax, nmax, alpha);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void MLIAPDescriptorSO3Kokkos<DeviceType>::compute_descriptors(class MLIAPData *data_)
{
  auto data = static_cast<MLIAPDataKokkos<DeviceType>*>(data_);
  so3ptr_kokkos->spectrum(data->nlistatoms, data->k_numneighs, data->k_jelems, this->k_wjelem, data->k_rij, data->k_ij,
                          nmax, lmax, rcutfac, alpha, data->npairs, data->ndescriptors);

  Kokkos::deep_copy(data->k_descriptors.template view<DeviceType>(), so3ptr_kokkos->m_plist_r);
  Kokkos::deep_copy(data->k_descriptors.h_view, so3ptr_kokkos->m_plist_r);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void MLIAPDescriptorSO3Kokkos<DeviceType>::compute_forces(class MLIAPData *data_)
{
  auto data = static_cast<MLIAPDataKokkos<DeviceType>*>(data_);
  int npairs = data->npairs;
  auto d_numneighs = data->k_numneighs.template view<DeviceType>();
  so3ptr_kokkos->spectrum_dxdr(data->nlistatoms, data->k_numneighs, data->k_jelems, this->k_wjelem, data->k_rij, data->k_ij,
                               nmax, lmax, rcutfac, alpha, npairs, data->ndescriptors);

  auto d_f = atomKK->k_f.view<DeviceType>();
  auto d_iatoms = data->k_iatoms.template view<DeviceType>();
  auto d_jatoms = data->k_jatoms.template view<DeviceType>();
  auto d_betas = data->k_betas.template view<DeviceType>();
  auto d_rij = data->k_rij.template view<DeviceType>();
  auto d_ij = data->k_ij.template view<DeviceType>();
  auto ndescriptors = data->ndescriptors;
  auto d_dplist_r = so3ptr_kokkos->k_dplist_r;
  auto vflag=data->vflag;
  int vflag_either=data->k_pairmliap->vflag_either, vflag_global=data->pairmliap->vflag_global, vflag_atom=data->pairmliap->vflag_atom;
  auto d_vatom = data->k_pairmliap->k_vatom.template view<DeviceType>();
  Kokkos::View<double[6], DeviceType> virial("virial");
  data->k_pairmliap->k_vatom.template modify<LMPHostType>();
  data->k_pairmliap->k_vatom.template sync<DeviceType>();
  Kokkos::parallel_for(data->nlistatoms, KOKKOS_LAMBDA(int ii) {
    double fij[3];
    const int i = d_iatoms(ii);

    // ensure rij, inside, wj, and rcutij are of size jnum

    const int jnum = d_numneighs(ii);
    int ij = d_ij(ii); // use precomputed ij to break loop dependency
    for (int jj = 0; jj < jnum; jj++) {
      int j = d_jatoms(ij);

      for (int ir = 0; ir < 3; ir++) {
        fij[ir] = 0.0;
        for (int icoeff = 0; icoeff < ndescriptors; icoeff++) {
          fij[ir] += d_betas(ii, icoeff) *
              d_dplist_r(ij,icoeff, ir);
        }
      }

      Kokkos::atomic_add(&d_f(i, 0),fij[0]);
      Kokkos::atomic_add(&d_f(i, 1),fij[1]);
      Kokkos::atomic_add(&d_f(i, 2),fij[2]);
      Kokkos::atomic_add(&d_f(j, 0),-fij[0]);
      Kokkos::atomic_add(&d_f(j, 1),-fij[1]);
      Kokkos::atomic_add(&d_f(j, 2),-fij[2]);

      if (vflag) {
        v_tally(vflag_either, vflag_global, vflag_atom, i, j, ij, fij, d_rij, virial, d_vatom);
      }
      ij++;
    }
  });

  if (vflag) {
    if (vflag_global) {
      Kokkos::View<double[6], LMPHostType> h_virial("h_virial");
      Kokkos::deep_copy(h_virial,virial);
      for (int i=0;i<6;++i)
        data->k_pairmliap->virial[i]+=h_virial[i];
    }
    if (vflag_atom) {
      data->k_pairmliap->k_vatom.template modify<DeviceType>();
      data->k_pairmliap->k_vatom.template sync<LMPHostType>();
    }
  }
}

/* ----------------------------------------------------------------------
   add virial contribution into global and per-atom accumulators
------------------------------------------------------------------------- */
template <class DeviceType>
template <typename ViewType>
KOKKOS_INLINE_FUNCTION
void MLIAPDescriptorSO3Kokkos<DeviceType>::v_tally(int vflag_either, int vflag_global, int vflag_atom, int i, int j, int ij,
    double *fij, ViewType rij, Kokkos::View<double[6],DeviceType> virial, ViewType vatom)
{
  double v[6];
  if (vflag_either) {
    v[0] = -rij(ij,0)*fij[0];
    v[1] = -rij(ij,1)*fij[1];
    v[2] = -rij(ij,2)*fij[2];
    v[3] = -rij(ij,0)*fij[1];
    v[4] = -rij(ij,0)*fij[2];
    v[5] = -rij(ij,1)*fij[2];
    if (vflag_global) {
      Kokkos::atomic_add(&virial[0], v[0]);
      Kokkos::atomic_add(&virial[1], v[1]);
      Kokkos::atomic_add(&virial[2], v[2]);
      Kokkos::atomic_add(&virial[3], v[3]);
      Kokkos::atomic_add(&virial[4], v[4]);
      Kokkos::atomic_add(&virial[5], v[5]);
    }
    if (vflag_atom) {
      Kokkos::atomic_add(&vatom(i,0), 0.5*v[0]);
      Kokkos::atomic_add(&vatom(i,1), 0.5*v[1]);
      Kokkos::atomic_add(&vatom(i,2), 0.5*v[2]);
      Kokkos::atomic_add(&vatom(i,3), 0.5*v[3]);
      Kokkos::atomic_add(&vatom(i,4), 0.5*v[4]);
      Kokkos::atomic_add(&vatom(i,5), 0.5*v[5]);

      Kokkos::atomic_add(&vatom(j,0), 0.5*v[0]);
      Kokkos::atomic_add(&vatom(j,1), 0.5*v[1]);
      Kokkos::atomic_add(&vatom(j,2), 0.5*v[2]);
      Kokkos::atomic_add(&vatom(j,3), 0.5*v[3]);
      Kokkos::atomic_add(&vatom(j,4), 0.5*v[4]);
      Kokkos::atomic_add(&vatom(j,5), 0.5*v[5]);
    }
  }
}
/* ---------------------------------------------------------------------- */

template <class DeviceType>
void MLIAPDescriptorSO3Kokkos<DeviceType>::compute_force_gradients(class MLIAPData *data_)
{
  error->all(FLERR,"This has not been tested in cuda/kokkos");

  auto data = static_cast<MLIAPDataKokkos<DeviceType>*>(data_);
  int npairs = data->npairs;
  so3ptr_kokkos->spectrum_dxdr(data->nlistatoms, data->k_numneighs, data->k_jelems, this->k_wjelem, data->k_rij, data->k_ij,
                               nmax, lmax, rcutfac, alpha, npairs, data->ndescriptors);
  auto d_dplist_r = so3ptr_kokkos->k_dplist_r;
  auto d_gradforce = data->k_gradforce.template view<DeviceType>();
  auto d_gamma = data->k_gamma.template view<DeviceType>();
  auto d_gamma_row_index = data->k_gamma_row_index.template view<DeviceType>();
  auto d_gamma_col_index = data->k_gamma_col_index.template view<DeviceType>();
  auto d_jatoms = data->k_jatoms.template view<DeviceType>();
  auto d_ij = data->k_ij.template view<DeviceType>();
  auto d_numneighs = data->k_numneighs.template view<DeviceType>();
  auto d_iatoms = data->k_iatoms.template view<DeviceType>();

  auto yoffset = data->yoffset, zoffset = data->zoffset, gamma_nnz = data->gamma_nnz;

  Kokkos::parallel_for (data->nlistatoms, KOKKOS_LAMBDA (int ii) {
    const int i = d_iatoms(ii);

    // ensure rij, inside, wj, and rcutij are of size jnum

    const int jnum = d_numneighs(ii);
    int ij = d_ij(ii);
    for (int jj = 0; jj < jnum; jj++) {
      int j = d_jatoms(ij);

      for (int inz = 0; inz < gamma_nnz; inz++) {
        const int l = d_gamma_row_index(ii, inz);
        const int k = d_gamma_col_index(ii, inz);

        Kokkos::atomic_add(&d_gradforce(i, l),  +
            d_gamma(ii, inz) * d_dplist_r(ij, k, 0));
        Kokkos::atomic_add(&d_gradforce(i, l + yoffset), +
            d_gamma(ii, inz) * d_dplist_r(ij, k, 1));
        Kokkos::atomic_add(&d_gradforce(i, l + zoffset), +
            d_gamma(ii, inz) * d_dplist_r(ij, k, 2));
        Kokkos::atomic_add(&d_gradforce(j, l), -
            d_gamma(ii, inz) * d_dplist_r(ij, k, 0));
        Kokkos::atomic_add(&d_gradforce(j, l + yoffset), -
            d_gamma(ii, inz) * d_dplist_r(ij, k, 1));
        Kokkos::atomic_add(&d_gradforce(j, l + zoffset), -
            d_gamma(ii, inz) * d_dplist_r(ij, k, 2));
      }
      ij++;
    }
  });
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void MLIAPDescriptorSO3Kokkos<DeviceType>::compute_descriptor_gradients(class MLIAPData *data_)
{
  auto data = static_cast<MLIAPDataKokkos<DeviceType>*>(data_);
  bigint npairs = data->npairs;
  so3ptr_kokkos->spectrum_dxdr(data->nlistatoms, data->k_numneighs, data->k_jelems, this->k_wjelem, data->k_rij, data->k_ij,
                               nmax, lmax, rcutfac, alpha, npairs, data->ndescriptors);
  auto graddesc = data->k_graddesc.template view<DeviceType>();
  Kokkos::deep_copy(graddesc, so3ptr_kokkos->k_dplist_r);
  Kokkos::deep_copy(data->k_graddesc.h_view, graddesc);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void MLIAPDescriptorSO3Kokkos<DeviceType>::init()
{
  so3ptr_kokkos->init();
  MLIAPDescriptorKokkos<DeviceType>::init_data();
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
double MLIAPDescriptorSO3Kokkos<DeviceType>::memory_usage()
{
  double bytes = MLIAPDescriptor::memory_usage();
  bytes += so3ptr_kokkos->memory_usage();

  return bytes;
}

namespace LAMMPS_NS {
template class MLIAPDescriptorSO3Kokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class MLIAPDescriptorSO3Kokkos<LMPHostType>;
#endif
}
