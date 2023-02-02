// clang-format off
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

#include "mliap_data_kokkos.h"

#include "atom_kokkos.h"
#include "kokkos_type.h"
#include "pair_mliap_kokkos.h"
#include "atom_masks.h"
#include "mliap_descriptor.h"
#include "lammps.h"
#include "kokkos.h"

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template<class DeviceType>
MLIAPDataKokkos<DeviceType>::MLIAPDataKokkos(LAMMPS *lmp_in, int gradgradflag_in, int *map_in,
    class MLIAPModel* model_in,
    class MLIAPDescriptor* descriptor_in,
    class PairMLIAPKokkos<DeviceType>* pairmliap_in) :
    MLIAPData(lmp_in, gradgradflag_in, map_in, model_in, descriptor_in, pairmliap_in),
    k_pairmliap(pairmliap_in),
    lmp(lmp_in)
{
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
MLIAPDataKokkos<DeviceType>::~MLIAPDataKokkos() {
  memoryKK->destroy_kokkos(k_gradforce,gradforce);
  memoryKK->destroy_kokkos(k_betas,betas);
  memoryKK->destroy_kokkos(k_descriptors,descriptors);
  memoryKK->destroy_kokkos(k_eatoms,eatoms);
  memoryKK->destroy_kokkos(k_gamma_row_index,gamma_row_index);
  memoryKK->destroy_kokkos(k_gamma_col_index,gamma_col_index);
  memoryKK->destroy_kokkos(k_gamma,gamma);
  memoryKK->destroy_kokkos(k_iatoms,iatoms);
  memoryKK->destroy_kokkos(k_ielems,ielems);
  memoryKK->destroy_kokkos(k_numneighs,numneighs);
  memoryKK->destroy_kokkos(k_jatoms,jatoms);
  memoryKK->destroy_kokkos(k_pair_i,pair_i);
  memoryKK->destroy_kokkos(k_jelems,jelems);
  memoryKK->destroy_kokkos(k_elems,elems);
  memoryKK->destroy_kokkos(k_ij);
  memoryKK->destroy_kokkos(k_rij,rij);
  memoryKK->destroy_kokkos(k_graddesc,graddesc);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void MLIAPDataKokkos<DeviceType>::generate_neighdata(class NeighList *list_in, int eflag_in, int vflag_in) {

  list = list_in;

  // grow nmax gradforce array if necessary

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    memoryKK->destroy_kokkos(k_gradforce,gradforce);
    memoryKK->create_kokkos(k_gradforce, gradforce, nmax, size_gradforce, "mliap_data:gradforce");
    memoryKK->destroy_kokkos(k_elems,elems);
    memoryKK->create_kokkos(k_elems, elems, nmax, "mliap_data:elems");  }

  // clear gradforce array

  int nall = atom->nlocal + atom->nghost;
  ntotal = nall;
  auto d_gradforce = k_gradforce.template view<DeviceType>();
  Kokkos::deep_copy(d_gradforce, 0.);
  auto d_elems = k_elems.template view<DeviceType>();
  Kokkos::deep_copy(d_elems, 0.);
  // grow arrays if necessary

  nlistatoms = list->inum;
  if (nlistatoms_max < nlistatoms) {
    memoryKK->destroy_kokkos(k_betas,betas);
    memoryKK->create_kokkos(k_betas, betas, nlistatoms, ndescriptors, "mliap_data:betas");
    memoryKK->destroy_kokkos(k_descriptors,descriptors);
    memoryKK->create_kokkos(k_descriptors, descriptors, nlistatoms, ndescriptors, "mliap_data:descriptors");
    memoryKK->destroy_kokkos(k_eatoms,eatoms);
    memoryKK->create_kokkos(k_eatoms, eatoms, nlistatoms, "mliap_data:eatoms");
    nlistatoms_max = nlistatoms;
  }

  // grow gamma arrays if necessary

  if (gradgradflag == 1) {
    if (natomgamma_max < nlistatoms) {
      memoryKK->destroy_kokkos(k_gamma_row_index,gamma_row_index);
      memoryKK->create_kokkos(k_gamma_row_index, gamma_row_index, nlistatoms, gamma_nnz, "mliap_data:gamma_row_index");
      memoryKK->destroy_kokkos(k_gamma_col_index,gamma_col_index);
      memoryKK->create_kokkos(k_gamma_col_index, gamma_col_index, nlistatoms, gamma_nnz, "mliap_data:gamma_col_index");
      memoryKK->destroy_kokkos(k_gamma,gamma);
      memoryKK->create_kokkos(k_gamma, gamma, nlistatoms, gamma_nnz, "mliap_data:");
      natomgamma_max = nlistatoms;
    }
  }
  atomKK->sync(execution_space,X_MASK | V_MASK | F_MASK | TYPE_MASK);
  k_pairmliap->x = atomKK->k_x.view<DeviceType>();
  k_pairmliap->v = atomKK->k_v.view<DeviceType>();
  k_pairmliap->f = atomKK->k_f.view<DeviceType>();

  grow_neigharrays();

  // Use the ielems memory for prefix scan and set it at the end to the i type

  auto d_iatoms = k_iatoms.template view<DeviceType>();
  auto d_ielems = k_ielems.template view<DeviceType>();
  auto d_ij = k_ij.template view<DeviceType>();
  auto d_numneighs = k_numneighs.template view<DeviceType>();
  auto d_jatoms = k_jatoms.template view<DeviceType>();
  auto d_pair_i= k_pair_i.template view<DeviceType>();
  auto d_jelems= k_jelems.template view<DeviceType>();
  auto d_rij= k_rij.template view<DeviceType>();

  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  auto d_numneigh = k_list->d_numneigh;
  auto d_neighbors = k_list->d_neighbors;
  auto d_ilist = k_list->d_ilist;
  auto d_cutsq=k_pairmliap->k_cutsq.template view<DeviceType>();

  AtomKokkos *atomKK = (AtomKokkos *) atom;
  auto x = atomKK->k_x.view<DeviceType>();
  auto type = atomKK->k_type.view<DeviceType>();
  auto map=k_pairmliap->k_map.template view<DeviceType>();

  Kokkos::parallel_scan(nlistatoms, KOKKOS_LAMBDA (int ii, int &update, const bool final) {
    if (final)
      d_ij(ii) = update;
    update += d_numneighs(ii);
  });

  Kokkos::parallel_for(nlistatoms, KOKKOS_LAMBDA (int ii)  {
    int ij = d_ij(ii);
    const int i = d_ilist[ii];
    const double xtmp = x(i, 0);
    const double ytmp = x(i, 1);
    const double ztmp = x(i, 2);
    const int itype = type(i);
    const int ielem = map(itype);
    const int jnum = d_numneigh(i);
    for (int jj = 0; jj < jnum; jj++) {
      int j = d_neighbors(i,jj);
      j &= NEIGHMASK;
      const double delx = x(j,0) - xtmp;
      const double dely = x(j,1) - ytmp;
      const double delz = x(j,2) - ztmp;
      const double rsq = delx*delx + dely*dely + delz*delz;
      int jtype = type(j);
      const int jelem = map(jtype);
      if (rsq < d_cutsq(itype,jtype)) {
        d_jatoms(ij) = j;
        d_pair_i(ij) = i;
        d_jelems(ij) = jelem;
        d_rij(ij, 0) = delx;
        d_rij(ij, 1) = dely;
        d_rij(ij, 2) = delz;
        ij++;
      }
    }
    d_iatoms[ii] = i;
    d_ielems[ii] = ielem;
  });
  Kokkos::parallel_for(nmax, KOKKOS_LAMBDA (int i)  {
    const int itype = type(i);
    d_elems(i) = map(itype);
  });
  modified(execution_space, NUMNEIGHS_MASK | IATOMS_MASK | IELEMS_MASK | ELEMS_MASK | JATOMS_MASK | PAIR_I_MASK | JELEMS_MASK | RIJ_MASK | IJ_MASK );
  eflag = eflag_in;
  vflag = vflag_in;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void MLIAPDataKokkos<DeviceType>::grow_neigharrays() {
  AtomKokkos *atomKK = (AtomKokkos *) atom;
  f = atom->f;
  f_device = atomKK->k_f.view<DeviceType>().data();
  // grow neighbor arrays if necessary

  if (natomneigh_max < nlistatoms) {
    natomneigh_max = nlistatoms;

    memoryKK->destroy_kokkos(k_iatoms,iatoms);
    memoryKK->create_kokkos(k_iatoms, iatoms, natomneigh_max, "mliap_data:iatoms");
    memoryKK->destroy_kokkos(k_ielems,ielems);
    memoryKK->create_kokkos(k_ielems, ielems, natomneigh_max, "mliap_data:ielems");
    memoryKK->destroy_kokkos(k_ij);
    memoryKK->create_kokkos(k_ij, natomneigh_max, "mliap_data:ij");
    memoryKK->destroy_kokkos(k_numneighs,numneighs);
    memoryKK->create_kokkos(k_numneighs, numneighs, natomneigh_max, "mliap_data:numneighs");
  }

  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  auto d_numneigh = k_list->d_numneigh;
  auto d_neighbors = k_list->d_neighbors;
  auto d_ilist = k_list->d_ilist;

  auto x = atomKK->k_x.view<DeviceType>();
  auto type = atomKK->k_type.view<DeviceType>();
  auto d_cutsq=k_pairmliap->k_cutsq.template view<DeviceType>();
  auto h_cutsq=k_pairmliap->k_cutsq.template view<LMPHostType>();
  auto d_numneighs = k_numneighs.template view<DeviceType>();
  Kokkos::parallel_reduce(nlistatoms, KOKKOS_LAMBDA (int ii, int &contrib) {
    const int i = d_ilist[ii];
    int count=0;
    const double xtmp = x(i, 0);
    const double ytmp = x(i, 1);
    const double ztmp = x(i, 2);
    const int itype = type(i);
    const int jnum = d_numneigh(i);
    for (int jj = 0; jj < jnum; jj++) {
      int j = d_neighbors(i,jj);
      j &= NEIGHMASK;
      const double delx = x(j,0) - xtmp;
      const double dely = x(j,1) - ytmp;
      const double delz = x(j,2) - ztmp;
      const double rsq = delx*delx + dely*dely + delz*delz;
      int jtype = type(j);
      if (rsq < d_cutsq(itype,jtype))
        count++;
    }
    d_numneighs(ii) = count;
    contrib += count;
  }, npairs);
  modified(execution_space, NUMNEIGHS_MASK);

  if (nneigh_max < npairs) {
    memoryKK->destroy_kokkos(k_jatoms,jatoms);
    memoryKK->create_kokkos(k_jatoms, jatoms, npairs, "mliap_data:jatoms");
    memoryKK->destroy_kokkos(k_pair_i,pair_i);
    memoryKK->create_kokkos(k_pair_i, pair_i, npairs, "mliap_data:pair_i");
    memoryKK->destroy_kokkos(k_jelems,jelems);
    memoryKK->create_kokkos(k_jelems, jelems, npairs, "mliap_data:jelems");
    memoryKK->destroy_kokkos(k_rij,rij);
    memoryKK->create_kokkos(k_rij, rij, npairs, 3, "mliap_data:rij");

    if (gradgradflag == 0){
      memoryKK->destroy_kokkos(k_graddesc,graddesc);
      memoryKK->create_kokkos(k_graddesc, graddesc, npairs, ndescriptors,3, "mliap_data:graddesc");
    }
    nneigh_max = npairs;
   }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void MLIAPDataKokkos<DeviceType>::modified(ExecutionSpace space, unsigned int mask, bool ignore_auto_sync) {
  if (space == Device) {
    if (mask & IATOMS_MASK      ) k_iatoms         .modify<LMPDeviceType>();
    if (mask & IELEMS_MASK      ) k_ielems         .modify<LMPDeviceType>();
    if (mask & JATOMS_MASK      ) k_jatoms         .modify<LMPDeviceType>();
    if (mask & PAIR_I_MASK      ) k_pair_i         .modify<LMPDeviceType>();
    if (mask & JELEMS_MASK      ) k_jelems         .modify<LMPDeviceType>();
    if (mask & ELEMS_MASK       ) k_elems          .modify<LMPDeviceType>();
    if (mask & IJ_MASK          ) k_ij             .modify<LMPDeviceType>();
    if (mask & BETAS_MASK       ) k_betas          .modify<LMPDeviceType>();
    if (mask & DESCRIPTORS_MASK ) k_descriptors    .modify<LMPDeviceType>();
    if (mask & EATOMS_MASK      ) k_eatoms         .modify<LMPDeviceType>();
    if (mask & RIJ_MASK         ) k_rij            .modify<LMPDeviceType>();
    if (mask & GRADFORCE_MASK   ) k_gradforce      .modify<LMPDeviceType>();
    if (mask & GRADDESC_MASK    ) k_graddesc       .modify<LMPDeviceType>();
    if (mask & NUMNEIGHS_MASK   ) k_numneighs      .modify<LMPDeviceType>();
    if (mask & GAMMA_MASK_MASK  ) k_gamma          .modify<LMPDeviceType>();
    if (mask & GAMMA_ROW_MASK   ) k_gamma_row_index.modify<LMPDeviceType>();
    if (mask & GAMMA_COL_MASK   ) k_gamma_col_index.modify<LMPDeviceType>();

    if (lmp->kokkos->auto_sync && !ignore_auto_sync) sync(Host, mask, true);
  } else {
    if (mask & IATOMS_MASK      ) k_iatoms         .modify<LMPHostType>();
    if (mask & IELEMS_MASK      ) k_ielems         .modify<LMPHostType>();
    if (mask & JATOMS_MASK      ) k_jatoms         .modify<LMPHostType>();
    if (mask & PAIR_I_MASK      ) k_pair_i         .modify<LMPHostType>();
    if (mask & JELEMS_MASK      ) k_jelems         .modify<LMPHostType>();
    if (mask & ELEMS_MASK       ) k_elems          .modify<LMPHostType>();
    if (mask & IJ_MASK          ) k_ij             .modify<LMPHostType>();
    if (mask & BETAS_MASK       ) k_betas          .modify<LMPHostType>();
    if (mask & DESCRIPTORS_MASK ) k_descriptors    .modify<LMPHostType>();
    if (mask & EATOMS_MASK      ) k_eatoms         .modify<LMPHostType>();
    if (mask & RIJ_MASK         ) k_rij            .modify<LMPHostType>();
    if (mask & GRADFORCE_MASK   ) k_gradforce      .modify<LMPHostType>();
    if (mask & GRADDESC_MASK    ) k_graddesc       .modify<LMPHostType>();
    if (mask & NUMNEIGHS_MASK   ) k_numneighs      .modify<LMPHostType>();
    if (mask & GAMMA_MASK_MASK  ) k_gamma          .modify<LMPHostType>();
    if (mask & GAMMA_ROW_MASK   ) k_gamma_row_index.modify<LMPHostType>();
    if (mask & GAMMA_COL_MASK   ) k_gamma_col_index.modify<LMPHostType>();
    if (lmp->kokkos->auto_sync && !ignore_auto_sync) sync(Device, mask, true);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void MLIAPDataKokkos<DeviceType>::sync(ExecutionSpace space, unsigned int mask, bool ignore_auto_sync) {

  if (space == Device) {
    if (lmp->kokkos->auto_sync && !ignore_auto_sync) modified(Host, mask, true);
    if (mask & IATOMS_MASK      ) k_iatoms         .sync<LMPDeviceType>();
    if (mask & IELEMS_MASK      ) k_ielems         .sync<LMPDeviceType>();
    if (mask & JATOMS_MASK      ) k_jatoms         .sync<LMPDeviceType>();
    if (mask & PAIR_I_MASK      ) k_pair_i         .sync<LMPDeviceType>();
    if (mask & JELEMS_MASK      ) k_jelems         .sync<LMPDeviceType>();
    if (mask & ELEMS_MASK       ) k_elems          .sync<LMPDeviceType>();
    if (mask & IJ_MASK          ) k_ij             .sync<LMPDeviceType>();
    if (mask & BETAS_MASK       ) k_betas          .sync<LMPDeviceType>();
    if (mask & DESCRIPTORS_MASK ) k_descriptors    .sync<LMPDeviceType>();
    if (mask & EATOMS_MASK      ) k_eatoms         .sync<LMPDeviceType>();
    if (mask & RIJ_MASK         ) k_rij            .sync<LMPDeviceType>();
    if (mask & GRADFORCE_MASK   ) k_gradforce      .sync<LMPDeviceType>();
    if (mask & GRADDESC_MASK    ) k_graddesc       .sync<LMPDeviceType>();
    if (mask & NUMNEIGHS_MASK   ) k_numneighs      .sync<LMPDeviceType>();
    if (mask & GAMMA_MASK_MASK  ) k_gamma          .sync<LMPDeviceType>();
    if (mask & GAMMA_ROW_MASK   ) k_gamma_row_index.sync<LMPDeviceType>();
    if (mask & GAMMA_COL_MASK   ) k_gamma_col_index.sync<LMPDeviceType>();
  } else {
    if (lmp->kokkos->auto_sync && !ignore_auto_sync) modified(Device, mask, true);
    if (mask & IATOMS_MASK      ) k_iatoms         .sync<LMPHostType>();
    if (mask & IELEMS_MASK      ) k_ielems         .sync<LMPHostType>();
    if (mask & JATOMS_MASK      ) k_jatoms         .sync<LMPHostType>();
    if (mask & PAIR_I_MASK      ) k_pair_i         .sync<LMPHostType>();
    if (mask & JELEMS_MASK      ) k_jelems         .sync<LMPHostType>();
    if (mask & ELEMS_MASK       ) k_elems          .sync<LMPHostType>();
    if (mask & IJ_MASK          ) k_ij             .sync<LMPHostType>();
    if (mask & BETAS_MASK       ) k_betas          .sync<LMPHostType>();
    if (mask & DESCRIPTORS_MASK ) k_descriptors    .sync<LMPHostType>();
    if (mask & EATOMS_MASK      ) k_eatoms         .sync<LMPHostType>();
    if (mask & RIJ_MASK         ) k_rij            .sync<LMPHostType>();
    if (mask & GRADFORCE_MASK   ) k_gradforce      .sync<LMPHostType>();
    if (mask & GRADDESC_MASK    ) k_graddesc       .sync<LMPHostType>();
    if (mask & NUMNEIGHS_MASK   ) k_numneighs      .sync<LMPHostType>();
    if (mask & GAMMA_MASK_MASK  ) k_gamma          .sync<LMPHostType>();
    if (mask & GAMMA_ROW_MASK   ) k_gamma_row_index.sync<LMPHostType>();
    if (mask & GAMMA_COL_MASK   ) k_gamma_col_index.sync<LMPHostType>();
  }
}

/* ---------------------------------------------------------------------- */

template class MLIAPDataKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class MLIAPDataKokkos<LMPHostType>;
#endif
}// namespace
