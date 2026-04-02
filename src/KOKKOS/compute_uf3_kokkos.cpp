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

#include "compute_uf3_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "kokkos.h"
#include "memory.h"
#include "my_page.h"
#include "neigh_list_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "update.h"
#include "utils.h"

#include <algorithm>

using namespace LAMMPS_NS;

namespace {
static constexpr int PGDELTA = 1;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
ComputeUF3Kokkos<DeviceType>::ComputeUF3Kokkos(LAMMPS *lmp, int narg, char **arg) :
    ComputeUF3(lmp, narg, arg), atomKK(nullptr), kk_mirror_nmax(0), kk_mirror_ilist(nullptr),
    kk_mirror_numneigh(nullptr), kk_mirror_firstneigh(nullptr), kk_mirror_ipage(nullptr)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | TYPE_MASK | MASK_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType> ComputeUF3Kokkos<DeviceType>::~ComputeUF3Kokkos()
{
  if (copymode) return;
  free_kk_mirror();
}

/* ---------------------------------------------------------------------- */

template <class DeviceType> void ComputeUF3Kokkos<DeviceType>::free_kk_mirror()
{
  memory->destroy(kk_mirror_ilist);
  memory->destroy(kk_mirror_numneigh);
  memory->sfree(kk_mirror_firstneigh);
  kk_mirror_ilist = nullptr;
  kk_mirror_numneigh = nullptr;
  kk_mirror_firstneigh = nullptr;
  kk_mirror_nmax = 0;
  delete[] kk_mirror_ipage;
  kk_mirror_ipage = nullptr;
}

/* ---------------------------------------------------------------------- */

template <class DeviceType> void ComputeUF3Kokkos<DeviceType>::mirror_kokkos_neighbors_to_cpu()
{
  auto *k_list = static_cast<NeighListKokkos<DeviceType> *>(list);
  k_list->k_ilist.sync_host();

  int inum_all = list->inum;
  if (list->ghost) inum_all += list->gnum;

  auto h_ilist = k_list->k_ilist.view_host();
  auto h_numneigh = Kokkos::create_mirror_view_and_copy(LMPHostType(), k_list->d_numneigh);
  auto h_neighbors = Kokkos::create_mirror_view_and_copy(LMPHostType(), k_list->d_neighbors);

  const int nmax_atoms = atom->nmax;
  const int need = std::max(nmax_atoms, inum_all + 1);
  if (need > kk_mirror_nmax) {
    memory->destroy(kk_mirror_ilist);
    memory->destroy(kk_mirror_numneigh);
    memory->sfree(kk_mirror_firstneigh);
    memory->create(kk_mirror_ilist, need, "compute_uf3_kk:mirror_ilist");
    memory->create(kk_mirror_numneigh, need, "compute_uf3_kk:mirror_numneigh");
    memory->create(kk_mirror_firstneigh, need, "compute_uf3_kk:mirror_firstneigh");
    kk_mirror_nmax = need;
  }

  if (kk_mirror_ipage == nullptr) {
    const int nmypage = comm->nthreads;
    kk_mirror_ipage = new MyPage<int>[nmypage];
    for (int i = 0; i < nmypage; i++)
      kk_mirror_ipage[i].init(neighbor->oneatom, neighbor->pgsize, PGDELTA);
  }

  for (int t = 0; t < comm->nthreads; t++) kk_mirror_ipage[t].reset();

  int *ilist = kk_mirror_ilist;
  int *numneigh = kk_mirror_numneigh;
  int **firstneigh = kk_mirror_firstneigh;
  MyPage<int> *ipage = kk_mirror_ipage;

  for (int ii = 0; ii < inum_all; ii++) {
    int *neighptr = ipage[0].vget();

    const int i = h_ilist[ii];
    ilist[ii] = i;

    const int jnum = h_numneigh[i];
    numneigh[i] = jnum;

    for (int jj = 0; jj < jnum; jj++) neighptr[jj] = h_neighbors(i, jj);

    firstneigh[i] = neighptr;
    ipage[0].vgot(jnum);
    if (ipage[0].status())
      error->one(FLERR, Error::NOLASTLINE,
                 "UF3 Kokkos neighbor mirror: list overflow, boost neigh_modify one" +
                     utils::errorurl(36));
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType> void ComputeUF3Kokkos<DeviceType>::init()
{
  ComputeUF3::init();

  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType, LMPHostType> &&
                           !std::is_same_v<DeviceType, LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType, LMPDeviceType>);
}

/* ---------------------------------------------------------------------- */

template <class DeviceType> void ComputeUF3Kokkos<DeviceType>::compute_array()
{
  atomKK->sync(execution_space, datamask_read);
  atomKK->sync(execution_space, F_MASK);

  if (execution_space != HostKK)
    error->all(FLERR,
               "compute uf3/kk: only Kokkos host execution is implemented; use compute ... uf3/kk/host");

  invoked_array = update->ntimestep;

  for (int irow = 0; irow < size_array_rows; irow++)
    for (int ic = 0; ic < size_array_cols; ic++) uf3local[irow][ic] = 0.0;

  int *save_ilist = nullptr;
  int *save_numneigh = nullptr;
  int **save_firstneigh = nullptr;
  int save_kokkos = 0;

  if (list->kokkos) {
    mirror_kokkos_neighbors_to_cpu();
    save_ilist = list->ilist;
    save_numneigh = list->numneigh;
    save_firstneigh = list->firstneigh;
    save_kokkos = list->kokkos;
    list->ilist = kk_mirror_ilist;
    list->numneigh = kk_mirror_numneigh;
    list->firstneigh = kk_mirror_firstneigh;
    list->kokkos = 0;
  }

  compute_array_core();

  if (save_kokkos) {
    list->ilist = save_ilist;
    list->numneigh = save_numneigh;
    list->firstneigh = save_firstneigh;
    list->kokkos = save_kokkos;
  }
}

namespace LAMMPS_NS {
template class ComputeUF3Kokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ComputeUF3Kokkos<LMPHostType>;
#endif
}    // namespace LAMMPS_NS
