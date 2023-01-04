// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Trimright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "npair_trim_kokkos.h"
#include "atom_kokkos.h"
#include "atom_masks.h"
#include "neigh_list_kokkos.h"
#include "my_page.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
NPairTrimKokkos<DeviceType>::NPairTrimKokkos(LAMMPS *lmp) : NPair(lmp) {}

/* ----------------------------------------------------------------------
   create list which is simply a copy of parent list
------------------------------------------------------------------------- */

template<class DeviceType>
void NPairTrimKokkos<DeviceType>::build(NeighList *list)
{
  NeighList *listcopy = list->listcopy;

  cutsq_custom = cutoff_custom*cutoff_custom;

  if (list->kokkos) {
    if (!listcopy->kokkos)
      error->all(FLERR,"Cannot trim non-Kokkos neighbor list to Kokkos neighbor list");
    trim_to_kokkos(list);
  } else {
    if (!listcopy->kokkos)
      error->all(FLERR,"Missing Kokkos neighbor list for trim");
    trim_to_cpu(list);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void NPairTrimKokkos<DeviceType>::trim_to_kokkos(NeighList *list)
{
  x = atomKK->k_x.view<DeviceType>();
  atomKK->sync(execution_space,X_MASK);

  cutsq_custom = cutoff_custom*cutoff_custom;

  NeighListKokkos<DeviceType>* k_list_copy = static_cast<NeighListKokkos<DeviceType>*>(list->listcopy);
  d_ilist_copy = k_list_copy->d_ilist;
  d_numneigh_copy = k_list_copy->d_numneigh;
  d_neighbors_copy = k_list_copy->d_neighbors;
  int inum_copy = list->listcopy->inum;
  if (list->ghost) inum_copy += list->listcopy->gnum;

  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  k_list->maxneighs = k_list_copy->maxneighs; // simple, but could be made more memory efficient
  k_list->grow(atom->nmax);
  d_ilist = k_list->d_ilist;
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;

  // loop over parent list and trim

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagNPairTrim>(0,inum_copy),*this);
  copymode = 0;

  list->inum = k_list_copy->inum;
  list->gnum = k_list_copy->gnum;

  k_list->k_ilist.template modify<DeviceType>();
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void NPairTrimKokkos<DeviceType>::operator()(TagNPairTrim, const int &ii) const {
  int n = 0;

  const int i = d_ilist_copy(ii);
  const double xtmp = x(i,0);
  const double ytmp = x(i,1);
  const double ztmp = x(i,2);

  // loop over copy neighbor list

  const int jnum = d_numneigh_copy(i);
  const AtomNeighbors neighbors_i = AtomNeighbors(&d_neighbors(i,0),d_numneigh(i),
                                                  &d_neighbors(i,1)-&d_neighbors(i,0));

  for (int jj = 0; jj < jnum; jj++) {
    const int joriginal = d_neighbors_copy(i,jj);
    const int j = joriginal & NEIGHMASK;

    const double delx = xtmp - x(j,0);
    const double dely = ytmp - x(j,1);
    const double delz = ztmp - x(j,2);
    const double rsq = delx*delx + dely*dely + delz*delz;

    if (rsq > cutsq_custom) continue;

    neighbors_i(n++) = joriginal;
  }

  d_numneigh(i) = n;
  d_ilist(ii) = i;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void NPairTrimKokkos<DeviceType>::trim_to_cpu(NeighList *list)
{
  NeighList *listcopy = list->listcopy;
  NeighListKokkos<DeviceType>* listcopy_kk = (NeighListKokkos<DeviceType>*) listcopy;

  listcopy_kk->k_ilist.template sync<LMPHostType>();

  double** x = atom->x;

  int inum = listcopy->inum;
  int gnum = listcopy->gnum;
  int inum_all = inum;
  if (list->ghost) inum_all += gnum;
  auto h_ilist = listcopy_kk->k_ilist.h_view;
  auto h_numneigh = Kokkos::create_mirror_view_and_copy(LMPHostType(),listcopy_kk->d_numneigh);
  auto h_neighbors = Kokkos::create_mirror_view_and_copy(LMPHostType(),listcopy_kk->d_neighbors);

  list->inum = inum;
  list->gnum = gnum;
  auto ilist = list->ilist;
  auto numneigh = list->numneigh;

  // Kokkos neighbor data is stored differently than regular CPU,
  //  must loop over lists

  int *neighptr;
  int **firstneigh = list->firstneigh;
  MyPage<int> *ipage = list->ipage;
  ipage->reset();

  for (int ii = 0; ii < inum_all; ii++) {
    int n = 0;
    neighptr = ipage->vget();

    const int i = h_ilist[ii];
    ilist[ii] = i;
    const double xtmp = x[i][0];
    const double ytmp = x[i][1];
    const double ztmp = x[i][2];

    // loop over Kokkos neighbor list

    const int jnum = h_numneigh[i];

    for (int jj = 0; jj < jnum; jj++) {
      const int joriginal = h_neighbors(i,jj);

      const int j = joriginal & NEIGHMASK;

      const double delx = xtmp - x[j][0];
      const double dely = ytmp - x[j][1];
      const double delz = ztmp - x[j][2];
      const double rsq = delx*delx + dely*dely + delz*delz;

      if (rsq > cutsq_custom) continue;

      neighptr[n++] = joriginal;
    }

    firstneigh[i] = neighptr;
    numneigh[i] = n;
    ipage->vgot(n);
    if (ipage->status())
      error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
  }
}

namespace LAMMPS_NS {
template class NPairTrimKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class NPairTrimKokkos<LMPHostType>;
#endif
}
