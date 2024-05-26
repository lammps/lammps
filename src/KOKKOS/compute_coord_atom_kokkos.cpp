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

#include "compute_coord_atom_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "compute_orientorder_atom_kokkos.h"
#include "error.h"
#include "memory_kokkos.h"
#include "neigh_request.h"
#include "neighbor_kokkos.h"
#include "update.h"

using namespace LAMMPS_NS;


/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeCoordAtomKokkos<DeviceType>::ComputeCoordAtomKokkos(LAMMPS *lmp, int narg, char **arg) :
  ComputeCoordAtom(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  d_typelo = typename AT::t_int_1d("coord/atom:typelo",ncol);
  d_typehi = typename AT::t_int_1d("coord/atom:typehi",ncol);

  auto h_typelo = Kokkos::create_mirror_view(d_typelo);
  auto h_typehi = Kokkos::create_mirror_view(d_typehi);

  for (int i = 0; i < ncol; i++) {
    h_typelo(i) = typelo[i];
    h_typehi(i) = typehi[i];
  }

  Kokkos::deep_copy(d_typelo,h_typelo);
  Kokkos::deep_copy(d_typehi,h_typehi);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeCoordAtomKokkos<DeviceType>::~ComputeCoordAtomKokkos()
{
  if (copymode) return;

  memoryKK->destroy_kokkos(k_cvec,cvec);
  memoryKK->destroy_kokkos(k_carray,carray);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeCoordAtomKokkos<DeviceType>::init()
{
  ComputeCoordAtom::init();

  // adjust neighbor list request for KOKKOS

  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeCoordAtomKokkos<DeviceType>::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow coordination array if necessary

  if (atom->nmax > nmax) {
    if (ncol == 1) {
      memoryKK->destroy_kokkos(k_cvec,cvec);
      nmax = atom->nmax;
      memoryKK->create_kokkos(k_cvec,cvec,nmax,"coord/atom:cvec");
      vector_atom = cvec;
      d_cvec = k_cvec.template view<DeviceType>();
    } else {
      memoryKK->destroy_kokkos(k_carray,carray);
      nmax = atom->nmax;
      memoryKK->create_kokkos(k_carray,carray,nmax,ncol,"coord/atom:carray");
      array_atom = carray;
      d_carray = k_carray.template view<DeviceType>();
    }
  }

  if (cstyle == ORIENT) {
    if (!(c_orientorder->invoked_flag & Compute::INVOKED_PERATOM)) {
      c_orientorder->compute_peratom();
      c_orientorder->invoked_flag |= Compute::INVOKED_PERATOM;
    }
    nqlist = c_orientorder->nqlist;
    normv = c_orientorder->array_atom;
    comm->forward_comm(this);

    if (!c_orientorder->kokkosable)
      error->all(FLERR,"Must use compute orientorder/atom/kk with compute coord/atom/kk");

    if (c_orientorder->execution_space == Host) {
      ComputeOrientOrderAtomKokkos<LMPHostType>* c_orientorder_kk;
      c_orientorder_kk = (ComputeOrientOrderAtomKokkos<LMPHostType>*) c_orientorder;
      c_orientorder_kk->k_qnarray.modify<LMPHostType>();
      c_orientorder_kk->k_qnarray.sync<DeviceType>();
      d_normv = c_orientorder_kk->k_qnarray.view<DeviceType>();
    } else {
      ComputeOrientOrderAtomKokkos<LMPDeviceType>* c_orientorder_kk;
      c_orientorder_kk = (ComputeOrientOrderAtomKokkos<LMPDeviceType>*) c_orientorder;
      c_orientorder_kk->k_qnarray.modify<LMPHostType>();
      c_orientorder_kk->k_qnarray.sync<DeviceType>();
      d_normv = c_orientorder_kk->k_qnarray.view<DeviceType>();
    }
  }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  inum = list->inum;
  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;

  // compute coordination number(s) for each atom in group
  // use full neighbor list to count atoms less than cutoff

  atomKK->sync(execution_space,X_MASK|TYPE_MASK|MASK_MASK);
  x = atomKK->k_x.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();

  copymode = 1;
  if (cstyle == CUTOFF) {
    if (ncol == 1) {
      typename Kokkos::RangePolicy<DeviceType, TagComputeCoordAtom<CUTOFF,1> > policy(0,inum);
      Kokkos::parallel_for("ComputeCoordAtom",policy,*this);
    } else {
      typename Kokkos::RangePolicy<DeviceType, TagComputeCoordAtom<CUTOFF,0> > policy(0,inum);
      Kokkos::parallel_for("ComputeCoordAtom",policy,*this);
    }
  } else if (cstyle == ORIENT) {
    typename Kokkos::RangePolicy<DeviceType, TagComputeCoordAtom<ORIENT,1> > policy(0,inum);
    Kokkos::parallel_for("ComputeCoordAtom",policy,*this);
  }
  copymode = 0;

  if (ncol == 1 || cstyle == ORIENT) {
    k_cvec.modify<DeviceType>();
    k_cvec.sync<LMPHostType>();
  } else {
    k_carray.modify<DeviceType>();
    k_carray.sync<LMPHostType>();
  }

}

template<class DeviceType>
template<int CSTYLE, int NCOL>
KOKKOS_INLINE_FUNCTION
void ComputeCoordAtomKokkos<DeviceType>::operator()(TagComputeCoordAtom<CSTYLE,NCOL>, const int &ii) const
{
  const int i = d_ilist[ii];
  if (NCOL == 1 || CSTYLE == ORIENT)
    d_cvec(i) = 0.0;
  else
    for (int m = 0; m < ncol; m++) d_carray(i,m) = 0.0;
  if (mask[i] & groupbit) {
    const X_FLOAT xtmp = x(i,0);
    const X_FLOAT ytmp = x(i,1);
    const X_FLOAT ztmp = x(i,2);
    const int jnum = d_numneigh[i];

    int n = 0;
    for (int jj = 0; jj < jnum; jj++) {
      int j = d_neighbors(i,jj);
      j &= NEIGHMASK;

      if (NCOL == 1)
        if (!(mask[j] & jgroupbit)) continue;

      const int jtype = type[j];
      const F_FLOAT delx = x(j,0) - xtmp;
      const F_FLOAT dely = x(j,1) - ytmp;
      const F_FLOAT delz = x(j,2) - ztmp;
      const F_FLOAT rsq = delx*delx + dely*dely + delz*delz;
      if (rsq < cutsq) {
        if (CSTYLE == CUTOFF) {
          if (NCOL == 1) {
            if (jtype >= d_typelo[0] && jtype <= d_typehi[0])
              n++;
          } else {
              for (int m = 0; m < ncol; m++)
                if (jtype >= d_typelo[m] && jtype <= d_typehi[m])
                  d_carray(i,m) += 1.0;
          }
        } else if (CSTYLE == ORIENT) {
            double dot_product = 0.0;
            for (int m=0; m < 2*(2*l+1); m++) {
              dot_product += d_normv(i,nqlist+m)*d_normv(j,nqlist+m);
            }
            if (dot_product > threshold) n++;
        }
      }
    }

    if (NCOL == 1 || CSTYLE == ORIENT)
      d_cvec[i] = n;
  }

}

namespace LAMMPS_NS {
template class ComputeCoordAtomKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ComputeCoordAtomKokkos<LMPHostType>;
#endif
}
