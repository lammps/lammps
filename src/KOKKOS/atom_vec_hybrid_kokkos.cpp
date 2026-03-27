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

#include "atom_vec_hybrid_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "domain.h"
#include "error.h"
#include "kokkos.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecHybridKokkos::AtomVecHybridKokkos(LAMMPS *lmp) : AtomVec(lmp),
AtomVecKokkos(lmp), AtomVecHybrid(lmp)
{
}

/* ----------------------------------------------------------------------
   process field strings to initialize data structs for all other methods
------------------------------------------------------------------------- */

void AtomVecHybridKokkos::init()
{
  AtomVecHybrid::init();

  set_atom_masks();
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::process_args(int narg, char **arg)
{
  AtomVecHybrid::process_args(narg,arg);

  // Seems process_args is first to be called, so set up stylesKK here
  stylesKK = new AtomVecKokkos*[nstyles];
  for (int k = 0; k < nstyles; k++) {
    stylesKK[k] = dynamic_cast<AtomVecKokkos*>(styles[k]);
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::grow(int n)
{
  for (int k = 0; k < nstyles; k++) styles[k]->grow(n);
  nmax = atomKK->k_x.view_host().extent(0);

  // ptrs for atom_vec.cpp
  tag = atom->tag;
  type = atom->type;
  mask = atom->mask;
  image = atom->image;
  x = atom->x;
  v = atom->v;
  f = atom->f;
}

/* ----------------------------------------------------------------------
   sort atom arrays on device
------------------------------------------------------------------------- */

void AtomVecHybridKokkos::sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter)
{
  for (int k = 0; k < nstyles; k++)
    stylesKK[k]->sort_kokkos(Sorter);
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::pack_comm_bonus_kokkos(const int &n, const DAT::tdual_int_1d &list,
                                                 const DAT::tdual_double_2d_lr &buf, int vel_flag)
{
  for (int k = 0; k < nstyles; k++)
    stylesKK[k]->pack_comm_bonus_kokkos(n,list,buf,vel_flag);
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::unpack_comm_bonus_kokkos(const int &n, const int &nfirst,
                                                   const DAT::tdual_double_2d_lr &buf, int vel_flag)
{
  for (int k = 0; k < nstyles; k++)
    stylesKK[k]->unpack_comm_bonus_kokkos(n,nfirst,buf,vel_flag);
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::pack_comm_self_bonus_kokkos(const int &n, const DAT::tdual_int_1d &list,
                                                      const int nfirst)
{
  for (int k = 0; k < nstyles; k++)
    stylesKK[k]->pack_comm_self_bonus_kokkos(n,list,nfirst);
}

/* ---------------------------------------------------------------------- */
void AtomVecHybridKokkos::pack_comm_self_fused_bonus_kokkos(const int &n,
                                         const DAT::tdual_int_2d_lr &list,
                                         const DAT::tdual_int_1d &sendnum_scan,
                                         const DAT::tdual_int_1d &firstrecv,
                                         const DAT::tdual_int_1d &g2l)
{
  for (int k = 0; k < nstyles; k++)
    stylesKK[k]->pack_comm_self_fused_bonus_kokkos(n,list,sendnum_scan,firstrecv,g2l);
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::pack_border_bonus_kokkos(int n, DAT::tdual_int_1d k_sendlist,
                                                   DAT::tdual_double_2d_lr &buf,
                                                   ExecutionSpace space, int vel_flag)
{
  for (int k = 0; k < nstyles; k++)
    stylesKK[k]->pack_border_bonus_kokkos(n,k_sendlist,buf,space,vel_flag);
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::unpack_border_bonus_kokkos(const int &n, const int &nfirst,
                                                     const DAT::tdual_double_2d_lr &buf,
                                                     ExecutionSpace space, int vel_flag)
{
  for (int k = 0; k < nstyles; k++)
    stylesKK[k]->unpack_border_bonus_kokkos(n,nfirst,buf,space,vel_flag);
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::pack_exchange_bonus_kokkos(const int &nsend, DAT::tdual_double_2d_lr &buf,
                                                     DAT::tdual_int_1d k_sendlist,
                                                     DAT::tdual_int_1d k_copylist,
                                                     DAT::tdual_int_1d k_copylist_bonus,
                                                     ExecutionSpace space)
{
  for (int k = 0; k < nstyles; k++)
    stylesKK[k]->pack_exchange_bonus_kokkos(nsend,buf,k_sendlist,k_copylist,
                                 k_copylist_bonus,space);
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::unpack_exchange_bonus_kokkos(DAT::tdual_double_2d_lr &k_buf,
                                                       int nrecv,
                                                       ExecutionSpace space,
                                                       DAT::tdual_int_1d &k_indices)
{
  for (int k = 0; k < nstyles; k++)
    stylesKK[k]->unpack_exchange_bonus_kokkos(k_buf,nrecv,space,k_indices);
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::set_size_exchange()
{
  AtomVecKokkos::set_size_exchange();

  size_exchange_bonus = 0;
  for (int k = 0; k < nstyles; k++)
    size_exchange_bonus += stylesKK[k]->size_exchange_bonus;

  size_exchange += size_exchange_bonus;

  for (int k = 0; k < nstyles; k++)
    stylesKK[k]->size_exchange = size_exchange;
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::sync(ExecutionSpace space, uint64_t h_mask)
{
  for (int k = 0; k < nstyles; k++) stylesKK[k]->sync(space,h_mask);
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::sync_pinned(ExecutionSpace space, uint64_t h_mask, int async_flag)
{
  for (int k = 0; k < nstyles; k++) stylesKK[k]->sync_pinned(space,h_mask,async_flag);
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::modified(ExecutionSpace space, uint64_t h_mask)
{
  for (int k = 0; k < nstyles; k++) stylesKK[k]->modified(space,h_mask);
}
