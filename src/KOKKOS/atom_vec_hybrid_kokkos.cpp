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
#include "fix.h"
#include "memory_kokkos.h"
#include "modify.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecHybridKokkos::AtomVecHybridKokkos(LAMMPS *lmp) : AtomVec(lmp),
AtomVecKokkos(lmp), AtomVecHybrid(lmp)
{
  no_comm_vel_flag = 1;
  no_border_vel_flag = 1;
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::grow(int n)
{
  for (int k = 0; k < nstyles; k++) styles[k]->grow(n);
  nmax = atomKK->k_x.h_view.extent(0);

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
    (dynamic_cast<AtomVecKokkos*>(styles[k]))->sort_kokkos(Sorter);
}

/* ---------------------------------------------------------------------- */

int AtomVecHybridKokkos::pack_comm_kokkos(const int &/*n*/, const DAT::tdual_int_2d &/*k_sendlist*/,
                                          const int & /*iswap*/,
                                          const DAT::tdual_xfloat_2d &/*buf*/,
                                          const int &/*pbc_flag*/, const int pbc[])
{
  error->all(FLERR,"AtomVecHybridKokkos doesn't yet support threaded comm");
  return 0;
}

void AtomVecHybridKokkos::unpack_comm_kokkos(const int &/*n*/, const int &/*nfirst*/,
                                             const DAT::tdual_xfloat_2d &/*buf*/)
{
  error->all(FLERR,"AtomVecHybridKokkos doesn't yet support threaded comm");
}

int AtomVecHybridKokkos::pack_comm_self(const int &/*n*/, const DAT::tdual_int_2d &/*list*/,
                                        const int & /*iswap*/, const int /*nfirst*/,
                                        const int &/*pbc_flag*/, const int pbc[])
{
  error->all(FLERR,"AtomVecHybridKokkos doesn't yet support threaded comm");
  return 0;
}

int AtomVecHybridKokkos::pack_border_kokkos(int /*n*/, DAT::tdual_int_2d /*k_sendlist*/,
                                            DAT::tdual_xfloat_2d /*buf*/,int /*iswap*/,
                                            int /*pbc_flag*/, int * /*pbc*/, ExecutionSpace /*space*/)
{
  error->all(FLERR,"AtomVecHybridKokkos doesn't yet support threaded comm");
  return 0;
}

void AtomVecHybridKokkos::unpack_border_kokkos(const int &/*n*/, const int &/*nfirst*/,
                                               const DAT::tdual_xfloat_2d &/*buf*/,
                                               ExecutionSpace /*space*/)
{
  error->all(FLERR,"AtomVecHybridKokkos doesn't yet support threaded comm");
}

int AtomVecHybridKokkos::pack_exchange_kokkos(const int &/*nsend*/,DAT::tdual_xfloat_2d &/*buf*/,
                                              DAT::tdual_int_1d /*k_sendlist*/,
                                              DAT::tdual_int_1d /*k_copylist*/,
                                              ExecutionSpace /*space*/)
{
  error->all(FLERR,"AtomVecHybridKokkos doesn't yet support threaded comm");
  return 0;
}

int AtomVecHybridKokkos::unpack_exchange_kokkos(DAT::tdual_xfloat_2d & /*k_buf*/, int /*nrecv*/,
                                                int /*nlocal*/, int /*dim*/, X_FLOAT /*lo*/,
                                                X_FLOAT /*hi*/, ExecutionSpace /*space*/,
                                                DAT::tdual_int_1d &k_indices)
{
  error->all(FLERR,"AtomVecHybridKokkos doesn't yet support threaded comm");
  return 0;
}

// TODO: move dynamic_cast into init

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::sync(ExecutionSpace space, unsigned int h_mask)
{
  for (int k = 0; k < nstyles; k++) (dynamic_cast<AtomVecKokkos*>(styles[k]))->sync(space,h_mask);
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::sync_overlapping_device(ExecutionSpace space, unsigned int h_mask)
{
  for (int k = 0; k < nstyles; k++) (dynamic_cast<AtomVecKokkos*>(styles[k]))->sync_overlapping_device(space,h_mask);
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::modified(ExecutionSpace space, unsigned int h_mask)
{
  for (int k = 0; k < nstyles; k++) (dynamic_cast<AtomVecKokkos*>(styles[k]))->modified(space,h_mask);
}
