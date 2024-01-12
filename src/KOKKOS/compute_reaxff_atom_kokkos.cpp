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
   Contributing author: Richard Berger (LANL)
------------------------------------------------------------------------- */

#include "compute_reaxff_atom_kokkos.h"
#include "atom.h"
#include "molecule.h"
#include "update.h"
#include "force.h"
#include "memory.h"
#include "error.h"
#include "neigh_list.h"

#include "memory_kokkos.h"
#include "pair_reaxff_kokkos.h"
#include "reaxff_api.h"

using namespace LAMMPS_NS;
using namespace ReaxFF;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeReaxFFAtomKokkos<DeviceType>::ComputeReaxFFAtomKokkos(LAMMPS *lmp, int narg, char **arg) :
  ComputeReaxFFAtom(lmp, narg, arg),
  nbuf(-1), buf(nullptr)
{
  kokkosable = 1;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeReaxFFAtomKokkos<DeviceType>::~ComputeReaxFFAtomKokkos()
{
  memoryKK->destroy_kokkos(k_buf, buf);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeReaxFFAtomKokkos<DeviceType>::init()
{
  ComputeReaxFFAtom::init();

  if (!reaxff || !reaxff->kokkosable) {
    error->all(FLERR,"Cannot use compute reaxff/atom/kk without "
                     "pair_style reaxff/kk");
  }
}

/* ---------------------------------------------------------------------- */
template<class DeviceType>
void ComputeReaxFFAtomKokkos<DeviceType>::compute_bonds()
{
  if (atom->nlocal > nlocal) {
    memory->destroy(array_atom);
    nlocal = atom->nlocal;
    memory->create(array_atom, nlocal, 3, "reaxff/atom:array_atom");
  }

  // retrieve bond information from kokkos pair style. the data potentially
  // lives on device. it is copied into buf on the host in a condensed format
  // compute_local and compute_atom then expand the data from this buffer into
  // appropiate arrays for consumption by others (e.g. dump local, dump custom
  // or library interface)

  int maxnumbonds = 0;
  if (reaxff->execution_space == Device)
    device_pair()->FindBond(maxnumbonds, groupbit);
  else
    host_pair()->FindBond(maxnumbonds, groupbit);

  nbuf = ((store_bonds ? maxnumbonds*2 : 0) + 3)*nlocal;

  if (!buf || ((int)k_buf.extent(0) < nbuf)) {
    memoryKK->destroy_kokkos(k_buf, buf);
    memoryKK->create_kokkos(k_buf, buf, nbuf, "reaxff/atom:buf");
  }

  // Pass information to buffer, will sync to host

  int nbuf_local;
  if (reaxff->execution_space == Device)
    device_pair()->PackReducedBondBuffer(k_buf, nbuf_local, store_bonds);
  else
    host_pair()->PackReducedBondBuffer(k_buf, nbuf_local, store_bonds);

  // Extract number of bonds from buffer

  nbonds = 0;
  int j = 0;
  for (int i = 0; i < nlocal; i++) {
    int numbonds = static_cast<int>(buf[j+2]);
    nbonds += numbonds;
    j += (store_bonds ? 2*numbonds : 0) + 3;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeReaxFFAtomKokkos<DeviceType>::compute_local()
{
  invoked_local = update->ntimestep;

  if (invoked_bonds < update->ntimestep)
    compute_bonds();

  if (nbonds > prev_nbonds) {
    // grow array_local
    memory->destroy(array_local);
    memory->create(array_local, nbonds, 3, "reaxff/atom:array_local");
    prev_nbonds = nbonds;
  }

  size_local_rows = nbonds;

  // extract local bond information from buffer

  int b = 0;
  int j = 0;
  auto tag = atom->tag;

  for (int i = 0; i < nlocal; ++i) {
    const int numbonds = static_cast<int>(buf[j+2]);
    const int neigh_offset = j + 3;
    const int bo_offset = neigh_offset + numbonds;
    for (int k = 0; k < numbonds; k++) {
      auto bond = array_local[b++];
      bond[0] = tag[i];
      bond[1] = static_cast<tagint> (buf[neigh_offset+k]);
      bond[2] = buf[bo_offset+k];
    }
    j += 2*numbonds + 3;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeReaxFFAtomKokkos<DeviceType>::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  if (invoked_bonds < update->ntimestep)
    compute_bonds();

  // extract peratom bond information from buffer

  int j = 0;
  for (int i = 0; i < nlocal; ++i) {
    auto ptr = array_atom[i];
    int numbonds = static_cast<int>(buf[j+2]);
    ptr[0] = buf[j]; // sbo
    ptr[1] = buf[j+1]; // nlp
    ptr[2] = numbonds;
    j += (store_bonds ? 2*numbonds : 0) + 3;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

template<class DeviceType>
double ComputeReaxFFAtomKokkos<DeviceType>::memory_usage()
{
  double bytes = (double)(nlocal*3) * sizeof(double);
  if (store_bonds)
    bytes += (double)(nbonds*3) * sizeof(double);
  bytes += (double)(nbuf > 0 ? nbuf * sizeof(double) : 0);
  return bytes;
}

namespace LAMMPS_NS {
template class ComputeReaxFFAtomKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ComputeReaxFFAtomKokkos<LMPHostType>;
#endif
}
