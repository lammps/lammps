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

#include "compute_reaxff_bonds_kokkos.h"
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
ComputeReaxFFBondsKokkos<DeviceType>::ComputeReaxFFBondsKokkos(LAMMPS *lmp, int narg, char **arg) :
  ComputeReaxFFBonds(lmp, narg, arg),
  nbuf(-1), buf(nullptr)
{
  kokkosable = 1;
}


/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeReaxFFBondsKokkos<DeviceType>::~ComputeReaxFFBondsKokkos()
{
  memoryKK->destroy_kokkos(k_buf, buf);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeReaxFFBondsKokkos<DeviceType>::init()
{
  reaxff = dynamic_cast<PairReaxFF*>(force->pair_match("^reax../kk",0));
  if (reaxff == nullptr) error->all(FLERR,"Cannot use compute reaxff/bonds without "
                  "pair_style reaxff/kk");
}

/* ---------------------------------------------------------------------- */
template<class DeviceType>
void ComputeReaxFFBondsKokkos<DeviceType>::compute_bonds()
{
  if (atom->nlocal > nlocal) {
    memory->destroy(array_atom);
    nlocal = atom->nlocal;
    memory->create(array_atom, nlocal, 3, "reaxff/bonds:array_atom");
  }

  // retrieve bond information from kokkos pair style. the data potentially
  // lives on device. it is copied into buf on the host in a condensed format
  // compute_local and compute_atom then expand the data from this buffer into
  // appropiate arrays for consumption by others (e.g. dump local, dump custom
  // or library interface)

  int maxnumbonds = 0;
  if (reaxff->execution_space == Device)
    device_pair()->FindBond(maxnumbonds);
  else
    host_pair()->FindBond(maxnumbonds);

  nbuf = (maxnumbonds*2 + 3)*nlocal;

  if(!buf || k_buf.extent(0) < nbuf) {
    memoryKK->destroy_kokkos(k_buf, buf);
    memoryKK->create_kokkos(k_buf, buf, nbuf, "reaxff/bonds:buf");
  }

  // Pass information to buffer, will sync to host

  int nbuf_local;
  if (reaxff->execution_space == Device)
    device_pair()->PackReducedBondBuffer(k_buf, nbuf_local);
  else
    host_pair()->PackReducedBondBuffer(k_buf, nbuf_local);

  // Extract number of bonds from buffer

  nbonds = 0;
  int j = 0;
  for (int i = 0; i < nlocal; i++) {
    int numbonds = static_cast<int>(buf[j+2]);
    nbonds += numbonds;
    j += 2*numbonds + 3;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeReaxFFBondsKokkos<DeviceType>::compute_local()
{
  invoked_local = update->ntimestep;

  if(invoked_bonds < update->ntimestep)
    compute_bonds();

  if(nbonds > prev_nbonds) {
    // grow array_local
    memory->destroy(array_local);
    memory->create(array_local, nbonds, 3, "reaxff/bonds:array_local");
    prev_nbonds = nbonds;
  }

  size_local_rows = nbonds;

  // extract local bond information from buffer

  int b = 0;
  int j = 0;

  for (int i = 0; i < nlocal; ++i) {
    const int numbonds = static_cast<int>(buf[j+2]);
    const int neigh_offset = j + 3;
    const int bo_offset = neigh_offset + numbonds;
    for (int k = 0; k < numbonds; k++) {
      auto bond = array_local[b++];
      bond[0] = i;
      bond[1] = static_cast<tagint> (buf[neigh_offset+k]);
      bond[2] = buf[bo_offset+k];
    }
    j += 2*numbonds + 3;
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeReaxFFBondsKokkos<DeviceType>::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  if(invoked_bonds < update->ntimestep)
    compute_bonds();

  // extract peratom bond information from buffer

  int j = 0;
  for (int i = 0; i < nlocal; ++i) {
    auto ptr = array_atom[i];
    int numbonds = static_cast<int>(buf[j+2]);
    ptr[0] = buf[j]; // sbo
    ptr[1] = buf[j+1]; // nlp
    ptr[2] = numbonds;
    j += 2*numbonds + 3;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

template<class DeviceType>
double ComputeReaxFFBondsKokkos<DeviceType>::memory_usage()
{
  double bytes = (double)(nbonds*3) * sizeof(double);
  bytes += (double)(nlocal*3) * sizeof(double);
  return bytes;
}

namespace LAMMPS_NS {
template class ComputeReaxFFBondsKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ComputeReaxFFBondsKokkos<LMPHostType>;
#endif
}
