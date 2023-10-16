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

#include "compute_reaxff_bonds_local_kokkos.h"
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
ComputeReaxFFBondsLocalKokkos<DeviceType>::ComputeReaxFFBondsLocalKokkos(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  alocal(nullptr), reaxff(nullptr)
{
  if (atom->tag_consecutive() == 0)
    error->all(FLERR,"Atom IDs must be consecutive for compute reaxff/bonds/local");

  local_flag = 1;

  nvalues = 7 + 2*MAXREAXBOND;
  prev_nvalues = 0;

  // initialize output

  nlocal = atom->nlocal;

  size_local_rows = atom->nlocal;
  size_local_cols = 7 + 2*MAXREAXBOND;
  printf("RUNNING KOKKOS VERSION\n");
}


/* ---------------------------------------------------------------------- */

template<class DeviceType>
ComputeReaxFFBondsLocalKokkos<DeviceType>::~ComputeReaxFFBondsLocalKokkos()
{
  memoryKK->destroy_kokkos(k_alocal, alocal);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeReaxFFBondsLocalKokkos<DeviceType>::init()
{
  reaxff = dynamic_cast<PairReaxFF*>(force->pair_match("^reax../kk",0));
  if (reaxff == nullptr) error->all(FLERR,"Cannot use fix reaxff/bonds without "
                  "pair_style reaxff/kk");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void ComputeReaxFFBondsLocalKokkos<DeviceType>::compute_local()
{
  invoked_local = update->ntimestep;

  int maxnumbonds = 0;
  if (reaxff->execution_space == Device)
    device_pair()->FindBond(maxnumbonds);
  else
    host_pair()->FindBond(maxnumbonds);
  nvalues = 7+2*maxnumbonds;

  if(atom->nlocal > nlocal || nvalues > prev_nvalues) {
    nlocal = atom->nlocal;
    memoryKK->destroy_kokkos(k_alocal, alocal);
    memoryKK->create_kokkos(k_alocal, alocal, atom->nlocal, nvalues,"reaxff/bonds/local:alocal");
    prev_nvalues = nvalues;
    array_local = alocal;
  }

  size_local_rows = nlocal;
  size_local_cols = nvalues;

  if (reaxff->execution_space == Device)
    device_pair()->PackBondInfo(k_alocal);
  else
    host_pair()->PackBondInfo(k_alocal);
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

template<class DeviceType>
double ComputeReaxFFBondsLocalKokkos<DeviceType>::memory_usage()
{
  double bytes = (double)nlocal*nvalues * sizeof(double);
  return bytes;
}

namespace LAMMPS_NS {
template class ComputeReaxFFBondsLocalKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class ComputeReaxFFBondsLocalKokkos<LMPHostType>;
#endif
}
