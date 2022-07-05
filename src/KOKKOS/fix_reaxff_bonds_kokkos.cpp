// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Stan Moore (Sandia)
------------------------------------------------------------------------- */

#include "fix_reaxff_bonds_kokkos.h"

#include "atom.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "memory_kokkos.h"
#include "pair_reaxff_kokkos.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixReaxFFBondsKokkos::FixReaxFFBondsKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixReaxFFBonds(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;

  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

void FixReaxFFBondsKokkos::init()
{
  Pair *pair_kk = force->pair_match("^reax../kk",0);
  if (pair_kk == nullptr) error->all(FLERR,"Cannot use fix reaxff/bonds without "
                  "pair_style reaxff/kk");

  FixReaxFFBonds::init();
}

/* ---------------------------------------------------------------------- */

void FixReaxFFBondsKokkos::Output_ReaxFF_Bonds()

{
  int nbuf_local;
  int nlocal_max, numbonds, numbonds_max;
  double *buf;
  DAT::tdual_ffloat_1d k_buf;

  int nlocal = atom->nlocal;
  int nlocal_tot = static_cast<int> (atom->natoms);

  numbonds = 0;
  if (reaxff->execution_space == Device)
    ((PairReaxFFKokkos<LMPDeviceType>*) reaxff)->FindBond(numbonds);
  else
    ((PairReaxFFKokkos<LMPHostType>*) reaxff)->FindBond(numbonds);

  // allocate a temporary buffer for the snapshot info
  MPI_Allreduce(&numbonds,&numbonds_max,1,MPI_INT,MPI_MAX,world);
  MPI_Allreduce(&nlocal,&nlocal_max,1,MPI_INT,MPI_MAX,world);

  nbuf = 1+(numbonds_max*2+10)*nlocal_max;
  memoryKK->create_kokkos(k_buf,buf,nbuf,"reaxff/bonds:buf");

  // Pass information to buffer
  if (reaxff->execution_space == Device)
    ((PairReaxFFKokkos<LMPDeviceType>*) reaxff)->PackBondBuffer(k_buf,nbuf_local);
  else
    ((PairReaxFFKokkos<LMPHostType>*) reaxff)->PackBondBuffer(k_buf,nbuf_local);
  buf[0] = nlocal;

  // Receive information from buffer for output
  RecvBuffer(buf, nbuf, nbuf_local, nlocal_tot, numbonds_max);

  memoryKK->destroy_kokkos(k_buf,buf);
}

/* ---------------------------------------------------------------------- */

double FixReaxFFBondsKokkos::memory_usage()
{
  double bytes;

  bytes = nbuf*sizeof(double);
  // These are accounted for in PairReaxFFKokkos:
  //bytes += nmax*sizeof(int);
  //bytes += 1.0*nmax*MAXREAXBOND*sizeof(double);
  //bytes += 1.0*nmax*MAXREAXBOND*sizeof(int);

  return bytes;
}
