/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
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

#include <cstdlib>
#include <cstring>
#include "fix_ave_atom.h"
#include "fix_reaxc_bonds_kokkos.h"
#include "atom.h"
#include "update.h"
#include "pair_reaxc_kokkos.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "comm.h"
#include "force.h"
#include "compute.h"
#include "input.h"
#include "variable.h"
#include "memory_kokkos.h"
#include "error.h"
#include "reaxc_list.h"
#include "reaxc_types.h"
#include "reaxc_defs.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixReaxCBondsKokkos::FixReaxCBondsKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixReaxCBonds(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;

  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
}

/* ---------------------------------------------------------------------- */

FixReaxCBondsKokkos::~FixReaxCBondsKokkos()
{

}

/* ---------------------------------------------------------------------- */

void FixReaxCBondsKokkos::init()
{
  Pair *pair_kk = force->pair_match("reax/c/kk",0);
  if (pair_kk == NULL) error->all(FLERR,"Cannot use fix reax/c/bonds without "
                  "pair_style reax/c/kk");

  FixReaxCBonds::init();
}

/* ---------------------------------------------------------------------- */

void FixReaxCBondsKokkos::Output_ReaxC_Bonds(bigint ntimestep, FILE *fp)

{
  int nbuf_local;
  int nlocal_max, numbonds, numbonds_max;
  double *buf;
  DAT::tdual_ffloat_1d k_buf;

  int nlocal = atom->nlocal;
  int nlocal_tot = static_cast<int> (atom->natoms);

  numbonds = 0;
  if (reaxc->execution_space == Device)
    ((PairReaxCKokkos<LMPDeviceType>*) reaxc)->FindBond(numbonds);
  else
    ((PairReaxCKokkos<LMPHostType>*) reaxc)->FindBond(numbonds);

  // allocate a temporary buffer for the snapshot info
  MPI_Allreduce(&numbonds,&numbonds_max,1,MPI_INT,MPI_MAX,world);
  MPI_Allreduce(&nlocal,&nlocal_max,1,MPI_INT,MPI_MAX,world);

  nbuf = 1+(numbonds_max*2+10)*nlocal_max;
  memoryKK->create_kokkos(k_buf,buf,nbuf,"reax/c/bonds:buf");

  // Pass information to buffer
  if (reaxc->execution_space == Device)
    ((PairReaxCKokkos<LMPDeviceType>*) reaxc)->PackBondBuffer(k_buf,nbuf_local);
  else
    ((PairReaxCKokkos<LMPHostType>*) reaxc)->PackBondBuffer(k_buf,nbuf_local);
  buf[0] = nlocal;

  // Receive information from buffer for output
  RecvBuffer(buf, nbuf, nbuf_local, nlocal_tot, numbonds_max);

  memoryKK->destroy_kokkos(k_buf,buf);
}

/* ---------------------------------------------------------------------- */

double FixReaxCBondsKokkos::memory_usage()
{
  double bytes;

  bytes = nbuf*sizeof(double);
  // These are accounted for in PairReaxCKokkos:
  //bytes += nmax*sizeof(int);
  //bytes += 1.0*nmax*MAXREAXBOND*sizeof(double);
  //bytes += 1.0*nmax*MAXREAXBOND*sizeof(int);

  return bytes;
}
