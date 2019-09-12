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

#include "fix_property_atom_kokkos.h"
#include <cstdlib>
#include <cstring>
#include "atom_kokkos.h"
#include "comm.h"
#include "memory_kokkos.h"
#include "error.h"
#include "update.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{MOLECULE,CHARGE,RMASS,INTEGER,DOUBLE};

/* ---------------------------------------------------------------------- */

FixPropertyAtomKokkos::FixPropertyAtomKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixPropertyAtom(lmp, narg, arg)
{
  atomKK = (AtomKokkos *) atom;
  grow_arrays(atom->nmax);
}

/* ----------------------------------------------------------------------
   allocate atom-based arrays
   initialize new values to 0,
   since AtomVec class won't do it as atoms are added,
   e.g. in create_atom() or data_atom()
------------------------------------------------------------------------- */

void FixPropertyAtomKokkos::grow_arrays(int nmax)
{
  for (int m = 0; m < nvalue; m++) {
    if (style[m] == MOLECULE) {
      memory->grow(atom->molecule,nmax,"atom:molecule");
      size_t nbytes = (nmax-nmax_old) * sizeof(tagint);
      memset(&atom->molecule[nmax_old],0,nbytes);
    } else if (style[m] == CHARGE) {
      memory->grow(atom->q,nmax,"atom:q");
      size_t nbytes = (nmax-nmax_old) * sizeof(double);
      memset(&atom->q[nmax_old],0,nbytes);
    } else if (style[m] == RMASS) {
      memory->grow(atom->rmass,nmax,"atom:rmass");
      size_t nbytes = (nmax-nmax_old) * sizeof(double);
      memset(&atom->rmass[nmax_old],0,nbytes);
    } else if (style[m] == INTEGER) {
      memory->grow(atom->ivector[index[m]],nmax,"atom:ivector");
      size_t nbytes = (nmax-nmax_old) * sizeof(int);
      memset(&atom->ivector[index[m]][nmax_old],0,nbytes);
    } else if (style[m] == DOUBLE) {
      atomKK->sync(Device,DVECTOR_MASK);
      memoryKK->grow_kokkos(atomKK->k_dvector,atomKK->dvector,atomKK->k_dvector.extent(0),nmax,
                          "atom:dvector");
      atomKK->modified(Device,DVECTOR_MASK);
      //memory->grow(atom->dvector[index[m]],nmax,"atom:dvector");
      //size_t nbytes = (nmax-nmax_old) * sizeof(double);
      //memset(&atom->dvector[index[m]][nmax_old],0,nbytes);
    }
  }

  nmax_old = nmax;
}
