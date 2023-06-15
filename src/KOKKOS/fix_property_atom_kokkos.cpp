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

#include "fix_property_atom_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "memory_kokkos.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

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
  for (int nv = 0; nv < nvalue; nv++) {
    if (styles[nv] == MOLECULE) {
      memory->grow(atom->molecule,nmax,"atom:molecule");
      size_t nbytes = (nmax-nmax_old) * sizeof(tagint);
      memset(&atom->molecule[nmax_old],0,nbytes);
    } else if (styles[nv] == CHARGE) {
      memory->grow(atom->q,nmax,"atom:q");
      size_t nbytes = (nmax-nmax_old) * sizeof(double);
      memset(&atom->q[nmax_old],0,nbytes);
    } else if (styles[nv] == RMASS) {
      memory->grow(atom->rmass,nmax,"atom:rmass");
      size_t nbytes = (nmax-nmax_old) * sizeof(double);
      memset(&atom->rmass[nmax_old],0,nbytes);
    } else if (styles[nv] == TEMPERATURE) {
      memory->grow(atom->temperature, nmax, "atom:temperature");
      size_t nbytes = (nmax - nmax_old) * sizeof(double);
      memset(&atom->temperature[nmax_old], 0, nbytes);
    } else if (styles[nv] == HEATFLOW) {
      memory->grow(atom->heatflow, nmax, "atom:heatflow");
      size_t nbytes = (nmax - nmax_old) * sizeof(double);
      memset(&atom->heatflow[nmax_old], 0, nbytes);
    } else if (styles[nv] == IVEC) {
      memory->grow(atom->ivector[index[nv]],nmax,"atom:ivector");
      size_t nbytes = (nmax-nmax_old) * sizeof(int);
      memset(&atom->ivector[index[nv]][nmax_old],0,nbytes);
    } else if (styles[nv] == DVEC) {
      atomKK->sync(Device,DVECTOR_MASK);
      memoryKK->grow_kokkos(atomKK->k_dvector,atomKK->dvector,atomKK->k_dvector.extent(0),nmax,
                          "atom:dvector");
      atomKK->modified(Device,DVECTOR_MASK);
    } else if (styles[nv] == IARRAY) {
      memory->grow(atom->iarray[index[nv]], nmax, cols[nv], "atom:iarray");
      size_t nbytes = (size_t) (nmax - nmax_old) * cols[nv] * sizeof(int);
      if (nbytes) memset(&atom->iarray[index[nv]][nmax_old][0], 0, nbytes);
    } else if (styles[nv] == DARRAY) {
      memory->grow(atom->darray[index[nv]], nmax, cols[nv], "atom:darray");
      size_t nbytes = (size_t) (nmax - nmax_old) * cols[nv] * sizeof(double);
      if (nbytes) memset(&atom->darray[index[nv]][nmax_old][0], 0, nbytes);
    }
  }
  nmax_old = nmax;
}
