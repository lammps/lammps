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

#include "nb_standard.h"
#include "atom.h"
#include "group.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NeighBinStandard::NeighBinStandard(LAMMPS *lmp) : NeighBin(lmp) {}

/* ----------------------------------------------------------------------
   bin owned and ghost atoms
------------------------------------------------------------------------- */

void NeighBinStandard::bin_atoms()
{
  int i,ibin;

  for (i = 0; i < mbins; i++) binhead[i] = -1;

  // bin in reverse order so linked list will be in forward order
  // also puts ghost atoms at end of list, which is necessary

  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  if (includegroup) {
    int bitmask = group->bitmask[includegroup];
    for (i = nall-1; i >= nlocal; i--) {
      if (mask[i] & bitmask) {
        ibin = coord2bin(x[i]);
        bins[i] = binhead[ibin];
        binhead[ibin] = i;
      }
    }
    for (i = atom->nfirst-1; i >= 0; i--) {
      ibin = coord2bin(x[i]);
      bins[i] = binhead[ibin];
      binhead[ibin] = i;
    }

  } else {
    for (i = nall-1; i >= 0; i--) {
      ibin = coord2bin(x[i]);
      bins[i] = binhead[ibin];
      binhead[ibin] = i;
    }
  }
}
