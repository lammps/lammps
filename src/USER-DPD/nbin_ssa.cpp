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
   Contributing authors:
   James Larentzos (ARL) and Timothy I. Mattox (Engility Corporation)
------------------------------------------------------------------------- */

#include "nbin_ssa.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NBinSSA::NBinSSA(LAMMPS *lmp) : NBinStandard(lmp)
{
  maxbin_ssa = 0;
  bins_ssa = NULL;
  maxhead_ssa = 0;
  binhead_ssa = NULL;
  for (int i = 0; i < 9; i++) gairhead_ssa[i] = -1;
}

NBinSSA::~NBinSSA()
{
  memory->destroy(bins_ssa);
  memory->destroy(binhead_ssa);
}

/* ----------------------------------------------------------------------
   bin owned and ghost atoms for the Shardlow Splitting Algorithm (SSA)
   local atoms are in distinct bins (binhead_ssa) from the ghosts
   ghost atoms are in distinct bins (gbinhead_ssa) from the locals
     ghosts which are not in an Active Interaction Region (AIR) are skipped
------------------------------------------------------------------------- */

void NBinSSA::bin_atoms()
{
  int i,ibin;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  if (includegroup) nlocal = atom->nfirst;
  double **x = atom->x;
  int *mask = atom->mask;
  int *ssaAIR = atom->ssaAIR;

  last_bin = update->ntimestep;

  bboxlo_[0] = bboxlo[0]; bboxlo_[1] = bboxlo[1]; bboxlo_[2] = bboxlo[2];
  bboxhi_[0] = bboxhi[0]; bboxhi_[1] = bboxhi[1]; bboxhi_[2] = bboxhi[2];

  for (i = 0; i < 9; i++) {
    gairhead_ssa[i] = -1;
  }

  for (i = 0; i < mbins; i++) {
    binhead_ssa[i] = -1;
  }

  // bin in reverse order so linked list will be in forward order

  if (includegroup) {
    int bitmask = group->bitmask[includegroup];
    int nowned = atom->nlocal; // NOTE: nlocal was set to atom->nfirst above
    for (i = nall-1; i >= nowned; i--) {
      ibin = ssaAIR[i];
      if (ibin < 2) continue; // skip ghost atoms not in AIR
      if (mask[i] & bitmask) {
        bins_ssa[i] = gairhead_ssa[ibin];
        gairhead_ssa[ibin] = i;
      }
    }
  } else {
    for (i = nall-1; i >= nlocal; i--) {
      ibin = ssaAIR[i];
      if (ibin < 2) continue; // skip ghost atoms not in AIR
      bins_ssa[i] = gairhead_ssa[ibin];
      gairhead_ssa[ibin] = i;
    }
  }
  for (i = nlocal-1; i >= 0; i--) {
    ibin = coord2bin(x[i][0], x[i][1], x[i][2]);
    bins_ssa[i] = binhead_ssa[ibin];
    binhead_ssa[ibin] = i;
  }
}

/* ---------------------------------------------------------------------- */

void NBinSSA::bin_atoms_setup(int nall)
{
  NBinStandard::bin_atoms_setup(nall); // Setup the parent class's data too

  if (mbins > maxhead_ssa) {
    maxhead_ssa = mbins;
    memory->destroy(binhead_ssa);
    memory->create(binhead_ssa,maxhead_ssa,"binhead_ssa");
  }

  if (nall > maxbin_ssa) {
    maxbin_ssa = nall;
    memory->destroy(bins_ssa);
    memory->create(bins_ssa,maxbin_ssa,"bins_ssa");
  }
}

/* ---------------------------------------------------------------------- */

bigint NBinSSA::memory_usage()
{
  bigint bytes = NBinStandard::memory_usage(); // Count the parent's usage too

  if (maxbin_ssa) bytes += memory->usage(bins_ssa,maxbin_ssa);
  if (maxhead_ssa) {
    bytes += memory->usage(binhead_ssa,maxhead_ssa);
  }
  return bytes;
}
