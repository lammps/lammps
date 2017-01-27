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
  binlist_ssa = NULL;
  binct_ssa = NULL;
  for (int i = 0; i < 9; i++) {
    gairhead_ssa[i] = -1;
  }
}

NBinSSA::~NBinSSA()
{
  memory->destroy(binlist_ssa);
  memory->destroy(binct_ssa);
}

/* ----------------------------------------------------------------------
   bin owned and ghost atoms for the Shardlow Splitting Algorithm (SSA)
   local atoms are in distinct bins (binhead[]) from the ghosts
   ghost atoms are "binned" in gairhead_ssa[] instead
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
  int xbin,ybin,zbin;

  last_bin = update->ntimestep;

  bboxlo_[0] = bboxlo[0]; bboxlo_[1] = bboxlo[1]; bboxlo_[2] = bboxlo[2];
  bboxhi_[0] = bboxhi[0]; bboxhi_[1] = bboxhi[1]; bboxhi_[2] = bboxhi[2];

  for (i = 0; i < 9; i++) {
    gairhead_ssa[i] = -1;
  }

  for (i = 0; i < mbins; i++) {
    binhead[i] = -1;
    binlist_ssa[i] = -1;
    binct_ssa[i] = 0;
  }

  // bin in reverse order so linked list will be in forward order

  if (includegroup) {
    int bitmask = group->bitmask[includegroup];
    int nowned = atom->nlocal; // NOTE: nlocal was set to atom->nfirst above
    for (i = nall-1; i >= nowned; i--) {
      ibin = ssaAIR[i];
      if (ibin < 2) continue; // skip ghost atoms not in AIR
      if (mask[i] & bitmask) {
        bins[i] = gairhead_ssa[ibin];
        gairhead_ssa[ibin] = i;
      }
    }
  } else {
    for (i = nall-1; i >= nlocal; i--) {
      ibin = ssaAIR[i];
      if (ibin < 2) continue; // skip ghost atoms not in AIR
      bins[i] = gairhead_ssa[ibin];
      gairhead_ssa[ibin] = i;
    }
  }
  for (i = nlocal-1; i >= 0; i--) {
    ibin = coord2bin(x[i][0], x[i][1], x[i][2], xbin, ybin, zbin);
    // Find the bounding box of the local atoms in the bins
    if (xbin < lbinxlo) lbinxlo = xbin;
    if (xbin >= lbinxhi) lbinxhi = xbin + 1;
    if (ybin < lbinylo) lbinylo = ybin;
    if (ybin >= lbinyhi) lbinyhi = ybin + 1;
    if (zbin < lbinzlo) lbinzlo = zbin;
    if (zbin >= lbinzhi) lbinzhi = zbin + 1;
    bins[i] = binhead[ibin];
    binhead[ibin] = i;
    ++(binct_ssa[ibin]);
  }

}

/* ---------------------------------------------------------------------- */

void NBinSSA::bin_atoms_setup(int nall)
{
  NBinStandard::bin_atoms_setup(nall); // Setup the parent class's data too

  if (mbins > maxbin_ssa) {
    maxbin_ssa = mbins;
    memory->destroy(binlist_ssa);
    memory->destroy(binct_ssa);
    memory->create(binlist_ssa,maxbin_ssa,"binlist_ssa");
    memory->create(binct_ssa,maxbin_ssa,"binct_ssa");
  }

  // Clear the local bin extent bounding box.
  lbinxlo = mbinx - 1; // Safe to = stencil->sx + 1
  lbinylo = mbiny - 1; // Safe to = stencil->sy + 1
  lbinzlo = mbinz - 1; // Safe to = stencil->sz + 1
  lbinxhi = 0; // Safe to = mbinx - stencil->sx - 1
  lbinyhi = 0; // Safe to = mbiny - stencil->sy - 1
  lbinzhi = 0; // Safe to = mbinz - stencil->sz - 1
}

/* ---------------------------------------------------------------------- */

bigint NBinSSA::memory_usage()
{
  bigint bytes = NBinStandard::memory_usage(); // Count the parent's usage too

  if (maxbin_ssa) {
    bytes += memory->usage(binlist_ssa,maxbin_ssa);
    bytes += memory->usage(binct_ssa,maxbin_ssa);
  }
  return bytes;
}
