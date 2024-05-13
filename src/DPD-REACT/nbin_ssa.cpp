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
   Contributing authors:
   James Larentzos (ARL) and Timothy I. Mattox (Engility Corporation)
------------------------------------------------------------------------- */

#include "nbin_ssa.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "domain.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NBinSSA::NBinSSA(LAMMPS *lmp) : NBinStandard(lmp)
{
  for (int i = 0; i < 8; i++) {
    gairhead_ssa[i] = -1;
  }
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
  int xbin,ybin,zbin;

  last_bin = update->ntimestep;

  bboxlo_[0] = bboxlo[0]; bboxlo_[1] = bboxlo[1]; bboxlo_[2] = bboxlo[2];
  bboxhi_[0] = bboxhi[0]; bboxhi_[1] = bboxhi[1]; bboxhi_[2] = bboxhi[2];

  for (i = 0; i < 8; i++) {
    gairhead_ssa[i] = -1;
  }

  for (i = 0; i < mbins; i++) {
    binhead[i] = -1;
  }

  // bin in reverse order so linked list will be in forward order

  if (includegroup) {
    int bitmask = group->bitmask[includegroup];
    int nowned = atom->nlocal; // NOTE: nlocal was set to atom->nfirst above
    for (i = nall-1; i >= nowned; i--) {
      ibin = coord2ssaAIR(x[i]);
      if (ibin < 1) continue; // skip ghost atoms not in AIR
      if (mask[i] & bitmask) {
        bins[i] = gairhead_ssa[ibin];
        gairhead_ssa[ibin] = i;
      }
    }
  } else {
    for (i = nall-1; i >= nlocal; i--) {
      ibin = coord2ssaAIR(x[i]);
      if (ibin < 1) continue; // skip ghost atoms not in AIR
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
  }

}

/* ---------------------------------------------------------------------- */

void NBinSSA::bin_atoms_setup(int nall)
{
  NBinStandard::bin_atoms_setup(nall); // Setup the parent class's data too

  // Clear the local bin extent bounding box.
  lbinxlo = mbinx - 1; // Safe to = stencil->sx + 1
  lbinylo = mbiny - 1; // Safe to = stencil->sy + 1
  lbinzlo = mbinz - 1; // Safe to = stencil->sz + 1
  lbinxhi = 0; // Safe to = mbinx - stencil->sx - 1
  lbinyhi = 0; // Safe to = mbiny - stencil->sy - 1
  lbinzhi = 0; // Safe to = mbinz - stencil->sz - 1
}

/* ---------------------------------------------------------------------- */

double NBinSSA::memory_usage()
{
  double bytes = NBinStandard::memory_usage(); // Count the parent's usage too

  return bytes;
}

/* ----------------------------------------------------------------------
   convert atom coords into the ssa active interaction region number
------------------------------------------------------------------------- */
int NBinSSA::coord2ssaAIR(const double *x)
{
  int ix, iy, iz;

  ix = iy = iz = 0;
  if (x[2] < domain->sublo[2]) iz = -1;
  if (x[2] >= domain->subhi[2]) iz = 1;
  if (x[1] < domain->sublo[1]) iy = -1;
  if (x[1] >= domain->subhi[1]) iy = 1;
  if (x[0] < domain->sublo[0]) ix = -1;
  if (x[0] >= domain->subhi[0]) ix = 1;

  if (iz < 0) {
    return -1;
  } else if (iz == 0) {
    if (iy<0) return -1; // bottom left/middle/right
    if ((iy==0) && (ix<0) ) return -1; // left atoms
    if ((iy==0) && (ix==0)) return 0; // Locally owned atoms
    if ((iy==0) && (ix>0) ) return 2; // Right atoms
    if ((iy>0)  && (ix==0)) return 1; // Top-middle atoms
    if ((iy>0)  && (ix!=0)) return 3; // Top-right and top-left atoms
  } else { // iz > 0
    if ((ix==0) && (iy==0)) return 4; // Back atoms
    if ((ix==0) && (iy!=0)) return 5; // Top-back and bottom-back atoms
    if ((ix!=0) && (iy==0)) return 6; // Left-back and right-back atoms
    if ((ix!=0) && (iy!=0)) return 7; // Back corner atoms
  }

  return -2;
}
