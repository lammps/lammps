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

#include "string.h"
#include "pair_hybrid_overlay.h"
#include "atom.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

PairHybridOverlay::PairHybridOverlay(LAMMPS *lmp) : PairHybrid(lmp) {}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairHybridOverlay::coeff(int narg, char **arg)
{
  if (narg < 3) error->all("Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  // 3rd arg = pair sub-style name
  // allow for "none" as valid sub-style name

  int m;
  for (m = 0; m < nstyles; m++)
    if (strcmp(arg[2],keywords[m]) == 0) break;

  int none = 0;
  if (m == nstyles) {
    if (strcmp(arg[2],"none") == 0) none = 1;
    else error->all("Pair coeff for hybrid has invalid style");
  }

  // move 1st/2nd args to 2nd/3rd args
  // just copy ptrs, since arg[] points into original input line

  arg[2] = arg[1];
  arg[1] = arg[0];

  // invoke sub-style coeff() starting with 1st arg

  if (!none) styles[m]->coeff(narg-1,&arg[1]);

  // set setflag and which type pairs map to which sub-style
  // if sub-style is none: set hybrid subflag, wipe out map
  // else: set hybrid setflag & map only if substyle setflag is set
  //       multiple mappings are allowed,
  //       but don't map again if already mapped to this sub-style

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      if (none) {
	setflag[i][j] = 1;
	nmap[i][j] = 0;
	count++;
      } else if (styles[m]->setflag[i][j] ) {
	int k;
	for (k = 0; k < nmap[i][j]; k++)
	  if (map[i][j][k] == m) break;
	if (k < nmap[i][j]) continue;
	setflag[i][j] = 1;
	map[i][j][nmap[i][j]++] = m;
	count++;
      }
    }
  }

  if (count == 0) error->all("Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   combine sub-style neigh list requests and create new ones if needed
------------------------------------------------------------------------- */

void PairHybridOverlay::modify_requests()
{
  int i,j;
  NeighRequest *irq,*jrq;

  // loop over pair requests only
  // if a previous list is same kind with same skip attributes
  // then make this one a copy list of that one
  // works whether both lists are no-skip or yes-skip
  // will not point a list at a copy list, but at copy list's parent

  for (i = 0; i < neighbor->nrequest; i++) {
    if (!neighbor->requests[i]->pair) continue;

    irq = neighbor->requests[i];
    for (j = 0; j < i; j++) {
      if (!neighbor->requests[j]->pair) continue;
      jrq = neighbor->requests[j];
      if (irq->same_kind(jrq) && irq->same_skip(jrq)) {
	irq->copy = 1;
	irq->otherlist = j;
	break;
      }
    }
  }

  // perform same operations on skip lists as pair style = hybrid

  PairHybrid::modify_requests();
}
