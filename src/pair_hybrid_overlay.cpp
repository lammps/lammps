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

#include "pair_hybrid_overlay.h"
#include <cstring>
#include <cctype>
#include "atom.h"
#include "force.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairHybridOverlay::PairHybridOverlay(LAMMPS *lmp) : PairHybrid(lmp) {}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairHybridOverlay::coeff(int narg, char **arg)
{
  if (narg < 3) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  // 3rd arg = pair sub-style name
  // 4th arg = pair sub-style index if name used multiple times
  // allow for "none" as valid sub-style name

  int multflag;
  int m;

  for (m = 0; m < nstyles; m++) {
    multflag = 0;
    if (strcmp(arg[2],keywords[m]) == 0) {
      if (multiple[m]) {
        multflag = 1;
        if (narg < 4) error->all(FLERR,"Incorrect args for pair coefficients");
        if (!isdigit(arg[3][0]))
          error->all(FLERR,"Incorrect args for pair coefficients");
        int index = force->inumeric(FLERR,arg[3]);
        if (index == multiple[m]) break;
        else continue;
      } else break;
    }
  }

  int none = 0;
  if (m == nstyles) {
    if (strcmp(arg[2],"none") == 0) none = 1;
    else error->all(FLERR,"Pair coeff for hybrid has invalid style");
  }

  // move 1st/2nd args to 2nd/3rd args
  // if multflag: move 1st/2nd args to 3rd/4th args
  // just copy ptrs, since arg[] points into original input line

  arg[2+multflag] = arg[1];
  arg[1+multflag] = arg[0];

  // invoke sub-style coeff() starting with 1st remaining arg

  if (!none) styles[m]->coeff(narg-1-multflag,&arg[1+multflag]);

  // set setflag and which type pairs map to which sub-style
  // if sub-style is none: set hybrid subflag, wipe out map
  // else: set hybrid setflag & map only if substyle setflag is set
  //       if sub-style is new for type pair, add as multiple mapping
  //       if sub-style exists for type pair, don't add, just update coeffs

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      if (none) {
        setflag[i][j] = 1;
        nmap[i][j] = 0;
        count++;
      } else if (styles[m]->setflag[i][j]) {
        int k;
        for (k = 0; k < nmap[i][j]; k++)
          if (map[i][j][k] == m) break;
        if (k == nmap[i][j]) map[i][j][nmap[i][j]++] = m;
        setflag[i][j] = 1;
        count++;
      }
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}
