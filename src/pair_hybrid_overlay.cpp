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

#include "pair_hybrid_overlay.h"

#include "atom.h"
#include "error.h"

#include <cstring>

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
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  // 3rd arg = pair sub-style name
  // 4th arg = pair sub-style index if name used multiple times
  // allow for "none" as valid sub-style name

  int multflag = 0;
  int m;

  for (m = 0; m < nstyles; m++) {
    multflag = 0;
    if (strcmp(arg[2],keywords[m]) == 0) {
      if (multiple[m]) {
        multflag = 1;
        if (narg < 4) error->all(FLERR,"Incorrect args for pair coefficients");
        if (multiple[m] == utils::inumeric(FLERR,arg[3],false,lmp)) break;
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

  // ensure that one_coeff flag is honored

  if (!none && styles[m]->one_coeff)
    if ((strcmp(arg[0],"*") != 0) || (strcmp(arg[1],"*") != 0))
      error->all(FLERR,"Incorrect args for pair coefficients");

  // invoke sub-style coeff() starting with 1st remaining arg

  if (!none) styles[m]->coeff(narg-1-multflag,arg+1+multflag);

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


/* ----------------------------------------------------------------------
   we need to handle Pair::svector special for hybrid/overlay
------------------------------------------------------------------------- */

void PairHybridOverlay::init_svector()
{
  // single_extra = list all sub-style single_extra
  // allocate svector

  single_extra = 0;
  for (int m = 0; m < nstyles; m++)
    single_extra += styles[m]->single_extra;

  if (single_extra) {
    delete [] svector;
    svector = new double[single_extra];
  }
}

/* ----------------------------------------------------------------------
   we need to handle Pair::svector special for hybrid/overlay
------------------------------------------------------------------------- */

void PairHybridOverlay::copy_svector(int itype, int jtype)
{
  int n=0;
  Pair *this_style = nullptr;

  // fill svector array.
  // copy data from active styles and use 0.0 for inactive ones
  for (int m = 0; m < nstyles; m++) {
    for (int k = 0; k < nmap[itype][jtype]; ++k) {
      if (m == map[itype][jtype][k]) {
        this_style = styles[m];
      } else {
        this_style = nullptr;
      }
    }
    for (int l = 0; l < styles[m]->single_extra; ++l) {
      if (this_style) {
        svector[n++] = this_style->svector[l];
      } else {
        svector[n++] = 0.0;
      }
    }
  }
}
