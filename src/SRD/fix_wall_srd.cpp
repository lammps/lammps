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

#include "stdlib.h"
#include "string.h"
#include "fix_wall_srd.h"
#include "atom.h"
#include "modify.h"
#include "fix.h"
#include "domain.h"
#include "lattice.h"
#include "input.h"
#include "modify.h"
#include "update.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{XLO,XHI,YLO,YHI,ZLO,ZHI};
enum{NONE,EDGE,CONSTANT,VARIABLE};

/* ---------------------------------------------------------------------- */

FixWallSRD::FixWallSRD(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4) error->all(FLERR,"Illegal fix wall/srd command");

  // parse args

  nwall = 0;
  int scaleflag = 1;

  int iarg = 3;
  while (iarg < narg) {
    if ((strcmp(arg[iarg],"xlo") == 0) || (strcmp(arg[iarg],"xhi") == 0) ||
        (strcmp(arg[iarg],"ylo") == 0) || (strcmp(arg[iarg],"yhi") == 0) ||
        (strcmp(arg[iarg],"zlo") == 0) || (strcmp(arg[iarg],"zhi") == 0)) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix wall/srd command");

      int newwall;
      if (strcmp(arg[iarg],"xlo") == 0) newwall = XLO;
      else if (strcmp(arg[iarg],"xhi") == 0) newwall = XHI;
      else if (strcmp(arg[iarg],"ylo") == 0) newwall = YLO;
      else if (strcmp(arg[iarg],"yhi") == 0) newwall = YHI;
      else if (strcmp(arg[iarg],"zlo") == 0) newwall = ZLO;
      else if (strcmp(arg[iarg],"zhi") == 0) newwall = ZHI;

      for (int m = 0; m < nwall; m++)
        if (newwall == wallwhich[m])
          error->all(FLERR,"Wall defined twice in fix wall/srd command");

      wallwhich[nwall] = newwall;
      if (strcmp(arg[iarg+1],"EDGE") == 0) {
        wallstyle[nwall] = EDGE;
        int dim = wallwhich[nwall] / 2;
        int side = wallwhich[nwall] % 2;
        if (side == 0) coord0[nwall] = domain->boxlo[dim];
        else coord0[nwall] = domain->boxhi[dim];
      } else if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
        wallstyle[nwall] = VARIABLE;
        int n = strlen(&arg[iarg+1][2]) + 1;
        varstr[nwall] = new char[n];
        strcpy(varstr[nwall],&arg[iarg+1][2]);
      } else {
        wallstyle[nwall] = CONSTANT;
        coord0[nwall] = force->numeric(FLERR,arg[iarg+1]);
      }

      nwall++;
      iarg += 2;

    } else if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal wall/srd command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all(FLERR,"Illegal fix wall/srd command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix wall/srd command");
  }

  // error check

  if (nwall == 0) error->all(FLERR,"Illegal fix wall command");

  for (int m = 0; m < nwall; m++) {
    if ((wallwhich[m] == XLO || wallwhich[m] == XHI) && domain->xperiodic)
      error->all(FLERR,"Cannot use fix wall/srd in periodic dimension");
    if ((wallwhich[m] == YLO || wallwhich[m] == YHI) && domain->yperiodic)
      error->all(FLERR,"Cannot use fix wall/srd in periodic dimension");
    if ((wallwhich[m] == ZLO || wallwhich[m] == ZHI) && domain->zperiodic)
      error->all(FLERR,"Cannot use fix wall/srd in periodic dimension");
  }

  for (int m = 0; m < nwall; m++)
    if ((wallwhich[m] == ZLO || wallwhich[m] == ZHI) && domain->dimension == 2)
      error->all(FLERR,"Cannot use fix wall/srd zlo/zhi for a 2d simulation");

  // setup wall force array

  array_flag = 1;
  size_array_rows = nwall;
  size_array_cols = 3;
  global_freq = 1;
  extarray = 1;

  memory->create(fwall,nwall,3,"wall/srd:fwall");
  memory->create(fwall_all,nwall,3,"wall/srd:fwall_all");

  // scale coord for CONSTANT walls

  int flag = 0;
  for (int m = 0; m < nwall; m++)
    if (wallstyle[m] == CONSTANT) flag = 1;

  if (flag) {
    double xscale,yscale,zscale;
    if (scaleflag) {
      xscale = domain->lattice->xlattice;
      yscale = domain->lattice->ylattice;
      zscale = domain->lattice->zlattice;
    }
    else xscale = yscale = zscale = 1.0;

    double scale;
    for (int m = 0; m < nwall; m++) {
      if (wallwhich[m] < YLO) scale = xscale;
      else if (wallwhich[m] < ZLO) scale = yscale;
      else scale = zscale;
      if (wallstyle[m] == CONSTANT) coord0[m] *= scale;
    }
  }

  // set overlap if walls exist in multiple dimensions

  int dimflag[3];
  dimflag[0] = dimflag[1] = dimflag[2] = 0;
  for (int m = 0; m < nwall; m++)
    dimflag[wallwhich[m]/2] = 1;
  if (dimflag[0] + dimflag[1] + dimflag[2] > 1) overlap = 1;
  else overlap = 0;

  // set varflag if any wall positions are variable

  varflag = 0;
  for (int m = 0; m < nwall; m++)
    if (wallstyle[m] == VARIABLE) varflag = 1;
  laststep = -1;
}

/* ---------------------------------------------------------------------- */

FixWallSRD::~FixWallSRD()
{
  for (int m = 0; m < nwall; m++)
    if (wallstyle[m] == VARIABLE) delete [] varstr[m];
  memory->destroy(fwall);
  memory->destroy(fwall_all);
}

/* ---------------------------------------------------------------------- */

int FixWallSRD::setmask()
{
  int mask = 0;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixWallSRD::init()
{
  int flag = 0;
  for (int m = 0; m < modify->nfix; m++)
    if (strcmp(modify->fix[m]->style,"srd") == 0) flag = 1;
  if (!flag) error->all(FLERR,"Cannot use fix wall/srd without fix srd");

  for (int m = 0; m < nwall; m++) {
    if (wallstyle[m] != VARIABLE) continue;
    varindex[m] = input->variable->find(varstr[m]);
    if (varindex[m] < 0)
      error->all(FLERR,"Variable name for fix wall/srd does not exist");
    if (!input->variable->equalstyle(varindex[m]))
      error->all(FLERR,"Variable for fix wall/srd is invalid style");
  }

  dt = update->dt;
}

/* ----------------------------------------------------------------------
   return force component on a wall
------------------------------------------------------------------------- */

double FixWallSRD::compute_array(int i, int j)
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(&fwall[0][0],&fwall_all[0][0],3*nwall,
                  MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return fwall_all[i][j];
}

/* ----------------------------------------------------------------------
   set wall position and velocity, zero forces on walls
   evaluate variable if necessary, wrap with clear/add
   if flag, then being called on reneighbor, so archive wall positions
------------------------------------------------------------------------- */

void FixWallSRD::wall_params(int flag)
{
  double xnew;

  if (varflag) modify->clearstep_compute();

  bigint ntimestep = update->ntimestep;

  for (int m = 0; m < nwall; m++) {
    if (wallstyle[m] == VARIABLE)
      xnew = input->variable->compute_equal(varindex[m]);
    else xnew = coord0[m];

    if (laststep < 0) {
      xwall[m] = xwalllast[m] = xnew;
      vwall[m] = 0.0;
    } else if (laststep < ntimestep) {
      xwalllast[m] = xwall[m];
      xwall[m] = xnew;
      vwall[m] = (xwall[m] - xwalllast[m]) / dt;
    }

    fwall[m][0] = fwall[m][1] = fwall[m][2] = 0.0;
  }

  laststep = ntimestep;

  if (varflag) modify->addstep_compute(update->ntimestep + 1);

  if (flag)
    for (int m = 0; m < nwall; m++)
      xwallhold[m] = xwall[m];

  force_flag = 0;
}
