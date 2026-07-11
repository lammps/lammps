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

#include "fix_wall_srd.h"

#include "domain.h"
#include "error.h"
#include "fix.h"
#include "fix_wall.h"
#include "graphics.h"
#include "input.h"
#include "lattice.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;

enum { XLO, XHI, YLO, YHI, ZLO, ZHI };
enum { NONE, EDGE, CONSTANT, VARIABLE };

/* ---------------------------------------------------------------------- */

FixWallSRD::FixWallSRD(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), nwall(0), fwall(nullptr), fwall_all(nullptr),
    imgobjs(nullptr), imgparms(nullptr)
{
  if (narg < 4) utils::missing_cmd_args(FLERR, "fix wall/srd", error);

  // parse args

  nwall = 0;
  int scaleflag = 1;

  int iarg = 3;
  while (iarg < narg) {
    const std::string thisarg(arg[iarg]);
    if ((thisarg == "xlo") || (thisarg == "ylo") || (thisarg == "zlo")
        || (thisarg == "xhi") || (thisarg == "yhi") || (thisarg == "zhi")) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix wall/srd " + thisarg, error);

      int newwall;
      if (thisarg == "xlo") newwall = XLO;
      else if (thisarg == "xhi") newwall = XHI;
      else if (thisarg == "ylo") newwall = YLO;
      else if (thisarg == "yhi") newwall = YHI;
      else if (thisarg == "zlo") newwall = ZLO;
      else if (thisarg == "zhi") newwall = ZHI;

      for (int m = 0; (m < nwall) && (m < 6); m++)
        if (newwall == wallwhich[m])
          error->all(FLERR, iarg, "Wall {} defined twice in fix wall/srd command", thisarg);

      wallwhich[nwall] = newwall;
      if (strcmp(arg[iarg+1], "EDGE") == 0) {
        wallstyle[nwall] = EDGE;
        int dim = wallwhich[nwall] / 2;
        int side = wallwhich[nwall] % 2;
        if (side == 0)
          coord0[nwall] = domain->boxlo[dim];
        else
          coord0[nwall] = domain->boxhi[dim];
      } else if (utils::strmatch(arg[iarg+1], "^v_")) {
        wallstyle[nwall] = VARIABLE;
        varstr[nwall] = utils::strdup(arg[iarg+1] + 2);
      } else {
        wallstyle[nwall] = CONSTANT;
        coord0[nwall] = utils::numeric(FLERR, arg[iarg+1], false, lmp);
      }

      nwall++;
      iarg += 2;

    } else if (thisarg == "units") {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "fix wall/srd units", error);
      if (strcmp(arg[iarg+1], "box") == 0)
        scaleflag = 0;
      else if (strcmp(arg[iarg+1], "lattice") == 0)
        scaleflag = 1;
      else
        error->all(FLERR, iarg+1, "Unknown fix wall/srd units setting {}", arg[iarg+1]);
      iarg += 2;
    } else
      error->all(FLERR, iarg, "Unknown fix wall/srd keyword {}", arg[iarg]);
  }

  // error check

  if (nwall == 0) error->all(FLERR, "Illegal fix wall command");

  for (int m = 0; m < nwall; m++) {
    if ((wallwhich[m] == XLO || wallwhich[m] == XHI) && domain->xperiodic)
      error->all(FLERR, "Cannot use fix wall/srd xlo or xhi with periodic x dimension");
    if ((wallwhich[m] == YLO || wallwhich[m] == YHI) && domain->yperiodic)
      error->all(FLERR, "Cannot use fix wall/srd ylo or yhi with periodic y dimension");
    if ((wallwhich[m] == ZLO || wallwhich[m] == ZHI) && domain->zperiodic)
      error->all(FLERR, "Cannot use fix wall/srd zlo or zhi with periodic z dimension");
  }

  for (int m = 0; m < nwall; m++)
    if ((wallwhich[m] == ZLO || wallwhich[m] == ZHI) && domain->dimension == 2)
      error->all(FLERR, "Cannot use fix wall/srd zlo/zhi for a 2d simulation");

  // setup wall force array

  array_flag = 1;
  size_array_rows = nwall;
  size_array_cols = 3;
  global_freq = 1;
  extarray = 1;

  memory->create(fwall, nwall, 3, "wall/srd:fwall");
  memory->create(fwall_all, nwall, 3, "wall/srd:fwall_all");

  // scale coord for CONSTANT walls

  int flag = 0;
  for (int m = 0; m < nwall; m++)
    if (wallstyle[m] == CONSTANT) flag = 1;

  if (flag) {
    double xscale, yscale, zscale;
    if (scaleflag) {
      xscale = domain->lattice->xlattice;
      yscale = domain->lattice->ylattice;
      zscale = domain->lattice->zlattice;
    } else
      xscale = yscale = zscale = 1.0;

    double scale;
    for (int m = 0; m < nwall; m++) {
      if (wallwhich[m] < YLO)
        scale = xscale;
      else if (wallwhich[m] < ZLO)
        scale = yscale;
      else
        scale = zscale;
      if (wallstyle[m] == CONSTANT) coord0[m] *= scale;
    }
  }

  // set overlap if walls exist in multiple dimensions

  int dimflag[3];
  dimflag[0] = dimflag[1] = dimflag[2] = 0;
  for (int m = 0; m < nwall; m++) dimflag[wallwhich[m] / 2] = 1;
  if (dimflag[0] + dimflag[1] + dimflag[2] > 1)
    overlap = 1;
  else
    overlap = 0;

  // set varflag if any wall positions are variable

  varflag = 0;
  for (int m = 0; m < nwall; m++)
    if (wallstyle[m] == VARIABLE) varflag = 1;
  laststep = -1;

  // for rendering walls with dump image.
  if (domain->dimension == 2) {
    // one cylinder object per wall to draw in 2d
    memory->create(imgobjs, nwall, "fix_wall:imgobjs");
    memory->create(imgparms, nwall, 8, "fix_wall:imgparms");
    for (int m = 0; m < nwall; ++m) {
      imgobjs[m] = Graphics::CYLINDER;
      imgparms[m][0] = 1;    // use color of first atom type by default
    }
  } else {
    // two triangle objects per wall to draw in 3d
    memory->create(imgobjs, 2 * nwall, "fix_wall:imgobjs");
    memory->create(imgparms, 2 * nwall, 10, "fix_wall:imgparms");
    for (int m = 0; m < nwall; ++m) {
      imgobjs[2 * m] = Graphics::TRIANGLE;
      imgobjs[2 * m + 1] = Graphics::TRIANGLE;
      imgparms[2 * m][0] = 1;        // use color of first atom type by default
      imgparms[2 * m + 1][0] = 1;    // use color of first atom type by default
    }
  }
}

/* ---------------------------------------------------------------------- */

FixWallSRD::~FixWallSRD()
{
  for (int m = 0; m < nwall; m++)
    if (wallstyle[m] == VARIABLE) delete[] varstr[m];
  memory->destroy(fwall);
  memory->destroy(fwall_all);

  memory->destroy(imgobjs);
  memory->destroy(imgparms);
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
  if (modify->get_fix_by_style("^srd").size() == 0)
    error->all(FLERR, Error::NOLASTLINE, "Cannot use fix wall/srd without fix srd");

  for (int m = 0; m < nwall; m++) {
    if (wallstyle[m] != VARIABLE) continue;
    varindex[m] = input->variable->find(varstr[m]);
    if (varindex[m] < 0)
      error->all(FLERR, Error::NOLASTLINE,
                 "Variable {} for fix wall/srd does not exist", varstr[m]);
    if (!input->variable->equalstyle(varindex[m]))
      error->all(FLERR, Error::NOLASTLINE,
                 "Variable {} for fix wall/srd is invalid style", varstr[m]);
  }

  dt = update->dt;
}

/* ----------------------------------------------------------------------
   return force component on a wall
------------------------------------------------------------------------- */

double FixWallSRD::compute_array(int i, int j)
{
  // only sum across procs one time

  if (!force_flag) {
    MPI_Allreduce(&fwall[0][0], &fwall_all[0][0], 3 * nwall, MPI_DOUBLE, MPI_SUM, world);
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
    else
      xnew = coord0[m];

    if (laststep < 0) {
      xwall[m] = xwalllast[m] = xnew;
      vwall[m] = 0.0;
    } else if (laststep < ntimestep) {
      xwalllast[m] = xwall[m];
      xwall[m] = xnew;
      vwall[m] = (xwall[m] - xwalllast[m]) / dt;
    }

    fwall[m][0] = fwall[m][1] = fwall[m][2] = 0.0;

    FixWall::update_image_plane(m, wallwhich[m], xnew, imgparms, domain);
  }

  laststep = ntimestep;

  if (varflag) modify->addstep_compute(update->ntimestep + 1);

  if (flag)
    for (int m = 0; m < nwall; m++) xwallhold[m] = xwall[m];

  force_flag = 0;
}

/* ----------------------------------------------------------------------
   provide graphics information to dump image to render wall as plane
   data has been copied to dedicated storage during fix indent execution
------------------------------------------------------------------------- */

int FixWallSRD::image(int *&objs, double **&parms)
{
  objs = imgobjs;
  parms = imgparms;
  if (domain->dimension == 2) {
    return nwall;
  } else {
    return 2 * nwall;
  }
}
