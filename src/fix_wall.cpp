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

#include "fix_wall.h"

#include "domain.h"
#include "error.h"
#include "input.h"
#include "lattice.h"
#include "modify.h"
#include "respa.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

enum { XLO = 0, XHI = 1, YLO = 2, YHI = 3, ZLO = 4, ZHI = 5 };
enum { NONE = 0, EDGE, CONSTANT, VARIABLE };

static const char *wallpos[] = {"xlo", "xhi", "ylo", "yhi", "zlo", "zhi"};

/* ---------------------------------------------------------------------- */

FixWall::FixWall(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg), nwall(0)
{
  scalar_flag = 1;
  vector_flag = 1;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;
  energy_global_flag = 1;
  virial_global_flag = virial_peratom_flag = 1;
  respa_level_support = 1;
  ilevel_respa = 0;

  // parse args

  int scaleflag = 1;
  fldflag = 0;
  int pbcflag = 0;

  for (int i = 0; i < 6; i++) xstr[i] = estr[i] = sstr[i] = lstr[i] = fstr[i] = kstr[i] = nullptr;

  int iarg = 3;
  if (utils::strmatch(style, "^wall/table")) iarg = 5;

  while (iarg < narg) {
    int wantargs = 5;
    if (utils::strmatch(style, "^wall/lepton")) wantargs = 4;
    if (utils::strmatch(style, "^wall/morse")) wantargs = 6;

    if ((strcmp(arg[iarg], "xlo") == 0) || (strcmp(arg[iarg], "xhi") == 0) ||
        (strcmp(arg[iarg], "ylo") == 0) || (strcmp(arg[iarg], "yhi") == 0) ||
        (strcmp(arg[iarg], "zlo") == 0) || (strcmp(arg[iarg], "zhi") == 0)) {
      if (iarg + wantargs > narg) error->all(FLERR, "Missing argument for fix {} command", style);

      int newwall;
      if (strcmp(arg[iarg], "xlo") == 0) {
        newwall = XLO;
      } else if (strcmp(arg[iarg], "xhi") == 0) {
        newwall = XHI;
      } else if (strcmp(arg[iarg], "ylo") == 0) {
        newwall = YLO;
      } else if (strcmp(arg[iarg], "yhi") == 0) {
        newwall = YHI;
      } else if (strcmp(arg[iarg], "zlo") == 0) {
        newwall = ZLO;
      } else if (strcmp(arg[iarg], "zhi") == 0) {
        newwall = ZHI;
      }
      for (int m = 0; (m < nwall) && (m < 6); m++) {
        if (newwall == wallwhich[m])
          error->all(FLERR, "{} wall defined twice in fix {} command", wallpos[newwall], style);
      }
      wallwhich[nwall] = newwall;

      if (strcmp(arg[iarg + 1], "EDGE") == 0) {
        xstyle[nwall] = EDGE;
        int dim = wallwhich[nwall] / 2;
        int side = wallwhich[nwall] % 2;
        if (side == 0) {
          coord0[nwall] = domain->boxlo[dim];
        } else {
          coord0[nwall] = domain->boxhi[dim];
        }
      } else if (utils::strmatch(arg[iarg + 1], "^v_")) {
        xstyle[nwall] = VARIABLE;
        xstr[nwall] = utils::strdup(arg[iarg + 1] + 2);
      } else {
        xstyle[nwall] = CONSTANT;
        coord0[nwall] = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      }

      if (utils::strmatch(style, "^wall/lepton")) {
        estyle[nwall] = sstyle[nwall] = CONSTANT;
        lstr[nwall] = utils::strdup(arg[iarg + 2]);
        cutoff[nwall] = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      } else if (utils::strmatch(style, "^wall/table")) {
        estyle[nwall] = sstyle[nwall] = CONSTANT;
        fstr[nwall] = utils::strdup(arg[iarg + 2]);
        kstr[nwall] = utils::strdup(arg[iarg + 3]);
        cutoff[nwall] = utils::numeric(FLERR, arg[iarg + 4], false, lmp);
      } else {
        if (iarg + 5 > narg) error->all(FLERR, "Missing argument for fix {} command", style);

        if (utils::strmatch(arg[iarg + 2], "^v_")) {
          estr[nwall] = utils::strdup(arg[iarg + 2] + 2);
          estyle[nwall] = VARIABLE;
        } else {
          epsilon[nwall] = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
          estyle[nwall] = CONSTANT;
        }

        if (utils::strmatch(style, "^wall/morse")) {
          if (utils::strmatch(arg[iarg + 3], "^v_")) {
            astr[nwall] = utils::strdup(arg[iarg + 3] + 2);
            astyle[nwall] = VARIABLE;
          } else {
            alpha[nwall] = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
            astyle[nwall] = CONSTANT;
          }
          // adjust so we can share the regular code path
          ++iarg;
          --wantargs;
        }

        if (utils::strmatch(arg[iarg + 3], "^v_")) {
          sstr[nwall] = utils::strdup(arg[iarg + 3] + 2);
          sstyle[nwall] = VARIABLE;
        } else {
          sigma[nwall] = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
          sstyle[nwall] = CONSTANT;
        }
        cutoff[nwall] = utils::numeric(FLERR, arg[iarg + 4], false, lmp);
      }
      nwall++;
      iarg += wantargs;
    } else if (strcmp(arg[iarg], "units") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix {} command", style);
      if (strcmp(arg[iarg + 1], "box") == 0)
        scaleflag = 0;
      else if (strcmp(arg[iarg + 1], "lattice") == 0)
        scaleflag = 1;
      else
        error->all(FLERR, "Illegal fix {} command", style);
      iarg += 2;
    } else if (strcmp(arg[iarg], "fld") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix {} command", style);
      fldflag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "pbc") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix {} command", style);
      pbcflag = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else
      error->all(FLERR, "Illegal fix {} command", style);
  }

  size_vector = nwall;

  // error checks

  if (nwall == 0) error->all(FLERR, "Illegal fix {} command: no walls defined", style);
  for (int m = 0; m < nwall; m++) {
    if (cutoff[m] <= 0.0)
      error->all(FLERR, "Fix {} cutoff <= 0.0 for {} wall", style, wallpos[wallwhich[m]]);
  }

  for (int m = 0; m < nwall; m++)
    if ((wallwhich[m] == ZLO || wallwhich[m] == ZHI) && domain->dimension == 2)
      error->all(FLERR, "Cannot use fix {} zlo/zhi for a 2d simulation", style);

  if (!pbcflag) {
    for (int m = 0; m < nwall; m++) {
      if ((wallwhich[m] == XLO || wallwhich[m] == XHI) && domain->xperiodic)
        error->all(FLERR, "Cannot use {} wall in periodic x dimension", wallpos[wallwhich[m]]);
      if ((wallwhich[m] == YLO || wallwhich[m] == YHI) && domain->yperiodic)
        error->all(FLERR, "Cannot use {} wall in periodic y dimension", wallpos[wallwhich[m]]);
      if ((wallwhich[m] == ZLO || wallwhich[m] == ZHI) && domain->zperiodic)
        error->all(FLERR, "Cannot use {} wall in periodic z dimension", wallpos[wallwhich[m]]);
    }
  }

  // scale factors for wall position for CONSTANT and VARIABLE walls

  int flag = 0;
  for (int m = 0; m < nwall; m++)
    if (xstyle[m] != EDGE) flag = 1;

  if (flag) {
    if (scaleflag) {
      xscale = domain->lattice->xlattice;
      yscale = domain->lattice->ylattice;
      zscale = domain->lattice->zlattice;
    } else
      xscale = yscale = zscale = 1.0;

    for (int m = 0; m < nwall; m++) {
      if (xstyle[m] != CONSTANT) continue;
      if (wallwhich[m] < YLO)
        coord0[m] *= xscale;
      else if (wallwhich[m] < ZLO)
        coord0[m] *= yscale;
      else
        coord0[m] *= zscale;
    }
  }

  // set xflag if any wall positions are variable
  // set varflag if any wall positions or parameters are variable
  // set wstyle to VARIABLE if either epsilon or sigma is a variable

  varflag = xflag = 0;
  for (int m = 0; m < nwall; m++) {
    if (xstyle[m] == VARIABLE) xflag = 1;
    if (xflag || estyle[m] == VARIABLE || sstyle[m] == VARIABLE) varflag = 1;
    if (estyle[m] == VARIABLE || sstyle[m] == VARIABLE)
      wstyle[m] = VARIABLE;
    else
      wstyle[m] = CONSTANT;
  }

  eflag = 0;
  for (int m = 0; m <= nwall; m++) ewall[m] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixWall::~FixWall()
{
  if (copymode) return;

  for (int m = 0; m < nwall; m++) {
    delete[] xstr[m];
    delete[] estr[m];
    delete[] sstr[m];
    delete[] lstr[m];
    delete[] fstr[m];
    delete[] kstr[m];
  }
}

/* ---------------------------------------------------------------------- */

int FixWall::setmask()
{
  int mask = 0;

  // FLD implicit needs to invoke wall forces before pair style

  if (fldflag)
    mask |= PRE_FORCE;
  else
    mask |= POST_FORCE;

  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixWall::init()
{
  for (int m = 0; m < nwall; m++) {
    if (xstyle[m] == VARIABLE) {
      xindex[m] = input->variable->find(xstr[m]);
      if (xindex[m] < 0) error->all(FLERR, "Variable name for fix wall does not exist");
      if (!input->variable->equalstyle(xindex[m]))
        error->all(FLERR, "Variable for fix wall is invalid style");
    }
    if (estyle[m] == VARIABLE) {
      eindex[m] = input->variable->find(estr[m]);
      if (eindex[m] < 0) error->all(FLERR, "Variable name for fix wall does not exist");
      if (!input->variable->equalstyle(eindex[m]))
        error->all(FLERR, "Variable for fix wall is invalid style");
    }
    if (sstyle[m] == VARIABLE) {
      sindex[m] = input->variable->find(sstr[m]);
      if (sindex[m] < 0) error->all(FLERR, "Variable name for fix wall does not exist");
      if (!input->variable->equalstyle(sindex[m]))
        error->all(FLERR, "Variable for fix wall is invalid style");
    }
  }

  // setup coefficients

  for (int m = 0; m < nwall; m++) precompute(m);

  if (utils::strmatch(update->integrate_style, "^respa")) {
    ilevel_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels - 1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level, ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixWall::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style, "^verlet")) {
    if (!fldflag) post_force(vflag);
  } else {
    (dynamic_cast<Respa *>(update->integrate))->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag, ilevel_respa, 0);
    (dynamic_cast<Respa *>(update->integrate))->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixWall::min_setup(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   only called if fldflag set, in place of post_force
------------------------------------------------------------------------- */

void FixWall::pre_force(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixWall::post_force(int vflag)
{
  // virial setup

  v_init(vflag);

  // energy intialize.
  // eflag is used to track whether wall energies have been communitcated.

  eflag = 0;
  for (int m = 0; m <= nwall; m++) ewall[m] = 0.0;

  // coord = current position of wall
  // evaluate variables if necessary, wrap with clear/add
  // for epsilon/sigma variables need to re-invoke precompute()

  if (varflag) modify->clearstep_compute();

  double coord;
  for (int m = 0; m < nwall; m++) {
    if (xstyle[m] == VARIABLE) {
      coord = input->variable->compute_equal(xindex[m]);
      if (wallwhich[m] < YLO)
        coord *= xscale;
      else if (wallwhich[m] < ZLO)
        coord *= yscale;
      else
        coord *= zscale;
    } else
      coord = coord0[m];
    if (wstyle[m] == VARIABLE) {
      if (estyle[m] == VARIABLE) {
        epsilon[m] = input->variable->compute_equal(eindex[m]);
        if (epsilon[m] < 0.0) error->all(FLERR, "Variable evaluation in fix wall gave bad value");
      }
      if (sstyle[m] == VARIABLE) {
        sigma[m] = input->variable->compute_equal(sindex[m]);
        if (sigma[m] < 0.0) error->all(FLERR, "Variable evaluation in fix wall gave bad value");
      }
      precompute(m);
    }

    wall_particle(m, wallwhich[m], coord);
  }

  if (varflag) modify->addstep_compute(update->ntimestep + 1);
}

/* ---------------------------------------------------------------------- */

void FixWall::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixWall::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   energy of wall interaction
------------------------------------------------------------------------- */

double FixWall::compute_scalar()
{
  // only sum across procs one time

  if (eflag == 0) {
    MPI_Allreduce(ewall, ewall_all, nwall + 1, MPI_DOUBLE, MPI_SUM, world);
    eflag = 1;
  }
  return ewall_all[0];
}

/* ----------------------------------------------------------------------
   components of force on wall
------------------------------------------------------------------------- */

double FixWall::compute_vector(int n)
{
  // only sum across procs one time

  if (eflag == 0) {
    MPI_Allreduce(ewall, ewall_all, nwall + 1, MPI_DOUBLE, MPI_SUM, world);
    eflag = 1;
  }
  return ewall_all[n + 1];
}
