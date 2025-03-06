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
   Joel Clemmer (SNL), Thomas O'Connor (CMU)
----------------------------------------------------------------------- */

#include "fix_rheo_viscosity.h"

#include "atom.h"
#include "comm.h"
#include "compute_rheo_grad.h"
#include "domain.h"
#include "error.h"
#include "fix_rheo.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace FixConst;
enum { NONE, CONSTANT, POWER };

/* ---------------------------------------------------------------------- */

FixRHEOViscosity::FixRHEOViscosity(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), eta(nullptr), npow(nullptr), K(nullptr), gd0(nullptr), tau0(nullptr),
    viscosity_style(nullptr), fix_rheo(nullptr), compute_grad(nullptr)
{
  if (narg < 4) error->all(FLERR, "Illegal fix command");

  comm_forward = 0;
  constant_flag = 0;
  evolve_flag = 0;

  // Currently can only have one instance of fix rheo/viscosity
  if (igroup != 0) error->all(FLERR, "fix rheo/viscosity command requires group all");

  int i, nlo, nhi;
  int n = atom->ntypes;
  memory->create(viscosity_style, n + 1, "rheo:viscosity_style");
  memory->create(eta, n + 1, "rheo:eta");
  memory->create(gd0, n + 1, "rheo:gd0");
  memory->create(K, n + 1, "rheo:K");
  memory->create(npow, n + 1, "rheo:npow");
  memory->create(tau0, n + 1, "rheo:tau0");
  for (i = 1; i <= n; i++) viscosity_style[i] = NONE;

  int iarg = 3;
  while (iarg < narg) {
    utils::bounds(FLERR, arg[iarg], 1, n, nlo, nhi, error);

    if (iarg + 1 >= narg) utils::missing_cmd_args(FLERR, "fix rheo/viscosity", error);

    if (strcmp(arg[iarg + 1], "constant") == 0) {
      if (iarg + 2 >= narg) utils::missing_cmd_args(FLERR, "fix rheo/viscosity constant", error);

      constant_flag = 1;
      double eta_one = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      if (eta_one < 0.0) error->all(FLERR, "The viscosity must be positive");
      iarg += 1;

      for (i = nlo; i <= nhi; i++) {
        viscosity_style[i] = CONSTANT;
        eta[i] = eta_one;
      }
    } else if (strcmp(arg[iarg + 1], "power") == 0) {
      if (iarg + 5 >= narg) utils::missing_cmd_args(FLERR, "fix rheo/viscosity power", error);

      comm_forward = 1;
      evolve_flag = 1;
      double eta_one = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      double gd0_one = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      double K_one = utils::numeric(FLERR, arg[iarg + 4], false, lmp);
      double npow_one = utils::numeric(FLERR, arg[iarg + 5], false, lmp);
      if (eta_one < 0.0) error->all(FLERR, "The viscosity must be positive");
      iarg += 4;

      for (i = nlo; i <= nhi; i++) {
        viscosity_style[i] = POWER;
        eta[i] = eta_one;
        gd0[i] = gd0_one;
        K[i] = K_one;
        npow[i] = npow_one;
        tau0[i] = eta[i] * gd0[i] - K[i] * pow(gd0[i], npow[i]);
      }
    } else {
      error->all(FLERR, "Illegal fix command, {}", arg[iarg]);
    }
    iarg += 2;
  }

  for (i = 1; i <= n; i++)
    if (viscosity_style[i] == NONE)
      error->all(FLERR, "Must specify viscosity for atom type {} in fix/rheo/viscosity", i);
}

/* ---------------------------------------------------------------------- */

FixRHEOViscosity::~FixRHEOViscosity()
{
  memory->destroy(viscosity_style);
  memory->destroy(eta);
  memory->destroy(gd0);
  memory->destroy(K);
  memory->destroy(npow);
  memory->destroy(tau0);
}

/* ---------------------------------------------------------------------- */

int FixRHEOViscosity::setmask()
{
  int mask = 0;
  mask |= POST_NEIGHBOR;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixRHEOViscosity::init()
{
  auto fixes = modify->get_fix_by_style("^rheo$");
  if (fixes.size() == 0) error->all(FLERR, "Need to define fix rheo to use fix rheo/viscosity");
  fix_rheo = dynamic_cast<FixRHEO *>(fixes[0]);

  compute_grad = fix_rheo->compute_grad;
}

/* ---------------------------------------------------------------------- */

void FixRHEOViscosity::setup_pre_force(int /*vflag*/)
{
  fix_rheo->viscosity_fix_defined = 1;

  if (modify->get_fix_by_style("rheo/viscosity").size() > 1)
    error->all(FLERR, "More than one fix rheo/viscosity defined");

  post_neighbor();
  pre_force(0);
}

/* ----------------------------------------------------------------------
  Only need to update non-evolving viscosity styles after atoms exchange
------------------------------------------------------------------------- */

void FixRHEOViscosity::post_neighbor()
{
  if (!constant_flag) return;

  int i, itype;
  int *type = atom->type;
  double *viscosity = atom->viscosity;

  int nall = atom->nlocal + atom->nghost;

  for (i = 0; i < nall; i++) {
    itype = type[i];
    if (viscosity_style[itype]) viscosity[i] = eta[itype];
  }
}

/* ----------------------------------------------------------------------
  Update (and forward) evolving viscosity styles every timestep
------------------------------------------------------------------------- */

void FixRHEOViscosity::pre_force(int /*vflag*/)
{
  if (!evolve_flag) return;

  int i, itype, a, b;
  double tmp, gdot;

  int *type = atom->type;
  double *viscosity = atom->viscosity;
  double **gradv = compute_grad->gradv;

  int nlocal = atom->nlocal;
  int dim = domain->dimension;

  for (i = 0; i < nlocal; i++) {
    itype = type[i];
    if (viscosity_style[itype] == POWER) {
      gdot = 0.0;
      for (a = 0; a < dim; a++) {
        for (b = a; b < dim; b++) {
          tmp = gradv[i][a * dim + b] + gradv[i][b * dim + a];
          tmp = tmp * tmp;
          if (a == b) tmp *= 0.5;
          gdot += tmp;
        }
      }
      gdot = sqrt(gdot);
      if (gdot <= gd0[itype]) {
        viscosity[i] = eta[itype];
      } else {
        viscosity[i] = K[itype] * pow(gdot, npow[itype] - 1) + tau0[itype] / gdot;
      }
    }
  }

  if (comm_forward) comm->forward_comm(this);
}

/* ---------------------------------------------------------------------- */

int FixRHEOViscosity::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/,
                                        int * /*pbc*/)
{
  double *viscosity = atom->viscosity;
  int m = 0;
  for (int i = 0; i < n; i++) {
    int j = list[i];
    buf[m++] = viscosity[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixRHEOViscosity::unpack_forward_comm(int n, int first, double *buf)
{
  double *viscosity = atom->viscosity;
  int m = 0;
  int last = first + n;
  for (int i = first; i < last; i++) { viscosity[i] = buf[m++]; }
}
