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
   Joel Clemmer (SNL), Thomas O'Connor (CMU), Eric Palermo (CMU)
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
enum {NONE, CONSTANT, TYPE, POWER};

/* ---------------------------------------------------------------------- */

FixRHEOViscosity::FixRHEOViscosity(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), fix_rheo(nullptr), eta_type(nullptr), viscosity(nullptr), compute_grad(nullptr)
{
  if (narg < 4) error->all(FLERR,"Illegal fix command");

  viscosity_style = NONE;

  comm_forward = 0;
  nmax_old = 0;

  int ntypes = atom->ntypes;
  int iarg = 3;
  if (strcmp(arg[iarg],"constant") == 0) {
    if (iarg + 1 >= narg) error->all(FLERR,"Insufficient arguments for viscosity option");
    viscosity_style = CONSTANT;
    eta = utils::numeric(FLERR,arg[iarg + 1],false,lmp);
    if (eta < 0.0) error->all(FLERR,"The viscosity must be positive");
    iarg += 1;
  } else if (strcmp(arg[iarg],"type") == 0) {
    if (iarg + ntypes >= narg) error->all(FLERR,"Insufficient arguments for viscosity option");
    viscosity_style = TYPE;
    memory->create(eta_type, ntypes + 1, "rheo_thermal:eta_type");
    for (int i = 1; i <= ntypes; i++) {
      eta_type[i] = utils::numeric(FLERR,arg[iarg + 1 + i], false, lmp);
      if (eta_type[i] < 0.0) error->all(FLERR,"The viscosity must be positive");
    }
    iarg += ntypes;
  } else if (strcmp(arg[iarg],"power") == 0) {
    if (iarg + 4 >= narg) error->all(FLERR,"Insufficient arguments for viscosity option");
    viscosity_style = POWER;
    comm_forward = 1;
    eta = utils::numeric(FLERR,arg[iarg + 1],false,lmp);
    gd0 = utils::numeric(FLERR,arg[iarg + 2],false,lmp);
    K = utils::numeric(FLERR,arg[iarg + 3],false,lmp);
    npow = utils::numeric(FLERR,arg[iarg + 4],false,lmp);
    tau0 = eta * gd0 - K * pow(gd0, npow);
    if (eta < 0.0) error->all(FLERR,"The viscosity must be positive");
    iarg += 5;
  } else {
    error->all(FLERR,"Illegal fix command, {}", arg[iarg]);
  }

  if (viscosity_style == NONE)
    error->all(FLERR,"Must specify viscosity style for fix/rheo/viscosity");
}

/* ---------------------------------------------------------------------- */

FixRHEOViscosity::~FixRHEOViscosity()
{
  // Remove custom property if it exists
  int tmp1, tmp2, index;
  index = atom->find_custom("rheo_viscosity", tmp1, tmp2);
  if (index != -1) atom->remove_custom(index_visc, 1, 0);

  memory->destroy(eta_type);
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
  auto fixes = modify->get_fix_by_style("rheo");
  if (fixes.size() == 0) error->all(FLERR, "Need to define fix rheo to use fix rheo/viscosity");
  fix_rheo = dynamic_cast<FixRHEO *>(fixes[0]);

  compute_grad = fix_rheo->compute_grad;
}

/* ---------------------------------------------------------------------- */

void FixRHEOViscosity::setup_pre_force(int /*vflag*/)
{
  fix_rheo->viscosity_fix_defined = 1;

  // Identify whether this is the first/last instance of fix viscosity
  // First will grow arrays, last will communicate
  first_flag = 0
  last_flag = 0;

  int i = 0;
  auto fixlist = modify->get_fix_by_style("rheo/viscosity");
  for (const auto &ifix : fixlist) {
    if (strcmp(ifix->id, id) == 0) break;
    i++;
  }

  if (i == 0) first_flag = 1;
  if ((i + 1) == fixlist.size()) last_flag = 1;

  // Create viscosity array if it doesn't already exist
  // Create a custom atom property so it works with compute property/atom
  // Do not create grow callback as there's no reason to copy/exchange data
  // Manually grow if nmax_old exceeded

  int tmp1, tmp2;
  int index = atom->find_custom("rheo_viscosity", tmp1, tmp2);
  if (index_visc == -1) {
    index = atom->add_custom("rheo_viscosity", 1, 0);
    nmax_old = atom->nmax;
  }
  viscosity = atom->dvector[index];

  post_neighbor();
  pre_force(0);
}

/* ----------------------------------------------------------------------
  Only need to update non-evolving viscosity styles after atoms exchange
------------------------------------------------------------------------- */

void FixRHEOViscosity::post_neighbor()
{
  int i;

  int *type = atom->type;
  int *mask = atom->mask;

  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  if (first_flag && (nmax_old < atom->nmax))
    memory->grow(viscosity, atom->nmax, "atom:rheo_viscosity");
  nmax_old = atom->nmax;

  if (viscosity_style == CONSTANT) {
    for (i = 0; i < nall; i++)
      if (mask[i] & groupbit) viscosity[i] = eta;
  } else if (viscosity_style == TYPE) {
    for (i = 0; i < nall; i++)
      if (mask[i] & groupbit) viscosity[i] = eta_type[type[i]];
    }
  }
}

/* ----------------------------------------------------------------------
  Update (and forward) evolving viscosity styles every timestep
------------------------------------------------------------------------- */

void FixRHEOViscosity::pre_force(int /*vflag*/)
{
  int i, a, b;
  double tmp, gdot;

  int *mask = atom->mask;
  double **gradv = compute_grad->gradv;

  int nlocal = atom->nlocal;
  int dim = domain->dimension;

  if (first_flag && (nmax_old < atom->nmax))
    memory->grow(viscosity, atom->nmax, "atom:rheo_viscosity");
  nmax_old = atom->nmax;

  if (viscosity_style == POWER) {
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
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
        if (gdot <= gd0) {
          viscosity[i] = eta;
        } else {
          viscosity[i] = K * pow(gdot, npow - 1) + tau0 / gdot;
        }
      }
    }
  }

  if (last_flag && comm_forward) comm->forward_comm(this);
}

/* ---------------------------------------------------------------------- */

int FixRHEOViscosity::pack_forward_comm(int n, int *list, double *buf,
                                        int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,k,m;
  m = 0;

  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = viscosity[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixRHEOViscosity::unpack_forward_comm(int n, int first, double *buf)
{
  int i, k, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    viscosity[i] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

double FixRHEOViscosity::memory_usage()
{
  double bytes = 0.0;
  bytes += (size_t) nmax_old * sizeof(double);
  return bytes;
}
