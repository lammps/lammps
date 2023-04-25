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

#include "fix_rheo_thermal.h"

#include "atom.h"
#include "comm.h"
#include "compute_rheo_grad.h"
#include "compute_rheo_vshift.h"
#include "domain.h"
#include "error.h"
#include "fix_rheo.h"
#include "force.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace RHEO_NS;
using namespace FixConst;
enum {NONE, CONSTANT, TYPE};

/* ---------------------------------------------------------------------- */

FixRHEOThermal::FixRHEOThermal(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), fix_rheo(nullptr), compute_grad(nullptr), compute_vshift(nullptr),
  Tc_type(nullptr), kappa_type(nullptr), cv_type(nullptr), conductivity(nullptr)
{
  if (narg < 4) error->all(FLERR,"Illegal fix command");

  Tc_style = NONE;
  cv_style = NONE;
  conductivity_style = NONE;

  comm_forward = 1;
  nmax_store = 0;

  int ntypes = atom->ntypes;
  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"conductivity") == 0) {
      // Conductivity arguments
      if (iarg + 1 >= narg) error->all(FLERR,"Insufficient arguments for conductivity option");
      if (strcmp(arg[iarg + 1],"constant") == 0) {
        if (iarg + 2 >= narg) error->all(FLERR,"Insufficient arguments for conductivity option");
        conductivity_style = CONSTANT;
        kappa = utils::numeric(FLERR,arg[iarg + 2],false,lmp);
        if (kappa < 0.0) error->all(FLERR,"The conductivity must be positive");
        iarg += 2;
      } else if (strcmp(arg[iarg + 1],"type") == 0) {
        if (iarg + 1 + ntypes >= narg) error->all(FLERR,"Insufficient arguments for conductivity option");
        conductivity_style = TYPE;
        memory->create(kappa_type,ntypes+1,"rheo_thermal:kappa_type");
        for (int i = 1; i <= ntypes; i++) {
          kappa_type[i] = utils::numeric(FLERR,arg[iarg + 1 + i],false,lmp);
          if (kappa_type[i] < 0.0) error->all(FLERR,"The conductivity must be positive");
        }
        iarg += 1 + ntypes;
      } else {
        error->all(FLERR,"Illegal fix command, {}", arg[iarg + 1]);
      }
    } else if (strcmp(arg[iarg],"cv") == 0) {
      // Cv arguments
      if (iarg + 1 >= narg) error->all(FLERR,"Insufficient arguments for cv option");
      if (strcmp(arg[iarg + 1],"constant") == 0) {
        if (iarg + 2 >= narg) error->all(FLERR,"Insufficient arguments for cv option");
        cv_style = CONSTANT;
        cv = utils::numeric(FLERR,arg[iarg + 2],false,lmp);
        if (cv < 0.0) error->all(FLERR,"The specific heat must be positive");
        iarg += 2;
      } else if (strcmp(arg[iarg + 1],"type") == 0) {
        if (iarg + 1 + ntypes >= narg) error->all(FLERR,"Insufficient arguments for cv option");
        cv_style = TYPE;
        memory->create(cv_type,ntypes + 1,"rheo_thermal:cv_type");
        for (int i = 1; i <= ntypes; i++) {
          cv_type[i] = utils::numeric(FLERR,arg[iarg + 1 + i],false,lmp);
          if (cv_type[i] < 0.0) error->all(FLERR,"The specific heat must be positive");
        }
        iarg += 1 + ntypes;
      } else {
        error->all(FLERR,"Illegal fix command, {}", arg[iarg + 1]);
      }
    } else if (strcmp(arg[iarg],"Tfreeze") == 0) {
      // T freeze arguments
      if (iarg+1 >= narg) error->all(FLERR,"Insufficient arguments for Tfreeze option");
      if (strcmp(arg[iarg + 1],"constant") == 0) {
        if (iarg + 2 >= narg) error->all(FLERR,"Insufficient arguments for Tfreeze option");
        Tc_style = CONSTANT;
        Tc = utils::numeric(FLERR,arg[iarg + 2],false,lmp);
        if (Tc < 0.0) error->all(FLERR,"The melting temperature must be positive");
        iarg += 2;
      } else if (strcmp(arg[iarg + 1],"type") == 0) {
        if (iarg + 1 + ntypes >= narg) error->all(FLERR,"Insufficient arguments for Tfreeze option");
        Tc_style = TYPE;
        memory->create(Tc_type,ntypes + 1,"rheo_thermal:Tc_type");
        for (int i = 1; i <= ntypes; i++) {
          Tc_type[i] = utils::numeric(FLERR,arg[iarg + 1 + i],false,lmp);
          if (Tc_type[i] < 0.0) error->all(FLERR,"The melting temperature must be positive");
        }
        iarg += 1 + ntypes;
      } else {
        error->all(FLERR,"Illegal fix command, {}", arg[iarg + 1]);
      }
    } else {
      error->all(FLERR,"Illegal fix command, {}", arg[iarg]);
    }
    iarg += 1;
  }

  if (cv_style == NONE || conductivity_style == NONE)
    error->all(FLERR, "Must specify specific heat and conductivity styles\n");
}

/* ---------------------------------------------------------------------- */

FixRHEOThermal::~FixRHEOThermal()
{
  // Remove custom property if it exists
  int tmp1, tmp2, index;
  index = atom->find_custom("rheo_conductivity", tmp1, tmp2);
  if (index != -1) atom->remove_custom(index, 1, 0);

  memory->destroy(cv_type);
  memory->destroy(Tc_type);
  memory->destroy(kappa_type);
}

/* ---------------------------------------------------------------------- */

int FixRHEOThermal::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= POST_INTEGRATE;
  mask |= POST_NEIGHBOR;
  mask |= PRE_FORCE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixRHEOThermal::init()
{
  auto fixes = modify->get_fix_by_style("rheo");
  if (fixes.size() == 0) error->all(FLERR, "Need to define fix rheo to use fix rheo/viscosity");
  fix_rheo = dynamic_cast<FixRHEO *>(fixes[0]);

  if (!fix_rheo->thermal_flag)
    error->all(FLERR, "Need to define thermal setting in fix rheo");
  compute_grad = fix_rheo->compute_grad;
  compute_vshift = fix_rheo->compute_vshift;

  dtf = 0.5 * update->dt * force->ftm2v;

  if (atom->temperature_flag != 1)
    error->all(FLERR,"fix rheo/thermal command requires atoms store temperature property");
  if (atom->heatflow_flag != 1)
    error->all(FLERR,"fix rheo/thermal command requires atoms store heatflow property");
}

/* ---------------------------------------------------------------------- */

void FixRHEOThermal::setup_pre_force(int /*vflag*/)
{
  fix_rheo->thermal_fix_defined = 1;

  // Identify whether this is the first/last instance of fix thermal
  // First will grow arrays, last will communicate
  first_flag = 0;
  last_flag = 0;

  int i = 0;
  auto fixlist = modify->get_fix_by_style("rheo/thermal");
  for (const auto &fix : fixlist) {
    if (strcmp(fix->id, id) == 0) break;
    i++;
  }

  if (i == 0) first_flag = 1;
  if ((i + 1) == fixlist.size()) last_flag = 1;

  // Create conductivity array if it doesn't already exist
  // Create a custom atom property so it works with compute property/atom
  // Do not create grow callback as there's no reason to copy/exchange data
  // Manually grow if nmax_store exceeded

  int tmp1, tmp2;
  int index = atom->find_custom("rheo_conductivity", tmp1, tmp2);
  if (index== -1) {
    index = atom->add_custom("rheo_conductivity", 1, 0);
    nmax_store = atom->nmax;
  }
  conductivity = atom->dvector[index];

  post_neighbor();
  pre_force(0);
}

/* ---------------------------------------------------------------------- */

void FixRHEOThermal::initial_integrate(int /*vflag*/)
{
  // update temperature from shifting
  if (!fix_rheo->shift_flag) return;
  int i, a;

  int *status = atom->status;
  int *mask = atom->mask;
  double *temperature = atom->temperature;
  double **gradt = compute_grad->gradt;
  double **vshift = compute_vshift->array_atom;

  int nlocal = atom->nlocal;
  int dim = domain->dimension;

  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  for (i = 0; i < nlocal; i++) {
    if (!(status[i] & STATUS_SHIFT)) continue;

    if (mask[i] & groupbit) {
      for (a = 0; a < dim; a++) {
        temperature[i] += dtv * vshift[i][a] * gradt[i][a];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixRHEOThermal::post_integrate()
{
  int *status = atom->status;
  double *temperature = atom->temperature;
  double *heatflow = atom->heatflow;
  double *rho = atom->rho;
  int *mask = atom->mask;
  int *type = atom->type;

  double cvi, Tci, Ti;

  //Integrate temperature and check status
  for (int i = 0; i < atom->nlocal; i++) {
    if (mask[i] & groupbit) {
      if (status[i] & STATUS_NO_FORCE) continue;

      cvi = calc_cv(i);
      temperature[i] += dtf * heatflow[i] / cvi;

      if (Tc_style != NONE) {
        Ti = temperature[i];
        if (Tc_style == CONSTANT) {
          Tci = Tc;
        } else if (Tc_style == TYPE) {
          Tci = Tc_type[type[i]];
        }

        if (Ti > Tci) {
          status[i] &= PHASEMASK;
          status[i] |= STATUS_FLUID;
        } else if (!(status[i] & STATUS_SOLID)) {
          status[i] &= PHASEMASK;
          status[i] |= STATUS_FREEZING;
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
  Only need to update non-evolving conductivity styles after atoms exchange
------------------------------------------------------------------------- */

void FixRHEOThermal::post_neighbor()
{
  int i;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  if (first_flag && (nmax_store < atom->nmax)) {
    memory->grow(conductivity, atom->nmax, "atom:rheo_conductivity");
    nmax_store = atom->nmax;
  }

  if (conductivity_style == CONSTANT) {
    for (i = 0; i < nall; i++)
      if (mask[i] & groupbit) conductivity[i] = kappa;
  } else if (conductivity_style == TYPE) {
    for (i = 0; i < nall; i++)
      if (mask[i] & groupbit) conductivity[i] = kappa_type[type[i]];
  }
}

/* ----------------------------------------------------------------------
  Update (and forward) evolving conductivity styles every timestep
  Zero heat flow
------------------------------------------------------------------------- */

void FixRHEOThermal::pre_force(int /*vflag*/)
{
  // send updated temperatures to ghosts if first instance of fix
  // then clear heatflow for next force calculation
  double *heatflow = atom->heatflow;
  if (first_flag) {
    comm->forward_comm(this);
    for (int i = 0; i < atom->nmax; i++) heatflow[i] = 0.0;
  }

  // Not needed yet, when needed add stage check for (un)pack_forward_comm() methods
  //int i;
  //double *conductivity = atom->dvector[index_cond];
  //int *mask = atom->mask;
  //int nlocal = atom->nlocal;

  //if (first_flag && (nmax_store < atom->nmax)) {
  //  memory->grow(conductivity, atom->nmax, "atom:rheo_conductivity");
  //  nmax_store = atom->nmax;
  //}

  //if (conductivity_style == TBD) {
  //  for (i = 0; i < nlocal; i++) {
  //    if (mask[i] & groupbit) {
  //    }
  //  }
  //}

  //if (last_flag && comm_forward) comm->forward_comm(this);
}

/* ---------------------------------------------------------------------- */

void FixRHEOThermal::final_integrate()
{
  double *temperature = atom->temperature;
  double *heatflow = atom->heatflow;
  int *status = atom->status;
  int *mask = atom->mask;

  double cvi;

  //Integrate temperature and check status
  for (int i = 0; i < atom->nlocal; i++) {
    if (mask[i] & groupbit) {
      if (status[i] & STATUS_NO_FORCE) continue;

      cvi = calc_cv(i);
      temperature[i] += dtf * heatflow[i] / cvi;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixRHEOThermal::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}

/* ---------------------------------------------------------------------- */

double FixRHEOThermal::calc_cv(int i)
{
  if (cv_style == CONSTANT) {
    return cv;
  } else if (cv_style == TYPE) {
    return(cv_type[atom->type[i]]);
  }
}

/* ---------------------------------------------------------------------- */

int FixRHEOThermal::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int i, j, m;

  double *temperature = atom->temperature;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = temperature[j];
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void FixRHEOThermal::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;

  double *temperature = atom->temperature;

  for (i = first; i < last; i++) temperature[i] = buf[m++];
}

/* ---------------------------------------------------------------------- */

int FixRHEOThermal::pack_reverse_comm(int n, int first, double *buf)
{
  int m = 0;
  int last = first + n;
  double *heatflow = atom->heatflow;

  for (int i = first; i < last; i++) {
    buf[m++] = heatflow[i];
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void FixRHEOThermal::unpack_reverse_comm(int n, int *list, double *buf)
{
  int m = 0;
  double *heatflow = atom->heatflow;

  for (int i = 0; i < n; i++)
    heatflow[list[i]] += buf[m++];
}

/* ---------------------------------------------------------------------- */

double FixRHEOThermal::memory_usage()
{
  double bytes = 0.0;
  bytes += (size_t) nmax_store * sizeof(double);
  return bytes;
}