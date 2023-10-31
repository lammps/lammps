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

#include "fix_rheo.h"

#include "atom.h"
#include "compute_rheo_grad.h"
#include "compute_rheo_interface.h"
#include "compute_rheo_surface.h"
#include "compute_rheo_kernel.h"
#include "compute_rheo_rho_sum.h"
#include "compute_rheo_vshift.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "modify.h"
#include "update.h"
#include "utils.h"

using namespace LAMMPS_NS;
using namespace RHEO_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixRHEO::FixRHEO(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), compute_grad(nullptr), compute_kernel(nullptr), compute_surface(nullptr),
  compute_interface(nullptr), compute_rhosum(nullptr), compute_vshift(nullptr)
{
  time_integrate = 1;

  viscosity_fix_defined = 0;
  pressure_fix_defined = 0;
  thermal_fix_defined = 0;

  thermal_flag = 0;
  rhosum_flag = 0;
  shift_flag = 0;
  interface_flag = 0;
  surface_flag = 0;

  rho0 = 1.0;
  csq = 1.0;

  if (igroup != 0)
    error->all(FLERR,"fix rheo command requires group all");

  if (atom->pressure_flag != 1)
    error->all(FLERR,"fix rheo command requires atom_style with pressure");
  if (atom->rho_flag != 1)
    error->all(FLERR,"fix rheo command requires atom_style with density");
  if (atom->viscosity_flag != 1)
    error->all(FLERR,"fix rheo command requires atom_style with viscosity");
  if (atom->status_flag != 1)
    error->all(FLERR,"fix rheo command requires atom_style with status");

  if (narg < 5)
    error->all(FLERR,"Insufficient arguments for fix rheo command");

  h = utils::numeric(FLERR,arg[3],false,lmp);
  cut = h;
  if (strcmp(arg[4],"quintic") == 0) {
      kernel_style = QUINTIC;
  } else if (strcmp(arg[4],"RK0") == 0) {
      kernel_style = RK0;
  } else if (strcmp(arg[4],"RK1") == 0) {
      kernel_style = RK1;
  } else if (strcmp(arg[4],"RK2") == 0) {
      kernel_style = RK2;
  } else error->all(FLERR,"Unknown kernel style {} in fix rheo", arg[4]);
  zmin_kernel = utils::numeric(FLERR,arg[5],false,lmp);

  int iarg = 6;
  while (iarg < narg){
    if (strcmp(arg[iarg],"shift") == 0) {
      shift_flag = 1;
    } else if (strcmp(arg[iarg],"thermal") == 0) {
      thermal_flag = 1;
    } else if (strcmp(arg[iarg],"surface/detection") == 0) {
      surface_flag = 1;
      if(iarg + 2 >= narg) error->all(FLERR,"Illegal surface/detection option in fix rheo");
      if (strcmp(arg[iarg + 1], "coordination") == 0) {
        surface_style = COORDINATION;
        zmin_surface = utils::inumeric(FLERR,arg[iarg + 2],false,lmp);
        zmin_splash = utils::inumeric(FLERR,arg[iarg + 3],false,lmp);
      } else if (strcmp(arg[iarg + 1], "divergence") == 0) {
        surface_style = DIVR;
        divr_surface = utils::numeric(FLERR,arg[iarg + 2],false,lmp);
        zmin_splash = utils::inumeric(FLERR,arg[iarg + 3],false,lmp);
      } else {
        error->all(FLERR,"Illegal surface/detection option in fix rheo, {}", arg[iarg + 1]);
      }

      iarg += 3;
    } else if (strcmp(arg[iarg],"interface/reconstruct") == 0) {
      interface_flag = 1;
    } else if (strcmp(arg[iarg],"rho/sum") == 0) {
      rhosum_flag = 1;
    } else if (strcmp(arg[iarg],"density") == 0) {
      if(iarg + 1 >= narg) error->all(FLERR,"Illegal rho0 option in fix rheo");
      rho0 = utils::numeric(FLERR,arg[iarg + 1],false,lmp);
      iarg += 1;
    } else if (strcmp(arg[iarg],"sound/squared") == 0) {
      if(iarg+1 >= narg) error->all(FLERR,"Illegal csq option in fix rheo");
      csq = utils::numeric(FLERR,arg[iarg + 1],false,lmp);
      iarg += 1;
    } else {
      error->all(FLERR, "Illegal fix rheo command: {}", arg[iarg]);
    }
    iarg += 1;
  }
}

/* ---------------------------------------------------------------------- */

FixRHEO::~FixRHEO()
{
  if (compute_kernel) modify->delete_compute("rheo_kernel");
  if (compute_grad) modify->delete_compute("rheo_grad");
  if (compute_interface) modify->delete_compute("rheo_interface");
  if (compute_surface) modify->delete_compute("compute_surface");
  if (compute_rhosum) modify->delete_compute("rheo_rhosum");
  if (compute_vshift) modify->delete_compute("rheo_vshift");
}


/* ----------------------------------------------------------------------
  Create necessary internal computes
------------------------------------------------------------------------- */

void FixRHEO::post_constructor()
{
  compute_kernel = dynamic_cast<ComputeRHEOKernel *>(modify->add_compute(
    fmt::format("rheo_kernel all RHEO/KERNEL {}", kernel_style)));
  compute_kernel->fix_rheo = this;

  std::string cmd = "rheo_grad all RHEO/GRAD velocity rho viscosity";
  if (thermal_flag) cmd += "temperature";
  compute_grad = dynamic_cast<ComputeRHEOGrad *>(modify->add_compute(cmd));
  compute_grad->fix_rheo = this;

  if (rhosum_flag) {
    compute_rhosum = dynamic_cast<ComputeRHEORhoSum *>(modify->add_compute(
      "rheo_rhosum all RHEO/RHO/SUM"));
    compute_rhosum->fix_rheo = this;
  }

  if (shift_flag) {
    compute_vshift = dynamic_cast<ComputeRHEOVShift *>(modify->add_compute(
      "rheo_vshift all RHEO/VSHIFT"));
    compute_vshift->fix_rheo = this;
  }

  if (interface_flag) {
    compute_interface = dynamic_cast<ComputeRHEOInterface *>(modify->add_compute(
      "rheo_interface all RHEO/INTERFACE"));
    compute_interface->fix_rheo = this;
  }

  if (surface_flag) {
    compute_surface = dynamic_cast<ComputeRHEOSurface *>(modify->add_compute(
      "rheo_surface all RHEO/SURFACE"));
    compute_surface->fix_rheo = this;
  }
}

/* ---------------------------------------------------------------------- */

int FixRHEO::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixRHEO::init()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;

  if (modify->get_fix_by_style("^rheo$").size() > 1)
    error->all(FLERR,"Can only specify one instance of fix rheo");
}

/* ---------------------------------------------------------------------- */

void FixRHEO::setup_pre_force(int /*vflag*/)
{
  // Check to confirm accessory fixes do not preceed FixRHEO
  // Note: fixes set this flag in setup_pre_force()
  if (viscosity_fix_defined || pressure_fix_defined || thermal_fix_defined)
    error->all(FLERR, "Fix RHEO must be defined before all other RHEO fixes");

  // Calculate surfaces
  if (surface_flag) compute_surface->compute_peratom();

  pre_force(0);
}

/* ---------------------------------------------------------------------- */

void FixRHEO::setup(int /*vflag*/)
{
  // Confirm all accessory fixes are defined
  // Note: fixes set this flag in setup_pre_force()
  if (!viscosity_fix_defined)
    error->all(FLERR, "Missing fix rheo/viscosity");

  if (!pressure_fix_defined)
    error->all(FLERR, "Missing fix rheo/pressure");

  if((!thermal_fix_defined) && thermal_flag)
    error->all(FLERR, "Missing fix rheo/thermal");

  // Reset to zero for future runs
  thermal_fix_defined = 0;
  viscosity_fix_defined = 0;
  pressure_fix_defined = 0;

  // Check fixes cover all atoms (may still fail if atoms are created)
  // FixRHEOPressure currently requires group all
  auto visc_fixes = modify->get_fix_by_style("rheo/viscosity");
  auto therm_fixes = modify->get_fix_by_style("rheo/thermal");

  int *mask = atom->mask;
  int v_coverage_flag = 1;
  int t_coverage_flag = 1;
  int covered;
  for (int i = 0; i < atom->nlocal; i++) {
    covered = 0;
    for (auto fix : visc_fixes)
      if (mask[i] & fix->groupbit) covered = 1;
    if (!covered) v_coverage_flag = 0;
    if (thermal_flag) {
      covered = 0;
      for (auto fix : therm_fixes)
        if (mask[i] & fix->groupbit) covered = 1;
      if (!covered) v_coverage_flag = 0;
    }
  }

  if (!v_coverage_flag)
    error->one(FLERR, "Fix rheo/viscosity does not fully cover all atoms");
  if (!t_coverage_flag)
    error->one(FLERR, "Fix rheo/thermal does not fully cover all atoms");
}

/* ---------------------------------------------------------------------- */

void FixRHEO::initial_integrate(int /*vflag*/)
{
  // update v, x and rho of atoms in group
  int i, a, b;
  double dtfm, divu;

  int *type = atom->type;
  int *mask = atom->mask;
  int *status = atom->status;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rho = atom->rho;
  double *drho = atom->drho;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  double **gradr = compute_grad->gradr;
  double **gradv = compute_grad->gradv;
  double **vshift;
  if (shift_flag)
    vshift = compute_vshift->vshift;

  int nlocal = atom->nlocal;
  int rmass_flag = atom->rmass_flag;
  int dim = domain->dimension;

  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  //Density Half-step
  for (i = 0; i < nlocal; i++) {
    if (status[i] & STATUS_NO_INTEGRATION) continue;

    if (mask[i] & groupbit) {
      if (rmass_flag) {
        dtfm = dtf / rmass[i];
      } else {
        dtfm = dtf / mass[type[i]];
      }

      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];
    }
  }

  // Update gradients and interpolate solid properties
  compute_grad->forward_fields(); // also forwards v and rho for chi
  if (interface_flag) {
    // Need to save, wiped in exchange
    compute_interface->store_forces();
    compute_interface->compute_peratom();
  }
  compute_grad->compute_peratom();

  // Position half-step
  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      for (a = 0; a < dim; a++) {
        x[i][a] += dtv * v[i][a];
      }
    }
  }

  // Update density using div(u)
  if (!rhosum_flag) {
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        if (status[i] & STATUS_NO_INTEGRATION) continue;
        if (status[i] & PHASECHECK) continue;

        divu = 0;
        for (a = 0; a < dim; a++) {
          divu += gradv[i][a * (1 + dim)];
        }
        rho[i] += dtf * (drho[i] - rho[i] * divu);
      }
    }
  }

  // Shifting atoms
  if (shift_flag) {
    compute_vshift->correct_surfaces(); // Could this be moved to preforce after the surface fix runs?
    for (i = 0; i < nlocal; i++) {

      if (status[i] & STATUS_NO_SHIFT) continue;

      if (mask[i] & groupbit) {
        for (a = 0; a < dim; a++) {
          x[i][a] += dtv * vshift[i][a];
          for (b = 0; b < dim; b++) {
            v[i][a] += dtv * vshift[i][b] * gradv[i][a * dim + b];
          }
        }

        if (!rhosum_flag) {
          if (status[i] & PHASECHECK) continue;
          for (a = 0; a < dim; a++) {
            rho[i] += dtv * vshift[i][a] * gradr[i][a];
          }
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixRHEO::pre_force(int /*vflag*/)
{
  compute_kernel->compute_coordination(); // Needed for rho sum

  if (rhosum_flag)
    compute_rhosum->compute_peratom();

  compute_kernel->compute_peratom();
  if (interface_flag) {
    // Note on first setup, have no forces for pressure to reference
    compute_interface->compute_peratom();
  }

  // No need to forward v, rho, or T for compute_grad since already done
  compute_grad->compute_peratom();
  compute_grad->forward_gradients();

  if (shift_flag)
    compute_vshift->compute_peratom();

  // Remove temporary options
  int *mask = atom->mask;
  int *status = atom->status;
  int nall = atom->nlocal + atom->nghost;
  for (int i = 0; i < nall; i++)
    if (mask[i] & groupbit)
      status[i] &= OPTIONSMASK;

  // Calculate surfaces, update status
  if (surface_flag) compute_surface->compute_peratom();
}

/* ---------------------------------------------------------------------- */

void FixRHEO::final_integrate()
{
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  double dtfm, divu;
  int i, a;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **gradv = compute_grad->gradv;
  double *rho = atom->rho;
  double *drho = atom->drho;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int *status = atom->status;

  int rmass_flag = atom->rmass_flag;
  int dim = domain->dimension;

  // Update velocity
  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (status[i] & STATUS_NO_INTEGRATION) continue;

      if (rmass_flag) {
        dtfm = dtf / rmass[i];
      } else {
        dtfm = dtf / mass[type[i]];
      }

      for (a = 0; a < dim; a++) {
        v[i][a] += dtfm * f[i][a];
      }
    }
  }

  // Update density using divu
  if (!rhosum_flag) {
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        if (status[i] & STATUS_NO_INTEGRATION) continue;
        if (status[i] & PHASECHECK) continue;

        divu = 0;
        for (a = 0; a < dim; a++) {
          divu += gradv[i][a * (1 + dim)];
        }
        rho[i] += dtf * (drho[i] - rho[i] * divu);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixRHEO::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}
