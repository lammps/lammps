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

#include "fix_rheo.h"

#include "atom.h"
#include "compute_rheo_grad.h"
#include "compute_rheo_interface.h"
#include "compute_rheo_kernel.h"
#include "compute_rheo_rhosum.h"
#include "compute_rheo_vshift.h"
#include "domain.h"
#include "error.h"
#include "fix_store_peratom.h"
#include "force.h"
#include "modify.h"
#include "update.h"
#include "utils.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixRHEO::FixRHEO(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), compute_grad(nullptr), compute_kernel(nullptr),
  compute_interface(nullptr), compute_rhosum(nullptr), compute_vshift(nullptr),
  fix_store_visc(nullptr), fix_store_pres(nullptr), fix_store_cond(nullptr),
  fix_store_surf(nullptr), fix_store_fp(nullptr), surface(nullptr), conductivity(nullptr),
  viscosity(nullptr), pressure(nullptr), f_pressure(nullptr)
{
  time_integrate = 1;

  viscosity_fix_defined = 0;
  pressure_fix_defined = 0;
  thermal_fix_defined = 0;
  surface_fix_defined = 0;

  thermal_flag = 0;
  rhosum_flag = 0;
  shift_flag = 0;
  interface_flag = 0;
  surface_flag = 0;

  rho0 = 1.0;
  csq = 1.0;

  if (igroup != 0)
    error->all(FLERR,"fix rheo command requires group all");

  if (atom->rho_flag != 1)
    error->all(FLERR,"fix rheo command requires atom_style with density");
  if (atom->status_flag != 1)
    error->all(FLERR,"fix rheo command requires atom_style with status");

  if (narg < 5)
    error->all(FLERR,"Insufficient arguments for fix rheo command");

  cut = utils::numeric(FLERR,arg[3],false,lmp);
  if (strcmp(arg[4],"Quintic") == 0) {
      kernel_style = QUINTIC;
  } else if (strcmp(arg[4],"CRK0") == 0) {
      kernel_style = CRK0;
  } else if (strcmp(arg[4],"CRK1") == 0) {
      kernel_style = CRK1;
  } else if (strcmp(arg[4],"CRK2") == 0) {
      kernel_style = CRK2;
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
    } else if (strcmp(arg[iarg],"interface/reconstruction") == 0) {
      interface_flag = 1;
    } else if (strcmp(arg[iarg],"rhosum") == 0) {
      rhosum_flag = 1;
      if(iarg + 1 >= narg) error->all(FLERR,"Illegal rhosum option in fix rheo");
      zmin_rhosum = utils::inumeric(FLERR,arg[iarg + 1],false,lmp);
      iarg += 1;
    } else if (strcmp(arg[iarg],"rho0") == 0) {
      if(iarg + 1 >= narg) error->all(FLERR,"Illegal rho0 option in fix rheo");
      rho0 = utils::numeric(FLERR,arg[iarg + 1],false,lmp);
      iarg += 1;
    } else if (strcmp(arg[iarg],"csq") == 0) {
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
  if (fix_store_visc) modify->delete_fix("rheo_store_visc");
  if (fix_store_pres) modify->delete_fix("rheo_store_pres");
  if (fix_store_surf) modify->delete_fix("rheo_store_surf");
  if (fix_store_cond) modify->delete_fix("rheo_store_cond");
  if (fix_store_fp) modify->delete_fix("rheo_store_fp");

  if (compute_kernel) modify->delete_compute("rheo_kernel");
  if (compute_grad) modify->delete_compute("rheo_grad");
  if (compute_interface) modify->delete_compute("rheo_interface");
  if (compute_rhosum) modify->delete_compute("rheo_rhosum");
  if (compute_vshift) modify->delete_compute("rheo_vshift");
}

/* ---------------------------------------------------------------------- */

void FixRHEO::post_constructor()
{
  compute_kernel = dynamic_cast<ComputeRHEOKernel *>(modify->add_compute(fmt::format("rheo_kernel all rheo/kernel {} {} {}", kernel_style, zmin_kernel, cut)));

  fix_store_visc = dynamic_cast<FixStorePeratom *>(modify->add_fix("rheo_store_visc all STORE/PERATOM 0 1"))
  fix_store_visc->disable = 1;
  viscosity = fix_store_visc->vstore;
  fix_store_pres = dynamic_cast<FixStorePeratom *>(modify->add_fix("rheo_store_pres all STORE/PERATOM 0 1"))
  pressure = fix_store_pres->vstore;
  fix_store_pres->disable = 1;


  std::string cmd = "rheo_grad all rheo/grad {} velocity rho viscosity";
  if (thermal_flag) cmd += "temperature";
  compute_grad = dynamic_cast<ComputeRHEOGrad *>(modify->add_compute(fmt::format(cmd, cut)));
  compute_grad->fix_rheo = this;

  if (rhosum_flag)
    compute_rhosum = dynamic_cast<ComputeRHEORhoSum *>(modify->add_compute(fmt::format("rheo_rhosum all rheo/rho/sum {} {}", cut, zmin_rhosum)));

  if (shift_flag)
    compute_vshift = dynamic_cast<ComputeRHEOVShift *>(modify->add_compute(fmt::format("rheo_vshift all rheo/vshift {}", cut)));

  if (surface_flag) {
    fix_store_surf = dynamic_cast<FixStorePeratom *>(modify->add_fix("rheo_store_surf all STORE/PERATOM 0 1"))
    surface = fix_store_surf->vstore;
    fix_store_surf->disable = 1;
  }

  if (thermal_flag) {
    fix_store_cond = dynamic_cast<FixStorePeratom *>(modify->add_fix("rheo_store_cond all STORE/PERATOM 0 1"))
    conductivity = fix_store_cond->vstore;
    fix_store_cond->disable = 1;
  }

  if (interface_flag) {
    compute_interface = dynamic_cast<ComputeRHEOInterface *>(modify->add_compute(fmt::format("rheo_interface all rheo/interface {}", cut)));

    fix_store_fp = dynamic_cast<FixStorePeratom *>(modify->add_fix("rheo_store_fp all STORE/PERATOM 0 3"))
    f_pressure = fix_store_fp->astore;
    fix_store_fp->disable = 1;
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

  if (modify->get_fix_by_style("rheo").size() > 1)
    error->all(FLERR,"Can only specify one instance of fix rheo");
}

/* ---------------------------------------------------------------------- */

void FixRHEO::setup_pre_force(int /*vflag*/)
{
  // Check to confirm no accessory fixes are yet defined
  // FixRHEO must be the first fix
  // Note: these fixes set this flag in setup_pre_force()
  if (viscosity_fix_defined || pressure_fix_defined || thermal_fix_defined || surface_fix_defined)
    error->all(FLERR, "Fix RHEO must be defined before all other RHEO fixes");

  pre_force(0);
}

/* ---------------------------------------------------------------------- */

void FixRHEO::setup()
{
  // Check to confirm all accessory fixes are defined
  // Does not ensure fixes correctly cover all atoms (could be a subset group)
  // Note: these fixes set this flag in setup_pre_force()
  if (!viscosity_fix_defined)
    error->all(FLERR, "Missing fix rheo/viscosity");

  if (!pressure_fix_defined)
    error->all(FLERR, "Missing fix rheo/pressure");

  if(!thermal_fix_defined && thermal_flag)
    error->all(FLERR, "Missing fix rheo/thermal");

  if(!surface_fix_defined && surface_flag)
    error->all(FLERR, "Missing fix rheo/surface");

  // Reset to zero for next run
  thermal_fix_defined = 0;
  viscosity_fix_defined = 0;
  pressure_fix_defined = 0;
  surface_fix_defined = 0;
}

/* ---------------------------------------------------------------------- */

void FixRHEO::initial_integrate(int /*vflag*/)
{
  // update v and x and rho of atoms in group
  int i, a, b;
  double dtfm, divu;
  int dim = domain->dimension;

  int *status = atom->status;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rho = atom->rho;
  double *drho = atom->drho;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int rmass_flag = atom->rmass_flag;

  double **gradr = compute_grad->gradr;
  double **gradv = compute_grad->gradv;
  double **vshift = compute_vshift->array_atom;

  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  //Density Half-step
  for (i = 0; i < nlocal; i++) {
    if (status[i] & STATUS_NO_FORCE) continue;

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
  compute_interface->store_forces(); // Need to save, wiped in exchange
  compute_interface->compute_peratom();
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
        if (status[i] & STATUS_NO_FORCE) continue;
        if (!(status[i] & STATUS_FLUID)) continue;

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
    compute_vshift->correct_surfaces();
    for (i = 0; i < nlocal; i++) {

      if (!(status[i] & STATUS_SHIFT)) continue;

      if (mask[i] & groupbit) {
        for (a = 0; a < dim; a++) {
          x[i][a] += dtv * vshift[i][a];
          for (b = 0; b < dim; b++) {
            v[i][a] += dtv * vshift[i][b] * gradv[i][a * dim + b];
          }
        }

        if (!rhosum_flag) {
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
  if (rhosum_flag)
    compute_rhosum->compute_peratom();

  compute_grad->forward_fields(); // also forwards v and rho for chi
  compute_kernel->compute_peratom();
  compute_interface->compute_peratom();

  compute_grad->compute_peratom();
  compute_grad->forward_gradients();

  if (shift_flag)
    compute_vshift->compute_peratom();

  // Remove extra shifting/no force options options
  int *status = atom->status;
  int nall = atom->nlocal + atom->nghost;
  for (int i = 0; i < nall; i++) {
    if (mask[i] & groupbit) {
      status[i] &= ~STATUS_NO_FORCE;

      if (status[i] & STATUS_FLUID)
        status[i] &= ~STATUS_SHIFT;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixRHEO::final_integrate() {
  int *status = atom->status;
  double **gradv = compute_grad->gradv;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;

  double *rho = atom->rho;
  double *drho = atom->drho;
  int *type = atom->type;
  int *mask = atom->mask;
  double *mass = atom->mass;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;
  double dtfm, divu;
  double *rmass = atom->rmass;
  int rmass_flag = atom->rmass_flag;
  int i, a;

  int dim = domain->dimension;

  // Update velocity
  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (status[i] & STATUS_NO_FORCE) continue;

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
        if (status[i] & STATUS_NO_FORCE) continue;
        if (!(status[i] & STATUS_FLUID)) continue;

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
