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

/* ----------------------------------------------------------------------
   Contributing author: Mark Stevens (SNL)
------------------------------------------------------------------------- */

#include "string.h"
#include "stdlib.h"
#include "math.h"
#include "fix_nph.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "domain.h"
#include "modify.h"
#include "fix_deform.h"
#include "compute.h"
#include "kspace.h"
#include "update.h"
#include "respa.h"
#include "domain.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

enum{XYZ,XY,YZ,XZ,ANISO};

/* ---------------------------------------------------------------------- */

FixNPH::FixNPH(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4) error->all("Illegal fix nph command");

  restart_global = 1;
  box_change = 1;
  time_integrate = 1;
  scalar_flag = 1;
  scalar_vector_freq = 1;
  extscalar = 1;

  double p_period[3];
  if (strcmp(arg[3],"xyz") == 0) {
    if (narg < 7) error->all("Illegal fix nph command");

    press_couple = XYZ;
    p_start[0] = p_start[1] = p_start[2] = atof(arg[4]);
    p_stop[0] = p_stop[1] = p_stop[2] = atof(arg[5]);
    p_period[0] = p_period[1] = p_period[2] = atof(arg[6]);
    p_flag[0] = p_flag[1] = p_flag[2] = 1;
    if (domain->dimension == 2) {
      p_start[2] = p_stop[2] = p_period[2] = 0.0;
      p_flag[2] = 0;
    }

  } else {
    if (strcmp(arg[3],"xy") == 0) press_couple = XY;
    else if (strcmp(arg[3],"yz") == 0) press_couple = YZ;
    else if (strcmp(arg[3],"xz") == 0) press_couple = XZ;
    else if (strcmp(arg[3],"aniso") == 0) press_couple = ANISO;
    else error->all("Illegal fix nph command");

    if (narg < 11) error->all("Illegal fix nph command");

    if (domain->dimension == 2 && 
	(press_couple == XY || press_couple == YZ || press_couple == XZ))
      error->all("Invalid fix nph command for a 2d simulation");

    if (strcmp(arg[4],"NULL") == 0) {
      p_start[0] = p_stop[0] = p_period[0] = 0.0;
      p_flag[0] = 0;
    } else {
      p_start[0] = atof(arg[4]);
      p_stop[0] = atof(arg[5]);
      p_flag[0] = 1;
    }
    if (strcmp(arg[6],"NULL") == 0) {
      p_start[1] = p_stop[1] = p_period[1] = 0.0;
      p_flag[1] = 0;
    } else {
      p_start[1] = atof(arg[6]);
      p_stop[1] = atof(arg[7]);
      p_flag[1] = 1;
    }
    if (strcmp(arg[8],"NULL") == 0) {
      p_start[2] = p_stop[2] = p_period[2] = 0.0;
      p_flag[2] = 0;
    } else {
      if (domain->dimension == 2)
	error->all("Invalid fix nph command for a 2d simulation");
      p_start[2] = atof(arg[8]);
      p_stop[2] = atof(arg[9]);
      p_flag[2] = 1;
    }

    double period = atof(arg[10]);
    if (p_flag[0]) p_period[0] = period;
    if (p_flag[1]) p_period[1] = period;
    if (p_flag[2]) p_period[2] = period;
  }

  // process extra keywords

  drag = 0.0;
  allremap = 1;

  int iarg;
  if (press_couple == XYZ) iarg = 7;
  else iarg = 11;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"drag") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix nph command");
      drag = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"dilate") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix nph command");
      if (strcmp(arg[iarg+1],"all") == 0) allremap = 1;
      else if (strcmp(arg[iarg+1],"partial") == 0) allremap = 0;
      else error->all("Illegal fix nph command");
      iarg += 2;
    } else error->all("Illegal fix nph command");
  }

  // error checks

  if (press_couple == XY && (p_flag[0] == 0 || p_flag[1] == 0))
    error->all("Invalid fix nph command pressure settings");
  if (press_couple == YZ && (p_flag[1] == 0 || p_flag[2] == 0))
    error->all("Invalid fix nph command pressure settings");
  if (press_couple == XZ && (p_flag[0] == 0 || p_flag[2] == 0))
    error->all("Invalid fix nph command pressure settings");

  if (press_couple == XY && 
      (p_start[0] != p_start[1] || p_stop[0] != p_stop[1]))
    error->all("Invalid fix nph command pressure settings");
  if (press_couple == YZ && 
      (p_start[1] != p_start[2] || p_stop[1] != p_stop[2]))
    error->all("Invalid fix nph command pressure settings");
  if (press_couple == XZ && 
      (p_start[0] != p_start[2] || p_stop[0] != p_stop[2]))
    error->all("Invalid fix nph command pressure settings");

  if (p_flag[0] && domain->xperiodic == 0)
    error->all("Cannot use fix nph on a non-periodic dimension");
  if (p_flag[1] && domain->yperiodic == 0)
    error->all("Cannot use fix nph on a non-periodic dimension");
  if (p_flag[2] && domain->zperiodic == 0)
    error->all("Cannot use fix nph on a non-periodic dimension");

  // convert input periods to frequencies

  if ((p_flag[0] && p_period[0] <= 0.0) || 
      (p_flag[1] && p_period[1] <= 0.0) || (p_flag[2] && p_period[2] <= 0.0))
    error->all("Fix nph periods must be > 0.0");

  p_freq[0] = p_freq[1] = p_freq[2] = 0.0;
  if (p_flag[0]) p_freq[0] = 1.0 / p_period[0];
  if (p_flag[1]) p_freq[1] = 1.0 / p_period[1];
  if (p_flag[2]) p_freq[2] = 1.0 / p_period[2];

  // create a new compute temp style
  // id = fix-ID + temp
  // compute group = all since pressure is always global (group all)
  //   and thus its KE/temperature contribution should use group all

  int n = strlen(id) + 6;
  id_temp = new char[n];
  strcpy(id_temp,id);
  strcat(id_temp,"_temp");

  char **newarg = new char*[3];
  newarg[0] = id_temp;
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "temp";
  modify->add_compute(3,newarg);
  delete [] newarg;
  tflag = 1;

  // create a new compute pressure style
  // id = fix-ID + press, compute group = all
  // pass id_temp as 4th arg to pressure constructor

  n = strlen(id) + 7;
  id_press = new char[n];
  strcpy(id_press,id);
  strcat(id_press,"_press");

  newarg = new char*[4];
  newarg[0] = id_press;
  newarg[1] = (char *) "all";
  newarg[2] = (char *) "pressure";
  newarg[3] = id_temp;
  modify->add_compute(4,newarg);
  delete [] newarg;
  pflag = 1;

  // Nose/Hoover pressure init

  omega[0] = omega[1] = omega[2] = 0.0;
  omega_dot[0] = omega_dot[1] = omega_dot[2] = 0.0;

  nrigid = 0;
  rfix = NULL;
}

/* ---------------------------------------------------------------------- */

FixNPH::~FixNPH()
{
  delete [] rfix;

  // delete temperature and pressure if fix created them

  if (tflag) modify->delete_compute(id_temp);
  if (pflag) modify->delete_compute(id_press);
  delete [] id_temp;
  delete [] id_press;
}

/* ---------------------------------------------------------------------- */

int FixNPH::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= THERMO_ENERGY;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNPH::init()
{
  if (domain->triclinic) error->all("Cannot use fix nph with triclinic box");

  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"deform") == 0) {
      int *dimflag = ((FixDeform *) modify->fix[i])->dimflag;
      if ((p_flag[0] && dimflag[0]) || (p_flag[1] && dimflag[1]) || 
	  (p_flag[2] && dimflag[2]))
	error->all("Cannot use fix npt and fix deform on same dimension");
    }

  // set temperature and pressure ptrs

  int icompute = modify->find_compute(id_temp);
  if (icompute < 0) error->all("Temperature ID for fix nph does not exist");
  temperature = modify->compute[icompute];

  icompute = modify->find_compute(id_press);
  if (icompute < 0) error->all("Pressure ID for fix nph does not exist");
  pressure = modify->compute[icompute];

  // set timesteps and frequencies

  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  dthalf = 0.5 * update->dt;

  double freq = MAX(p_freq[0],p_freq[1]);
  freq = MAX(freq,p_freq[2]);
  drag_factor = 1.0 - (update->dt * freq * drag);

  boltz = force->boltz;
  nktv2p = force->nktv2p;
  dimension = domain->dimension;
  if (dimension == 3) vol0 = domain->xprd * domain->yprd * domain->zprd;
  else vol0 = domain->xprd * domain->yprd;

  if (force->kspace) kspace_flag = 1;
  else kspace_flag = 0;

  if (strcmp(update->integrate_style,"respa") == 0) {
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
    step_respa = ((Respa *) update->integrate)->step;
  }

  // detect if any rigid fixes exist so rigid bodies move when box is remapped
  // rfix[] = indices to each fix rigid

  delete [] rfix;
  nrigid = 0;
  rfix = NULL;

  for (int i = 0; i < modify->nfix; i++)
    if (modify->fix[i]->rigid_flag) nrigid++;
  if (nrigid) {
    rfix = new int[nrigid];
    nrigid = 0;
    for (int i = 0; i < modify->nfix; i++)
      if (modify->fix[i]->rigid_flag) rfix[nrigid++] = i;
  }
}

/* ----------------------------------------------------------------------
   compute T,P before integrator starts 
------------------------------------------------------------------------- */

void FixNPH::setup(int vflag)
{
  // Nkt = initial value for piston mass and energy conservation
  // cannot be done in init() b/c temperature cannot be called there
  // is b/c Modify::init() inits computes after fixes due to dof dependence
  // guesstimate a unit-dependent t_initial if actual T = 0.0

  double t_initial = temperature->compute_scalar();
  if (t_initial == 0.0) {
    if (strcmp(update->unit_style,"lj") == 0) t_initial = 1.0;
    else t_initial = 300.0;
  }
  nkt = atom->natoms * boltz * t_initial;

  p_target[0] = p_start[0];                 // used by compute_scalar()
  p_target[1] = p_start[1];
  p_target[2] = p_start[2];

  if (press_couple == XYZ) {
    double tmp = temperature->compute_scalar();
    tmp = pressure->compute_scalar();
  } else {
    temperature->compute_vector();
    pressure->compute_vector();
  }
  couple();

  // trigger virial computation on next timestep

  pressure->addstep(update->ntimestep+1);
}

/* ----------------------------------------------------------------------
   1st half of Verlet update 
------------------------------------------------------------------------- */

void FixNPH::initial_integrate(int vflag)
{
  int i;
  double dtfm;

  double delta = update->ntimestep - update->beginstep;
  delta /= update->endstep - update->beginstep;

  // update omega_dot
  // for non-varying dims, p_freq is 0.0, so omega doesn't change

  double f_omega,volume;
  if (dimension == 3) volume = domain->xprd*domain->yprd*domain->zprd;
  else volume = domain->xprd*domain->yprd;
  double denskt = nkt / volume * nktv2p;

  for (i = 0; i < 3; i++) {
    p_target[i] = p_start[i] + delta * (p_stop[i]-p_start[i]);
    f_omega = p_freq[i]*p_freq[i] * (p_current[i]-p_target[i])/denskt;
    omega_dot[i] += f_omega*dthalf;
    omega_dot[i] *= drag_factor;
    omega[i] += dtv*omega_dot[i];
    factor[i] = exp(-dthalf*omega_dot[i]);
    dilation[i] = exp(dthalf*omega_dot[i]);
  }

  // v update only for atoms in NPH group

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (rmass) {
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	dtfm = dtf / rmass[i];
	v[i][0] = v[i][0]*factor[0] + dtfm*f[i][0];
	v[i][1] = v[i][1]*factor[1] + dtfm*f[i][1];
	v[i][2] = v[i][2]*factor[2] + dtfm*f[i][2];
      }
    }
  } else {
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	dtfm = dtf / mass[type[i]];
	v[i][0] = v[i][0]*factor[0] + dtfm*f[i][0];
	v[i][1] = v[i][1]*factor[1] + dtfm*f[i][1];
	v[i][2] = v[i][2]*factor[2] + dtfm*f[i][2];
      }
    }
  }

  // remap simulation box and all owned atoms by 1/2 step

  remap(0);

  // x update by full step only for atoms in NPH group

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];
    }
  }

  // remap simulation box and all owned atoms by 1/2 step
  // redo KSpace coeffs since volume has changed

  remap(0);
  if (kspace_flag) force->kspace->setup();
}

/* ----------------------------------------------------------------------
   2nd half of Verlet update 
------------------------------------------------------------------------- */

void FixNPH::final_integrate()
{
  int i;
  double dtfm;

  // v update only for atoms in NPH group

  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (rmass) {
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	dtfm = dtf / rmass[i];
	v[i][0] = (v[i][0] + dtfm*f[i][0]) * factor[0];
	v[i][1] = (v[i][1] + dtfm*f[i][1]) * factor[1];
	v[i][2] = (v[i][2] + dtfm*f[i][2]) * factor[2];
      }
    }
  } else {
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	dtfm = dtf / mass[type[i]];
	v[i][0] = (v[i][0] + dtfm*f[i][0]) * factor[0];
	v[i][1] = (v[i][1] + dtfm*f[i][1]) * factor[1];
	v[i][2] = (v[i][2] + dtfm*f[i][2]) * factor[2];
      }
    }
  }

  // compute new pressure

  if (press_couple == XYZ) {
    double tmp = temperature->compute_scalar();
    tmp = pressure->compute_scalar();
  } else {
    temperature->compute_vector();
    pressure->compute_vector();
  }
  couple();

  // trigger virial computation on next timestep

  pressure->addstep(update->ntimestep+1);

  // update omega_dot
  // for non-varying dims, p_freq is 0.0, so omega_dot doesn't change

  double f_omega,volume;
  if (dimension == 3) volume = domain->xprd*domain->yprd*domain->zprd;
  else volume = domain->xprd*domain->yprd;
  double denskt = nkt / volume * nktv2p;

  for (i = 0; i < 3; i++) {
    f_omega = p_freq[i]*p_freq[i] * (p_current[i]-p_target[i])/denskt;
    omega_dot[i] += f_omega*dthalf;
    omega_dot[i] *= drag_factor;
  }
}

/* ---------------------------------------------------------------------- */

void FixNPH::initial_integrate_respa(int vflag, int ilevel, int flag)
{
  // if flag = 1, then is 2nd call at outermost level from rRESPA
  // perform 2nd half of box remap on own + ghost atoms and return
  // redo KSpace coeffs since volume has changed

  if (flag == 1) {
    remap(1);
    if (kspace_flag) force->kspace->setup();
    return;
  }

  // set timesteps by level

  double dtfm;
  dtv = step_respa[ilevel];
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  dthalf = 0.5 * step_respa[ilevel];

  // atom quantities

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // outermost level - update omega_dot, apply to v, remap box
  // all other levels - NVE update of v
  // x,v updates only performed for atoms in NPH group

  if (ilevel == nlevels_respa-1) {

    double delta = update->ntimestep - update->beginstep;
    delta /= update->endstep - update->beginstep;

    // update omega_dot
    // for non-varying dims, p_freq is 0.0, so omega doesn't change

    double f_omega,volume;
    if (dimension == 3) volume = domain->xprd*domain->yprd*domain->zprd;
    else volume = domain->xprd*domain->yprd;
    double denskt = nkt / volume * nktv2p;

    for (int i = 0; i < 3; i++) {
      p_target[i] = p_start[i] + delta * (p_stop[i]-p_start[i]);
      f_omega = p_freq[i]*p_freq[i] * (p_current[i]-p_target[i])/denskt;
      omega_dot[i] += f_omega*dthalf;
      omega_dot[i] *= drag_factor;
      omega[i] += dtv*omega_dot[i];
      factor[i] = exp(-dthalf*omega_dot[i]);
      dilation[i] = exp(dthalf*omega_dot[i]);
    }

    // v update only for atoms in NPH group

    if (rmass) {
      for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
	  dtfm = dtf / rmass[i];
	  v[i][0] = v[i][0]*factor[0] + dtfm*f[i][0];
	  v[i][1] = v[i][1]*factor[1] + dtfm*f[i][1];
	  v[i][2] = v[i][2]*factor[2] + dtfm*f[i][2];
	}
      }
    } else {
      for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
	  dtfm = dtf / mass[type[i]];
	  v[i][0] = v[i][0]*factor[0] + dtfm*f[i][0];
	  v[i][1] = v[i][1]*factor[1] + dtfm*f[i][1];
	  v[i][2] = v[i][2]*factor[2] + dtfm*f[i][2];
	}
      }
    }

    // remap simulation box and all owned atoms by 1/2 step

    remap(0);

  } else {

    // v update only for atoms in NPH group

    if (rmass) {
      for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
	  dtfm = dtf / rmass[i];
	  v[i][0] += dtfm*f[i][0];
	  v[i][1] += dtfm*f[i][1];
	  v[i][2] += dtfm*f[i][2];
	}
      }
    } else {
      for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
	  dtfm = dtf / mass[type[i]];
	  v[i][0] += dtfm*f[i][0];
	  v[i][1] += dtfm*f[i][1];
	  v[i][2] += dtfm*f[i][2];
	}
      }
    }
  }

  // innermost level - also update x only for atoms in NPH group

  if (ilevel == 0) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	x[i][0] += dtv * v[i][0];
	x[i][1] += dtv * v[i][1];
	x[i][2] += dtv * v[i][2];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNPH::final_integrate_respa(int ilevel)
{
  double dtfm;

  // set timesteps by level

  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  dthalf = 0.5 * step_respa[ilevel];

  // outermost level - update omega_dot, apply to v via final_integrate()
  // all other levels - NVE update of v
  // v update only performed for atoms in NPH group

  if (ilevel == nlevels_respa-1) final_integrate();
  else {

    // v update only for atoms in NPH group

    double **v = atom->v;
    double **f = atom->f;
    double *rmass = atom->rmass;
    double *mass = atom->mass;
    int *type = atom->type;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    if (rmass) {
      for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
	  dtfm = dtf / rmass[i];
	  v[i][0] += dtfm*f[i][0];
	  v[i][1] += dtfm*f[i][1];
	  v[i][2] += dtfm*f[i][2];
	}
      }
    } else {
      for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
	  dtfm = dtf / mass[type[i]];
	  v[i][0] += dtfm*f[i][0];
	  v[i][1] += dtfm*f[i][1];
	  v[i][2] += dtfm*f[i][2];
	}
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNPH::couple()
{
  double *tensor = pressure->vector;

  if (press_couple == XYZ)
    p_current[0] = p_current[1] = p_current[2] = pressure->scalar;
  else if (press_couple == XY) {
    double ave = 0.5 * (tensor[0] + tensor[1]);
    p_current[0] = p_current[1] = ave;
    p_current[2] = tensor[2];
  } else if (press_couple == YZ) {
    double ave = 0.5 * (tensor[1] + tensor[2]);
    p_current[1] = p_current[2] = ave;
    p_current[0] = tensor[0];
  } else if (press_couple == XZ) {
    double ave = 0.5 * (tensor[0] + tensor[2]);
    p_current[0] = p_current[2] = ave;
    p_current[1] = tensor[1];
  } else if (press_couple == ANISO) {
    p_current[0] = tensor[0];
    p_current[1] = tensor[1];
    p_current[2] = tensor[2];
  }
}

/* ----------------------------------------------------------------------
   change box size
   remap owned or owned+ghost atoms depending on flag
   remap all atoms or fix group atoms depending on allremap flag
   if rigid bodies exist, scale rigid body centers-of-mass
------------------------------------------------------------------------- */

void FixNPH::remap(int flag)
{
  int i,n;
  double oldlo,oldhi,ctr;

  double **x = atom->x;
  int *mask = atom->mask;
  if (flag) n = atom->nlocal + atom->nghost;
  else n = atom->nlocal;

  // convert pertinent atoms and rigid bodies to lamda coords

  if (allremap) domain->x2lamda(n);
  else {
    for (i = 0; i < n; i++)
      if (mask[i] & groupbit)
	domain->x2lamda(x[i],x[i]);
  }

  if (nrigid)
    for (i = 0; i < nrigid; i++)
      modify->fix[rfix[i]]->deform(0);

  // reset global and local box to new size/shape

  for (i = 0; i < 3; i++) {
    if (p_flag[i]) {
      oldlo = domain->boxlo[i];
      oldhi = domain->boxhi[i];
      ctr = 0.5 * (oldlo + oldhi);
      domain->boxlo[i] = (oldlo-ctr)*dilation[i] + ctr;
      domain->boxhi[i] = (oldhi-ctr)*dilation[i] + ctr;
    }
  }

  domain->set_global_box();
  domain->set_local_box();

  // convert pertinent atoms and rigid bodies back to box coords

  if (allremap) domain->lamda2x(n);
  else {
    for (i = 0; i < n; i++)
      if (mask[i] & groupbit)
	domain->lamda2x(x[i],x[i]);
  }

  if (nrigid)
    for (i = 0; i < nrigid; i++)
      modify->fix[rfix[i]]->deform(1);
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write 
------------------------------------------------------------------------- */

void FixNPH::write_restart(FILE *fp)
{
  int n = 0;
  double list[6];
  list[n++] = omega[0];
  list[n++] = omega[1];
  list[n++] = omega[2];
  list[n++] = omega_dot[0];
  list[n++] = omega_dot[1];
  list[n++] = omega_dot[2];

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix 
------------------------------------------------------------------------- */

void FixNPH::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;
  omega[0] = list[n++];
  omega[1] = list[n++];
  omega[2] = list[n++];
  omega_dot[0] = list[n++];
  omega_dot[1] = list[n++];
  omega_dot[2] = list[n++];
}

/* ---------------------------------------------------------------------- */

int FixNPH::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"temp") == 0) {
    if (narg < 2) error->all("Illegal fix_modify command");
    if (tflag) {
      modify->delete_compute(id_temp);
      tflag = 0;
    }
    delete [] id_temp;
    int n = strlen(arg[1]) + 1;
    id_temp = new char[n];
    strcpy(id_temp,arg[1]);

    int icompute = modify->find_compute(id_temp);
    if (icompute < 0) error->all("Could not find fix_modify temperature ID");
    temperature = modify->compute[icompute];

    if (temperature->tempflag == 0)
      error->all("Fix_modify temperature ID does not compute temperature");
    if (temperature->igroup != 0 && comm->me == 0)
      error->warning("Temperature for NPH is not for group all");

    // reset id_temp of pressure to new temperature ID
    
    icompute = modify->find_compute(id_press);
    if (icompute < 0) error->all("Pressure ID for fix npt does not exist");
    modify->compute[icompute]->reset_extra_compute_fix(id_temp);

    return 2;

  } else if (strcmp(arg[0],"press") == 0) {
    if (narg < 2) error->all("Illegal fix_modify command");
    if (pflag) {
      modify->delete_compute(id_press);
      pflag = 0;
    }
    delete [] id_press;
    int n = strlen(arg[1]) + 1;
    id_press = new char[n];
    strcpy(id_press,arg[1]);

    int icompute = modify->find_compute(id_press);
    if (icompute < 0) error->all("Could not find fix_modify pressure ID");
    pressure = modify->compute[icompute];

    if (pressure->pressflag == 0)
      error->all("Fix_modify pressure ID does not compute pressure");
    return 2;
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

double FixNPH::compute_scalar()
{
  double volume;
  if (dimension == 3) volume = domain->xprd * domain->yprd * domain->zprd;
  else volume = domain->xprd * domain->yprd;

  int pdim = p_flag[0] + p_flag[1] + p_flag[2];

  double energy = 0.0;
  for (int i = 0; i < 3; i++)
    if (p_freq[i] > 0.0)
      energy += 0.5*nkt*omega_dot[i]*omega_dot[i] / 
	(p_freq[i]*p_freq[i]) + p_target[i]*(volume-vol0) / (pdim*nktv2p);

  return energy;
}
