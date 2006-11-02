/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

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
#include "fix_npt.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "output.h"
#include "modify.h"
#include "kspace.h"
#include "update.h"
#include "respa.h"
#include "temperature.h"
#include "pressure.h"
#include "domain.h"
#include "error.h"

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

/* ---------------------------------------------------------------------- */

FixNPT::FixNPT(int narg, char **arg) : Fix(narg, arg)
{
  if (narg < 7) error->all("Illegal fix npt command");

  restart_global = 1;

  t_start = atof(arg[3]);
  t_stop = atof(arg[4]);
  double t_period = atof(arg[5]);

  if (t_start < 0.0 || t_stop <= 0.0)
    error->all("Target T for fix npt cannot be 0.0");

  double p_period[3];
  if (strcmp(arg[6],"xyz") == 0) {
    if (narg < 10) error->all("Illegal fix npt command");

    press_couple = 0;
    p_start[0] = p_start[1] = p_start[2] = atof(arg[7]);
    p_stop[0] = p_stop[1] = p_stop[2] = atof(arg[8]);
    p_period[0] = p_period[1] = p_period[2] = atof(arg[9]);
    p_flag[0] = p_flag[1] = p_flag[2] = 1;

  } else {
    if (strcmp(arg[6],"xy") == 0) press_couple = 1;
    else if (strcmp(arg[6],"yz") == 0) press_couple = 2;
    else if (strcmp(arg[6],"xz") == 0) press_couple = 3;
    else if (strcmp(arg[6],"aniso") == 0) press_couple = 4;
    else error->all("Illegal fix npt command");

    if (narg < 14) error->all("Illegal fix npt command");

    if (strcmp(arg[7],"NULL") == 0) {
      p_start[0] = p_stop[0] = p_period[0] = 0.0;
      p_flag[0] = 0;
    } else {
      p_start[0] = atof(arg[7]);
      p_stop[0] = atof(arg[8]);
      p_flag[0] = 1;
    }
    if (strcmp(arg[9],"NULL") == 0) {
      p_start[1] = p_stop[1] = p_period[1] = 0.0;
      p_flag[1] = 0;
    } else {
      p_start[1] = atof(arg[9]);
      p_stop[1] = atof(arg[10]);
      p_flag[1] = 1;
    }
    if (strcmp(arg[11],"NULL") == 0) {
      p_start[2] = p_stop[2] = p_period[2] = 0.0;
      p_flag[2] = 0;
    } else {
      p_start[2] = atof(arg[11]);
      p_stop[2] = atof(arg[12]);
      p_flag[2] = 1;
    }

    double period = atof(arg[13]);
    if (p_flag[0]) p_period[0] = period;
    if (p_flag[1]) p_period[1] = period;
    if (p_flag[2]) p_period[2] = period;
  }

  // process extra keywords

  drag = 0.0;
  dilate_partial = 0;

  int iarg;
  if (press_couple == 0) iarg = 10;
  else iarg = 14;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"drag") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix npt command");
      drag = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"dilate") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix npt command");
      if (strcmp(arg[iarg+1],"all") == 0) dilate_partial = 0;
      else if (strcmp(arg[iarg+1],"partial") == 0) dilate_partial = 1;
      else error->all("Illegal fix npt command");
      iarg += 2;
    } else error->all("Illegal fix npt command");
  }

  // check for periodicity in controlled dimensions

  if (p_flag[0] && domain->xperiodic == 0)
    error->all("Cannot fix npt on a non-periodic dimension");
  if (p_flag[1] && domain->yperiodic == 0)
    error->all("Cannot fix npt on a non-periodic dimension");
  if (p_flag[2] && domain->zperiodic == 0)
    error->all("Cannot fix npt on a non-periodic dimension");

  // create a new temperature full style with fix ID and group all
  // pressure is always global (group all) and thus its
  //   KE/temperature contribution should use group all

  eta = eta_dot = 0.0;

  char **newarg = new char*[3];
  newarg[0] = id;
  newarg[1] = "all";
  newarg[2] = "full";
  force->add_temp(3,newarg,1);
  delete [] newarg;

  temperature = force->find_temp(id);

  // convert input periods to frequencies

  if (t_period <= 0.0 || (p_flag[0] && p_period[0] <= 0.0) || 
      (p_flag[1] && p_period[1] <= 0.0) || (p_flag[2] && p_period[2] <= 0.0))
    error->all("Fix npt periods must be > 0.0");

  t_freq = 1.0 / t_period;
  p_freq[0] = p_freq[1] = p_freq[2] = 0.0;
  if (p_flag[0]) p_freq[0] = 1.0 / p_period[0];
  if (p_flag[1]) p_freq[1] = 1.0 / p_period[1];
  if (p_flag[2]) p_freq[2] = 1.0 / p_period[2];

  // pressure init

  pressure = force->pressure;
  omega[0] = omega[1] = omega[2] = 0.0;
  omega_dot[0] = omega_dot[1] = omega_dot[2] = 0.0;

  nrigid = 0;
  rfix = NULL;
}

/* ---------------------------------------------------------------------- */

FixNPT::~FixNPT()
{
  delete [] rfix;
}

/* ---------------------------------------------------------------------- */

int FixNPT::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= THERMO;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNPT::init()
{
  if (atom->mass_require == 0)
    error->all("Cannot use fix npt with no per-type mass defined");

  // set timesteps and frequencies

  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
  dthalf = 0.5 * update->dt;

  double freq = MAX(p_freq[0],p_freq[1]);
  freq = MAX(freq,p_freq[2]);
  drag_factor = 1.0 - (update->dt * freq * drag);

  boltz = force->boltz;
  nktv2p = force->nktv2p;
  vol0 = domain->xprd * domain->yprd * domain->zprd;

  double mass = 0.0;
  for (int i = 0; i < atom->nlocal; i++) mass += atom->mass[atom->type[i]];
  MPI_Allreduce(&mass,&total_mass,1,MPI_DOUBLE,MPI_SUM,world);

  if (force->kspace) kspace_flag = 1;
  else kspace_flag = 0;

  if (strcmp(update->integrate_style,"respa") == 0) {
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
    step_respa = ((Respa *) update->integrate)->step;
  }

  // detect if any fix rigid exist so rigid bodies move when box is dilated
  // rfix[] = indices to each fix rigid

  delete [] rfix;
  nrigid = 0;
  rfix = NULL;

  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"rigid") == 0 ||
	strcmp(modify->fix[i]->style,"poems") == 0) nrigid++;
  if (nrigid) {
    rfix = new int[nrigid];
    nrigid = 0;
    for (int i = 0; i < modify->nfix; i++)
      if (strcmp(modify->fix[i]->style,"rigid") == 0 ||
	  strcmp(modify->fix[i]->style,"poems") == 0) rfix[nrigid++] = i;
  }
}

/* ----------------------------------------------------------------------
   compute T,P before integrator starts 
------------------------------------------------------------------------- */

void FixNPT::setup()
{
  t_target = t_start;                      // used by thermo_compute()
  p_target[0] = p_start[0];
  p_target[1] = p_start[1];
  p_target[2] = p_start[2];

  t_current = temperature->compute();
  pressure->compute(temperature);
  couple();
}

/* ----------------------------------------------------------------------
   1st half of Verlet update 
------------------------------------------------------------------------- */

void FixNPT::initial_integrate()
{
  int i;

  double delta = update->ntimestep - update->beginstep;
  delta /= update->endstep - update->beginstep;

  // update eta_dot

  t_target = t_start + delta * (t_stop-t_start);
  f_eta = t_freq*t_freq * (t_current/t_target - 1.0);
  eta_dot += f_eta*dthalf;
  eta_dot *= drag_factor;
  eta += dtv*eta_dot;

  // update omega_dot
  // for non-varying dims, p_freq is 0.0, so omega_dot doesn't change

  double f_omega;
  double denskt = (atom->natoms*boltz*t_target) / 
    (domain->xprd*domain->yprd*domain->zprd) * nktv2p;

  for (i = 0; i < 3; i++) {
    p_target[i] = p_start[i] + delta * (p_stop[i]-p_start[i]);
    f_omega = p_freq[i]*p_freq[i] * (p_current[i]-p_target[i])/denskt;
    omega_dot[i] += f_omega*dthalf;
    omega_dot[i] *= drag_factor;
    omega[i] += dtv*omega_dot[i];
    factor[i] = exp(-dthalf*(eta_dot+omega_dot[i]));
    dilation[i] = exp(dthalf*omega_dot[i]);
  }

  // v update only for atoms in NPT group

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double dtfm;
  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      dtfm = dtf / mass[type[i]];
      v[i][0] = v[i][0]*factor[0] + dtfm*f[i][0];
      v[i][1] = v[i][1]*factor[1] + dtfm*f[i][1];
      v[i][2] = v[i][2]*factor[2] + dtfm*f[i][2];
    }
  }

  // rescale simulation box and all owned atoms by 1/2 step

  box_dilate(0);

  // x update by full step only for atoms in NPT group

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];
    }
  }

  // rescale simulation box and all owned atoms by 1/2 step
  // redo KSpace coeffs since volume has changed

  box_dilate(0);
  if (kspace_flag) force->kspace->setup();
}

/* ----------------------------------------------------------------------
   2nd half of Verlet update 
------------------------------------------------------------------------- */

void FixNPT::final_integrate()
{
  int i;

  // v update only for atoms in NPT group

  double **v = atom->v;
  double **f = atom->f;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double dtfm;
  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      dtfm = dtf / mass[type[i]];
      v[i][0] = (v[i][0] + dtfm*f[i][0]) * factor[0];
      v[i][1] = (v[i][1] + dtfm*f[i][1]) * factor[1];
      v[i][2] = (v[i][2] + dtfm*f[i][2]) * factor[2];
    }
  }

  // compute new T,P

  t_current = temperature->compute();
  pressure->compute(temperature);
  couple();

  // update eta_dot

  f_eta = t_freq*t_freq * (t_current/t_target - 1.0);
  eta_dot += f_eta*dthalf;
  eta_dot *= drag_factor;

  // update omega_dot
  // for non-varying dims, p_freq is 0.0, so omega_dot doesn't change

  double f_omega;
  double denskt = (atom->natoms*boltz*t_target) / 
    (domain->xprd*domain->yprd*domain->zprd) * nktv2p;

  for (i = 0; i < 3; i++) {
    f_omega = p_freq[i]*p_freq[i] * (p_current[i]-p_target[i])/denskt;
    omega_dot[i] += f_omega*dthalf;
    omega_dot[i] *= drag_factor;
  }
}

/* ---------------------------------------------------------------------- */

void FixNPT::initial_integrate_respa(int ilevel, int flag)
{
  // if flag = 1, then is 2nd call at outermost level from rRESPA
  // perform 2nd half of box dilate on own + ghost atoms and return
  // redo KSpace coeffs since volume has changed

  if (flag == 1) {
    box_dilate(1);
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
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // outermost level - update eta_dot and omega_dot, apply to v, dilate box
  // all other levels - NVE update of v
  // x,v updates only performed for atoms in NPT group

  if (ilevel == nlevels_respa-1) {

    double delta = update->ntimestep - update->beginstep;
    delta /= update->endstep - update->beginstep;

    // update eta_dot

    t_target = t_start + delta * (t_stop-t_start);
    f_eta = t_freq*t_freq * (t_current/t_target - 1.0);
    eta_dot += f_eta*dthalf;
    eta_dot *= drag_factor;
    eta += dtv*eta_dot;

    // update omega_dot
    // for non-varying dims, p_freq is 0.0, so omega_dot doesn't change

    double f_omega;
    double denskt = (atom->natoms*boltz*t_target) / 
      (domain->xprd*domain->yprd*domain->zprd) * nktv2p;

    for (int i = 0; i < 3; i++) {
      p_target[i] = p_start[i] + delta * (p_stop[i]-p_start[i]);
      f_omega = p_freq[i]*p_freq[i] * (p_current[i]-p_target[i])/denskt;
      omega_dot[i] += f_omega*dthalf;
      omega_dot[i] *= drag_factor;
      omega[i] += dtv*omega_dot[i];
      factor[i] = exp(-dthalf*(eta_dot+omega_dot[i]));
      dilation[i] = exp(dthalf*omega_dot[i]);
    }

    // v update only for atoms in NPT group

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	dtfm = dtf / mass[type[i]];
	v[i][0] = v[i][0]*factor[0] + dtfm*f[i][0];
	v[i][1] = v[i][1]*factor[1] + dtfm*f[i][1];
	v[i][2] = v[i][2]*factor[2] + dtfm*f[i][2];
      }
    }

    // rescale simulation box and all owned atoms by 1/2 step

    box_dilate(0);

  } else {

    // v update only for atoms in NPT group

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	dtfm = dtf / mass[type[i]];
	v[i][0] += dtfm*f[i][0];
	v[i][1] += dtfm*f[i][1];
	v[i][2] += dtfm*f[i][2];
      }
    }
  }

  // innermost level - also update x only for atoms in NPT group

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

void FixNPT::final_integrate_respa(int ilevel)
{
  double dtfm;

  // set timesteps by level

  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  dthalf = 0.5 * step_respa[ilevel];

  // outermost level - update eta_dot and omega_dot,
  //   apply to v via final_integrate()
  // all other levels - NVE update of v
  // v update only performed for atoms in NPT group

  if (ilevel == nlevels_respa-1) final_integrate();
  else {

    // v update only for atoms in NPT group

    double **v = atom->v;
    double **f = atom->f;
    double *mass = atom->mass;
    int *type = atom->type;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

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

/* ---------------------------------------------------------------------- */

void FixNPT::couple()
{
  double *p_tensor = pressure->p_tensor;

  if (press_couple == 0)
    p_current[0] = p_current[1] = p_current[2] = pressure->p_total;
  else if (press_couple == 1) {
    double ave = 0.5 * (p_tensor[0] + p_tensor[1]);
    p_current[0] = p_current[1] = ave;
    p_current[2] = p_tensor[2];
  } else if (press_couple == 2) {
    double ave = 0.5 * (p_tensor[1] + p_tensor[2]);
    p_current[1] = p_current[2] = ave;
    p_current[0] = p_tensor[0];
  } else if (press_couple == 3) {
    double ave = 0.5 * (p_tensor[0] + p_tensor[2]);
    p_current[0] = p_current[2] = ave;
    p_current[1] = p_tensor[1];
  } if (press_couple == 4) {
    p_current[0] = p_tensor[0];
    p_current[1] = p_tensor[1];
    p_current[2] = p_tensor[2];
  }
}

/* ----------------------------------------------------------------------
   dilate the box around center of box
------------------------------------------------------------------------- */

void FixNPT::box_dilate(int flag)
{
  int i,n;

  // ctr = geometric center of box in a dimension
  // scale owned or owned+ghost atoms depending on flag
  // re-define simulation box via xprd/yprd/zprd
  // scale atom coords for all atoms or only for fix group atoms
  // if fix rigid exists, scale rigid body centers-of-mass
  // don't do anything if non-periodic or press style is constant volume

  double **x = atom->x;
  int *mask = atom->mask;
  if (flag) n = atom->nlocal + atom->nghost;
  else n = atom->nlocal;

  double oldlo,oldhi,ctr;

  if (domain->xperiodic && p_flag[0]) {
    oldlo = domain->boxxlo;
    oldhi = domain->boxxhi;
    ctr = 0.5 * (oldlo + oldhi);
    domain->boxxlo = (oldlo-ctr)*dilation[0] + ctr;
    domain->boxxhi = (oldhi-ctr)*dilation[0] + ctr;
    domain->xprd = domain->boxxhi - domain->boxxlo;
    if (dilate_partial) {
      for (i = 0; i < n; i++)
	if (mask[i] & groupbit)
	  x[i][0] = ctr + (x[i][0]-ctr)*dilation[0];
    } else {
      for (i = 0; i < n; i++)
	x[i][0] = ctr + (x[i][0]-ctr)*dilation[0];
    }
    if (nrigid)
      for (i = 0; i < nrigid; i++)
	modify->fix[rfix[i]]->
	  dilate(0,oldlo,oldhi,domain->boxxlo,domain->boxxhi);
  }

  if (domain->yperiodic && p_flag[1]) {
    oldlo = domain->boxylo;
    oldhi = domain->boxyhi;
    ctr = 0.5 * (oldlo + oldhi);
    domain->boxylo = (oldlo-ctr)*dilation[1] + ctr;
    domain->boxyhi = (oldhi-ctr)*dilation[1] + ctr;
    domain->yprd = domain->boxyhi - domain->boxylo;
    if (dilate_partial) {
      for (i = 0; i < n; i++)
	if (mask[i] & groupbit)
	  x[i][1] = ctr + (x[i][1]-ctr)*dilation[1];
    } else {
      for (i = 0; i < n; i++)
	x[i][1] = ctr + (x[i][1]-ctr)*dilation[1];
    }
    if (nrigid)
      for (i = 0; i < nrigid; i++)
	modify->fix[rfix[i]]->
	  dilate(1,oldlo,oldhi,domain->boxylo,domain->boxyhi);
  }

  if (domain->zperiodic && p_flag[2]) {
    oldlo = domain->boxzlo;
    oldhi = domain->boxzhi;
    ctr = 0.5 * (oldlo + oldhi);
    domain->boxzlo = (oldlo-ctr)*dilation[2] + ctr;
    domain->boxzhi = (oldhi-ctr)*dilation[2] + ctr;
    domain->zprd = domain->boxzhi - domain->boxzlo;
    if (dilate_partial) {
      for (i = 0; i < n; i++)
	if (mask[i] & groupbit)
	  x[i][2] = ctr + (x[i][2]-ctr)*dilation[2];
    } else {
      for (i = 0; i < n; i++)
	x[i][2] = ctr + (x[i][2]-ctr)*dilation[2];
    }
    if (nrigid)
      for (i = 0; i < nrigid; i++)
	modify->fix[rfix[i]]->
	  dilate(2,oldlo,oldhi,domain->boxzlo,domain->boxzhi);
  }
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write 
------------------------------------------------------------------------- */

void FixNPT::write_restart(FILE *fp)
{
  int n = 0;
  double list[8];
  list[n++] = eta;
  list[n++] = eta_dot;
  list[n++] = omega[0];
  list[n++] = omega[1];
  list[n++] = omega[2];
  list[n++] = omega_dot[0];
  list[n++] = omega_dot[1];
  list[n++] = omega_dot[2];

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(&list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix 
------------------------------------------------------------------------- */

void FixNPT::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;
  eta = list[n++];
  eta_dot = list[n++];
  omega[0] = list[n++];
  omega[1] = list[n++];
  omega[2] = list[n++];
  omega_dot[0] = list[n++];
  omega_dot[1] = list[n++];
  omega_dot[2] = list[n++];
}

/* ---------------------------------------------------------------------- */

int FixNPT::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"temp") == 0) {
    if (narg < 2) error->all("Illegal fix_modify command");
    temperature = force->find_temp(arg[1]);
    if (temperature == NULL)
      error->all("Could not find fix_modify temperature ID");
    if (temperature->igroup != 0 && comm->me == 0)
      error->warning("Temperature for NPT is not for group all");
    if (strcmp(temperature->style,"region") == 0 && comm->me == 0)
      error->warning("Temperature for NPT is style region");
    return 2;
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

int FixNPT::thermo_fields(int n, int *flags, char **keywords)
{
  if (n == 0) return 1;
  flags[0] = 3;
  strcpy(keywords[0],"EngNPT");
  return 1;
}

/* ---------------------------------------------------------------------- */

int FixNPT::thermo_compute(double *values)
{
  double ke = temperature->dof * boltz * t_target;
  double keplus = atom->natoms * boltz * t_target;
  double volume = domain->xprd * domain->yprd * domain->zprd;
  int pdim = p_flag[0] + p_flag[1] + p_flag[2];

  values[0] = ke * (eta + 0.5*eta_dot*eta_dot/(t_freq*t_freq));
  for (int i = 0; i < 3; i++)
    if (p_freq[i] > 0.0)
      values[0] += 0.5*keplus*omega_dot[i]*omega_dot[i] / 
	(p_freq[i]*p_freq[i]) + p_target[i]*(volume-vol0) / (pdim*nktv2p);
  return 1;
}
