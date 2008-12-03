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
   Contributing author: Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_spring.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "domain.h"
#include "error.h"
#include "group.h"

using namespace LAMMPS_NS;

enum{TETHER,COUPLE};

/* ---------------------------------------------------------------------- */

FixSpring::FixSpring(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 9) error->all("Illegal fix spring command");

  vector_flag = 1;
  size_vector = 3;
  scalar_vector_freq = 1;
  extvector = 1;

  if (strcmp(arg[3],"tether") == 0) {
    if (narg != 9) error->all("Illegal fix spring command");
    styleflag = TETHER;
    k_spring = atof(arg[4]);
    xflag = yflag = zflag = 1;
    if (strcmp(arg[5],"NULL") == 0) xflag = 0;
    else xc = atof(arg[5]);
    if (strcmp(arg[6],"NULL") == 0) yflag = 0;
    else yc = atof(arg[6]);
    if (strcmp(arg[7],"NULL") == 0) zflag = 0;
    else zc = atof(arg[7]);
    r0 = atof(arg[8]);
    if (r0 < 0) error->all("R0 < 0 for fix spring command");

  } else if (strcmp(arg[3],"couple") == 0) {
    if (narg != 10) error->all("Illegal fix spring command");
    styleflag = COUPLE;
    igroup2 = group->find(arg[4]);
    if (igroup2 == -1) 
      error->all("Could not find fix spring couple group ID"); 
    if (igroup2 == igroup) 
      error->all("Two groups cannot be the same in fix spring couple"); 
    group2bit = group->bitmask[igroup2];
    k_spring = atof(arg[5]);
    xflag = yflag = zflag = 1;
    if (strcmp(arg[6],"NULL") == 0) xflag = 0;
    else xc = atof(arg[6]);
    if (strcmp(arg[7],"NULL") == 0) yflag = 0;
    else yc = atof(arg[7]);
    if (strcmp(arg[8],"NULL") == 0) zflag = 0;
    else zc = atof(arg[8]);
    r0 = atof(arg[9]);
    if (r0 < 0) error->all("R0 < 0 for fix spring command");

  } else error->all("Illegal fix spring command");

  force_flag = 0;
  ftotal[0] = ftotal[1] = ftotal[2] = 0.0;
}

/* ---------------------------------------------------------------------- */

int FixSpring::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSpring::init()
{
  masstotal = group->mass(igroup);
  if (styleflag == COUPLE) masstotal2 = group->mass(igroup2);
  
  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixSpring::setup(int vflag)
{
  if (strcmp(update->integrate_style,"verlet") == 0)
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixSpring::post_force(int vflag)
{
  if (styleflag == TETHER) spring_tether();
  else spring_couple();
}

/* ---------------------------------------------------------------------- */

void FixSpring::spring_tether()
{
  double xcm[3];
  group->xcm(igroup,masstotal,xcm);

  // fx,fy,fz = components of k * (r-r0)
  
  double dx,dy,dz,fx,fy,fz,r,dr;
  
  dx = xcm[0] - xc;
  dy = xcm[1] - yc;
  dz = xcm[2] - zc;
  if (!xflag) dx = 0.0;
  if (!yflag) dy = 0.0;
  if (!zflag) dz = 0.0;
  r = sqrt(dx*dx + dy*dy + dz*dz);
  dr = r - r0;
  fx = k_spring*dx*dr/r;
  fy = k_spring*dy*dr/r;
  fz = k_spring*dz*dr/r;
  
  // apply restoring force to atoms in group
  // f = -k*(r-r0)*mass/masstotal

  double **f = atom->f;
  int *mask = atom->mask;
  int *type = atom->type;
  double *mass = atom->mass;
  int nlocal = atom->nlocal;
  
  double massfrac;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      massfrac = mass[type[i]]/masstotal;
      f[i][0] -= fx*massfrac;
      f[i][1] -= fy*massfrac;
      f[i][2] -= fz*massfrac;
    }

  ftotal[0] = -fx;
  ftotal[1] = -fy;
  ftotal[2] = -fz;
  force_flag = 0;
}

/* ---------------------------------------------------------------------- */

void FixSpring::spring_couple()
{
  double xcm[3],xcm2[3];
  group->xcm(igroup,masstotal,xcm);
  group->xcm(igroup2,masstotal2,xcm2);
  
  // fx,fy,fz = components of k * (r-r0)
  
  double dx,dy,dz,fx,fy,fz,r,dr;
  
  dx = xcm2[0] - xcm[0] - xc;
  dy = xcm2[1] - xcm[1] - yc;
  dz = xcm2[2] - xcm[2] - zc;  
  if (!xflag) dx = 0.0;
  if (!yflag) dy = 0.0;
  if (!zflag) dz = 0.0;
  r = sqrt(dx*dx + dy*dy + dz*dz);
  dr = r - r0;
  
  fx = k_spring*dx*dr/r;
  fy = k_spring*dy*dr/r;
  fz = k_spring*dz*dr/r;
  
  // apply restoring force to atoms in group
  // f = -k*(r-r0)*mass/masstotal

  double **f = atom->f;
  int *mask = atom->mask;
  int *type = atom->type;
  double *mass = atom->mass;
  int nlocal = atom->nlocal;

  double massfrac;
  for (int i = 0; i < nlocal; i++) {         
    if (mask[i] & groupbit) {
      massfrac = mass[type[i]]/masstotal;
      f[i][0] += fx*massfrac;
      f[i][1] += fy*massfrac;
      f[i][2] += fz*massfrac;
    }
    if (mask[i] & group2bit) {
      massfrac = mass[type[i]]/masstotal2;
      f[i][0] -= fx*massfrac;
      f[i][1] -= fy*massfrac;
      f[i][2] -= fz*massfrac;
    }
  }

  ftotal[0] = -fx;
  ftotal[1] = -fy;
  ftotal[2] = -fz;
  force_flag = 0;
}

/* ---------------------------------------------------------------------- */

void FixSpring::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ----------------------------------------------------------------------
   return components of total spring force on fix group
------------------------------------------------------------------------- */

double FixSpring::compute_vector(int n)
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(ftotal,ftotal_all,3,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return ftotal_all[n];
}
