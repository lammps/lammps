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
   Contributing author: Laurent Joly (U Lyon, France), ljoly.ulyon@gmail.com
------------------------------------------------------------------------- */

#include "string.h"
#include "stdlib.h"
#include "fix_addtorque.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "group.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

FixAddTorque::FixAddTorque(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 6) error->all(FLERR,"Illegal fix addtorque command");

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;

  xstr = ystr = zstr = NULL;

  if (strstr(arg[3],"v_") == arg[3]) {
    int n = strlen(&arg[3][2]) + 1;
    xstr = new char[n];
    strcpy(xstr,&arg[3][2]);
  } else {
    xvalue = force->numeric(FLERR,arg[3]);
    xstyle = CONSTANT;
  }
  if (strstr(arg[4],"v_") == arg[4]) {
    int n = strlen(&arg[4][2]) + 1;
    ystr = new char[n];
    strcpy(ystr,&arg[4][2]);
  } else {
    yvalue = force->numeric(FLERR,arg[4]);
    ystyle = CONSTANT;
  }
  if (strstr(arg[5],"v_") == arg[5]) {
    int n = strlen(&arg[5][2]) + 1;
    zstr = new char[n];
    strcpy(zstr,&arg[5][2]);
  } else {
    zvalue = force->numeric(FLERR,arg[5]);
    zstyle = CONSTANT;
  }

  force_flag = 0;
  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixAddTorque::~FixAddTorque()
{
  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
}

/* ---------------------------------------------------------------------- */

int FixAddTorque::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAddTorque::init()
{
  // check variables

  if (xstr) {
    xvar = input->variable->find(xstr);
    if (xvar < 0) 
      error->all(FLERR,"Variable name for fix addtorque does not exist");
    if (input->variable->equalstyle(xvar)) xstyle = EQUAL;
    else error->all(FLERR,"Variable for fix addtorque is invalid style");
  }
  if (ystr) {
    yvar = input->variable->find(ystr);
    if (yvar < 0) 
      error->all(FLERR,"Variable name for fix addtorque does not exist");
    if (input->variable->equalstyle(yvar)) ystyle = EQUAL;
    else error->all(FLERR,"Variable for fix addtorque is invalid style");
  }
  if (zstr) {
    zvar = input->variable->find(zstr);
    if (zvar < 0) 
      error->all(FLERR,"Variable name for fix addtorque does not exist");
    if (input->variable->equalstyle(zvar)) zstyle = EQUAL;
    else error->all(FLERR,"Variable for fix addtorque is invalid style");
  }

  if (xstyle == EQUAL || ystyle == EQUAL || zstyle == EQUAL)
    varflag = EQUAL;
  else varflag = CONSTANT;

  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixAddTorque::setup(int vflag)
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

void FixAddTorque::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAddTorque::post_force(int vflag)
{
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int *type = atom->type;
  tagint *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;
  double mvv2e = force->mvv2e;

  double dx,dy,dz,vx,vy,vz,fx,fy,fz,massone,omegadotr;
  double tcm[3],xcm[3],angmom[3],omega[3],itorque[3],domegadt[3],tlocal[3];
  double inertia[3][3];
  double unwrap[3];

  // foriginal[0] = "potential energy" for added force
  // foriginal[123] = torque on atoms before extra force added

  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;
  force_flag = 0;

  if (varflag == EQUAL) {
    // variable torque, wrap with clear/add
    modify->clearstep_compute();
    if (xstyle == EQUAL) xvalue = input->variable->compute_equal(xvar);
    if (ystyle == EQUAL) yvalue = input->variable->compute_equal(yvar);
    if (zstyle == EQUAL) zvalue = input->variable->compute_equal(zvar);
    modify->addstep_compute(update->ntimestep + 1);
  }

  atom->check_mass();
  double masstotal = group->mass(igroup);
  group->xcm(igroup,masstotal,xcm);
  group->inertia(igroup,xcm,inertia);
  group->angmom(igroup,xcm,angmom);
  group->omega(angmom,inertia,omega);

  tlocal[0] = tlocal[1] = tlocal[2] = 0.0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - xcm[0];
      dy = unwrap[1] - xcm[1];
      dz = unwrap[2] - xcm[2];
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      omegadotr = omega[0]*dx+omega[1]*dy+omega[2]*dz;
      tlocal[0] += massone * omegadotr * (dy*omega[2] - dz*omega[1]);
      tlocal[1] += massone * omegadotr * (dz*omega[0] - dx*omega[2]);
      tlocal[2] += massone * omegadotr * (dx*omega[1] - dy*omega[0]);
    }
  MPI_Allreduce(tlocal,itorque,3,MPI_DOUBLE,MPI_SUM,world);

  tcm[0] = xvalue - mvv2e*itorque[0];
  tcm[1] = yvalue - mvv2e*itorque[1];
  tcm[2] = zvalue - mvv2e*itorque[2];
  group->omega(tcm,inertia,domegadt);

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      domain->unmap(x[i],image[i],unwrap);
      dx = unwrap[0] - xcm[0];
      dy = unwrap[1] - xcm[1];
      dz = unwrap[2] - xcm[2];
      vx = mvv2e*(dz*omega[1]-dy*omega[2]);
      vy = mvv2e*(dx*omega[2]-dz*omega[0]);
      vz = mvv2e*(dy*omega[0]-dx*omega[1]);
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      fx = massone * (dz*domegadt[1]-dy*domegadt[2] + vz*omega[1]-vy*omega[2]);
      fy = massone * (dx*domegadt[2]-dz*domegadt[0] + vx*omega[2]-vz*omega[0]);
      fz = massone * (dy*domegadt[0]-dx*domegadt[1] + vy*omega[0]-vx*omega[1]);

      // potential energy = - x dot f
      foriginal[0] -= fx*x[i][0] + fy*x[i][1] + fz*x[i][2];
      foriginal[1] += dy*f[i][2] - dz*f[i][1];
      foriginal[2] += dz*f[i][0] - dx*f[i][2];
      foriginal[3] += dx*f[i][1] - dy*f[i][0];
      f[i][0] += fx;
      f[i][1] += fy;
      f[i][2] += fz;
    }
}

/* ---------------------------------------------------------------------- */

void FixAddTorque::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAddTorque::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   potential energy of added torque
------------------------------------------------------------------------- */

double FixAddTorque::compute_scalar()
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return foriginal_all[0];
}

/* ----------------------------------------------------------------------
   return components of total torque on fix group before torque was changed
------------------------------------------------------------------------- */

double FixAddTorque::compute_vector(int n)
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return foriginal_all[n+1];
}
