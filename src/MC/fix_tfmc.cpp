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
   Contributing author: Kristof Bal (University of Antwerp, Belgium)
------------------------------------------------------------------------- */

#include "fix_tfmc.h"
#include <mpi.h>
#include <cstring>
#include <cmath>
#include <cfloat>
#include "atom.h"
#include "force.h"
#include "update.h"
#include "group.h"
#include "random_mars.h"
#include "comm.h"
#include "domain.h"
#include "memory.h"
#include "modify.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixTFMC::FixTFMC(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  xd(NULL), rotflag(0), random_num(NULL)
{
  if (narg < 6) error->all(FLERR,"Illegal fix tfmc command");

  // although we are not doing MD, we would like to use tfMC as an MD "drop in"
  time_integrate = 1;

  d_max = force->numeric(FLERR,arg[3]);
  T_set = force->numeric(FLERR,arg[4]);
  seed = force->inumeric(FLERR,arg[5]);

  if (d_max <= 0) error->all(FLERR,"Fix tfmc displacement length must be > 0");
  if (T_set <= 0) error->all(FLERR,"Fix tfmc temperature must be > 0");
  if (seed <= 0) error->all(FLERR,"Illegal fix tfmc random seed");

  // additional keywords

  comflag = 0;
  rotflag = 0;

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"com") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix tfmc command");
      comflag = 1;
      xflag = force->inumeric(FLERR,arg[iarg+1]);
      yflag = force->inumeric(FLERR,arg[iarg+2]);
      zflag = force->inumeric(FLERR,arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"rot") == 0) {
      if (iarg+1 > narg) error->all(FLERR,"Illegal fix tfmc command");
      rotflag = 1;
      iarg += 1;
    } else error->all(FLERR,"Illegal fix tfmc command");
  }

  // error checks
  if (comflag)
    if (xflag < 0 || xflag > 1 || yflag < 0 || yflag > 1 ||
        zflag < 0 || zflag > 1)
      error->all(FLERR,"Illegal fix tfmc command");

  if (xflag + yflag + zflag == 0)
    comflag = 0;

  if (rotflag) {
    xd = NULL;
    nmax = -1;
  }

  random_num = new RanMars(lmp,seed + comm->me);
}

/* ---------------------------------------------------------------------- */

FixTFMC::~FixTFMC()
{
  delete random_num;
  if (rotflag) {
    memory->destroy(xd);
    xd = NULL;
    nmax = -1;
  }
}

/* ---------------------------------------------------------------------- */

int FixTFMC::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTFMC::init()
{
  // shake cannot be handled because it requires velocities
  // (and real MD in general)
  int has_shake = 0;
  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"shake") == 0) ++has_shake;

  if (has_shake > 0)
    error->all(FLERR,"Fix tfmc is not compatible with fix shake");

  // obtain lowest mass in the system
  // We do this here, in init(), rather than in initial_integrate().
  // This might seem somewhat odd: after all, another atom could be added with a
  // mass smaller than mass_min (in the case of a per-particle mass), so mass_min
  // should change during the run. However, this would imply that the overall
  // meaning of the input Delta is not very well-defined, because its meaning
  // can change during the run. So we'll assume all particle types (in terms of
  // possible masses) are defined before the run starts

  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double mass_min_local = DBL_MAX;
  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (mass_min_local > rmass[i]) mass_min_local = rmass[i];
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (mass_min_local > mass[type[i]]) mass_min_local = mass[type[i]];
      }
  }
  MPI_Allreduce(&mass_min_local,&mass_min,1,MPI_DOUBLE,MPI_MIN,world);
}

/* ---------------------------------------------------------------------- */

void FixTFMC::initial_integrate(int /*vflag*/)
{
  double boltz = force->boltz;
  double **x = atom->x;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  double massone;
  double masstotal;
  double xcm_d[3], xcm_dall[3];
  double d_i, xi;
  double gamma, gamma_exp, gamma_expi;
  double P_acc, P_ran;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // in case we wish to track (and zero) the com movement
  if (comflag) {
    xcm_d[0] = 0.0;
    xcm_d[1] = 0.0;
    xcm_d[2] = 0.0;
  }

  // displacement vector, needed to calculate (and zero) rotation
  if (rotflag && nmax < nlocal) {
    nmax = nlocal + 1;
    memory->destroy(xd);
    memory->create(xd,nmax,3,"tfmc:xd");
  }

  // generate displacements for each atom
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      d_i = d_max * pow(mass_min/massone, 0.25);
      for (int j = 0; j < 3; j++) {
        P_acc = 0.0;
        P_ran = 1.0;
        gamma = f[i][j] * d_i / (2.0*boltz*T_set);
        gamma_exp = exp(gamma);
        gamma_expi = 1.0/gamma_exp;
        // generate displacements according to the tfMC distribution
        while (P_acc < P_ran) {
          xi = 2.0*random_num->uniform() - 1.0;
          P_ran = random_num->uniform();
          if (xi < 0) {
            P_acc = exp(2.0*xi*gamma) * gamma_exp - gamma_expi;
            P_acc = P_acc / (gamma_exp - gamma_expi);
          } else if (xi > 0) {
            P_acc = gamma_exp - exp(2.0*xi*gamma) * gamma_expi;
            P_acc = P_acc / (gamma_exp - gamma_expi);
          } else {
            P_acc = 1.0;
          }
        }
        // displace
        x[i][j] += xi * d_i;
        if (comflag) xcm_d[j] += xi * d_i * massone;
        if (rotflag) xd[i][j] = xi * d_i;
      }
    }
  }

  // if post factum zeroing of linear or rotational motion
  if (comflag || rotflag) masstotal = group->mass(igroup);

  // zero com motion
  if (comflag == 1 && group->count(igroup) != 0) {
    MPI_Allreduce(xcm_d,xcm_dall,3,MPI_DOUBLE,MPI_SUM,world);
    if (masstotal > 0.0) {
      xcm_dall[0] /= masstotal;
      xcm_dall[1] /= masstotal;
      xcm_dall[2] /= masstotal;
    } else xcm_dall[0] = xcm_dall[1] = xcm_dall[2] = 0.0;

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        if (xflag) x[i][0] -= xcm_dall[0];
        if (yflag) x[i][1] -= xcm_dall[1];
        if (zflag) x[i][2] -= xcm_dall[2];
      }
    }
  }

  // zero rotation
  if (rotflag == 1 && group->count(igroup) != 0) {

    double dx, dy, dz;
    double unwrap[3];
    double cm[3], angmom[3], inertia[3][3], omega[3];
    tagint *image = atom->image;
    group->xcm(igroup,masstotal,cm);

    // to zero rotations, we can employ the same principles the
        // velocity command uses to zero the angular momentum. of course,
        // there is no (conserved) momentum in MC, but we can substitute
        // "velocities" by a displacement vector and proceed from there.
        // this of course requires "forking" group->angmom(), which is
        // what we do here.

    double p[3];
    p[0] = p[1] = p[2] = 0.0;

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        domain->unmap(x[i],image[i],unwrap);
        dx = unwrap[0] - cm[0];
        dy = unwrap[1] - cm[1];
        dz = unwrap[2] - cm[2];
        if (rmass) massone = rmass[i];
        else massone = mass[type[i]];
        p[0] += massone * (dy*xd[i][2] - dz*xd[i][1]);
        p[1] += massone * (dz*xd[i][0] - dx*xd[i][2]);
        p[2] += massone * (dx*xd[i][1] - dy*xd[i][0]);
      }
    }
    MPI_Allreduce(p,angmom,3,MPI_DOUBLE,MPI_SUM,world);
   // end "angmom" calculation

    group->inertia(igroup,cm,inertia);
    group->omega(angmom,inertia,omega);

    // now, get rid of the rotation
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        domain->unmap(x[i],image[i],unwrap);
        dx = unwrap[0] - cm[0];
        dy = unwrap[1] - cm[1];
        dz = unwrap[2] - cm[2];
        x[i][0] -= omega[1]*dz - omega[2]*dy;
        x[i][1] -= omega[2]*dx - omega[0]*dz;
        x[i][2] -= omega[0]*dy - omega[1]*dx;
      }
    }
  }
}
