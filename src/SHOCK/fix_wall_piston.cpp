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

#include "fix_wall_piston.h"
#include <cmath>
#include <cstring>
#include "atom.h"
#include "domain.h"
#include "lattice.h"
#include "update.h"
#include "error.h"
#include "random_mars.h"
#include "force.h"
#include "comm.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixWallPiston::FixWallPiston(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), randomt(NULL), gfactor1(NULL), gfactor2(NULL)
{
  force_reneighbor = 1;
  next_reneighbor = -1;

  if (narg < 4) error->all(FLERR,"Illegal fix wall/piston command");

  randomt = NULL;
  gfactor1 = gfactor2 = NULL;
  tempflag = 0;
  scaleflag = 1;
  roughflag = 0;
  rampflag = 0;
  rampNL1flag = 0;
  rampNL2flag = 0;
  rampNL3flag = 0;
  rampNL4flag = 0;
  rampNL5flag = 0;
  t_target = z0 = vx = vy = vz = 0.0;
  xloflag = xhiflag = yloflag = yhiflag = zloflag = zhiflag = 0;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"xlo") == 0)
      error->all(FLERR,"Fix wall/piston command only available at zlo");
    else if (strcmp(arg[iarg],"ylo") == 0)
      error->all(FLERR,"Fix wall/piston command only available at zlo");
    else if (strcmp(arg[iarg],"zlo") == 0) {
      zloflag = 1;
      iarg++;
      if (domain->boundary[2][0] != 2)
        error->all(FLERR,"Must shrink-wrap piston boundary");
    } else if (strcmp(arg[iarg],"xhi") == 0)
      error->all(FLERR,"Fix wall/piston command only available at zlo");
    else if (strcmp(arg[iarg],"yhi") == 0)
      error->all(FLERR,"Fix wall/piston command only available at zlo");
    else if (strcmp(arg[iarg],"zhi") == 0)
      error->all(FLERR,"Fix wall/piston command only available at zlo");
    else if (strcmp(arg[iarg],"vel") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix wall/piston command");
      vz = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"pos") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix wall/piston command");
      z0 = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"temp") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal fix wall/piston command");
      tempflag = 1;
      t_target = force->numeric(FLERR,arg[iarg+1]);
      t_period = force->numeric(FLERR,arg[iarg+2]);
      tseed    = force->inumeric(FLERR,arg[iarg+3]);
      t_extent = force->numeric(FLERR,arg[iarg+4]);
      if (t_target <= 0) error->all(FLERR,"Illegal fix wall/piston command");
      if (t_period <= 0) error->all(FLERR,"Illegal fix wall/piston command");
      if (t_extent <= 0) error->all(FLERR,"Illegal fix wall/piston command");
      if (tseed <= 0) error->all(FLERR,"Illegal fix wall/piston command");
      randomt = new RanMars(lmp,tseed + comm->me);
      gfactor1 = new double[atom->ntypes+1];
      gfactor2 = new double[atom->ntypes+1];
      iarg += 5;
    } else if (strcmp(arg[iarg],"rough") == 0) {
      roughflag = 1;
      roughdist = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"ramp") == 0) {
      rampflag = 1;
      iarg++;
    } else if (strcmp(arg[iarg],"rampNL1") == 0) {
      rampNL1flag = 1;
      iarg++;
    } else if (strcmp(arg[iarg],"rampNL2") == 0) {
      rampNL2flag = 1;
      iarg++;
    } else if (strcmp(arg[iarg],"rampNL3") == 0) {
      rampNL3flag = 1;
      iarg++;
    } else if (strcmp(arg[iarg],"rampNL4") == 0) {
      rampNL4flag = 1;
      iarg++;
    } else if (strcmp(arg[iarg],"rampNL5") == 0) {
      rampNL5flag = 1;
      iarg++;
    } else if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix wall/piston command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all(FLERR,"Illegal fix wall/piston command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix wall/piston command");
  }

  if (vx < 0.0 || vy < 0.0 || vz < 0.0)
    error->all(FLERR,"Illegal fix wall/piston velocity");
  if ((xloflag || xhiflag) && domain->xperiodic)
    error->all(FLERR,"Cannot use wall in periodic dimension");
  if ((yloflag || yhiflag) && domain->yperiodic)
    error->all(FLERR,"Cannot use wall in periodic dimension");
  if ((zloflag || zhiflag) && domain->zperiodic)
    error->all(FLERR,"Cannot use wall in periodic dimension");

  // setup scaling

  const double zscale = (scaleflag) ? domain->lattice->zlattice : 1.0;
  vz *= zscale;
  z0 *= zscale;
  roughdist *= zscale;

  if (rampflag || rampNL1flag || rampNL2flag || rampNL3flag ||
      rampNL4flag || rampNL5flag) {
    maxvx = vx;
    maxvy = vy;
    maxvz = vz;
  }
}

/* ---------------------------------------------------------------------- */

FixWallPiston::~FixWallPiston()
{
  delete[] gfactor2;
  delete[] gfactor1;
  delete randomt;
}

/* ---------------------------------------------------------------------- */

int FixWallPiston::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= POST_INTEGRATE;
  mask |= POST_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixWallPiston::initial_integrate(int /*vflag*/)
{
  next_reneighbor = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void FixWallPiston::post_integrate()
{
  double zlo;

  double **x = atom->x;
  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double t    = (update->ntimestep - update->beginstep) * update->dt;
  double tott = (update->endstep - update->beginstep) * update->dt;
  double tt = t * t;
  double ttt = tt * t;
  double tttt = tt * tt;
  double t0p5 = sqrt(t/tott);
  double t1p5 = t0p5*t0p5*t0p5;
  double t2p5 = t1p5*t0p5*t0p5;

  if (rampflag) {
    paccelx = maxvx / tott;
    paccely = maxvy / tott;
    paccelz = maxvz / tott;
    if (zloflag) {
      zlo = z0 + 0.5 * paccelz * tt; vz =  paccelz * t;
    }
  }
  else if (rampNL1flag) {
    paccelz = maxvz / tott;
    angfreq = MY_2PI / (0.5 * tott);

    if (zloflag) {
      zlo = z0 + paccelz * (0.5*tt + 1.0/(angfreq*angfreq) -
                            1.0/(angfreq*angfreq)*cos(angfreq*t));
      vz =  paccelz * (t + 1.0/angfreq*sin(angfreq*t));
    }
    else error->all(FLERR,
                    "NL ramp in wall/piston only implemented in zlo for now");
  }
  else if (rampNL2flag) {
    paccelz = maxvz / tott;
    angfreq = 3.0*MY_2PI / tott;

    if (zloflag) {
      zlo = z0 + paccelz * (0.5*tt + 4.0/(3.0*angfreq*angfreq)*
                            (1.0-cos(angfreq*t)) +
                            1.0/(6.0*angfreq*angfreq)*(1.0-cos(2.0*angfreq*t)));
      vz =  paccelz * (t + 4.0/(3.0*angfreq)*sin(angfreq*t) +
                       1.0/(3.0*angfreq)*sin(2.0*angfreq*t));
    }
    else error->all(FLERR,
                    "NL ramp in wall/piston only implemented in zlo for now");
  }
  else if (rampNL3flag) {
    paccelz = maxvz / tott;

    if (zloflag) {
      zlo = z0 + paccelz*tott*tott/2.5 * (t2p5 );
      vz =  paccelz * tott * (t1p5 );
    }
    else error->all(FLERR,
                    "NL ramp in wall/piston only implemented in zlo for now");
  }
  else if (rampNL4flag) {
    paccelz = maxvz / tott;

    if (zloflag) {
      zlo = z0 + paccelz/tott/3.0 * (ttt);
      vz =  paccelz / tott * (tt);
    }
    else error->all(FLERR,
                    "NL ramp in wall/piston only implemented in zlo for now");
  }
  else if (rampNL5flag) {
    paccelz = maxvz / tott;

    if (zloflag) {
      zlo = z0 + paccelz/tott/tott/4.0 * (tttt);
      vz =  paccelz / tott / tott * (ttt);
    }
    else error->all(FLERR,
                    "NL ramp in wall/piston only implemented in zlo for now");
  }
  else {
    if (zloflag) { zlo = z0 + vz * t; }
  }

  if (update->ntimestep % 1000 == 0)
    if (comm->me == 0) {
      if (screen)
        fprintf(screen,"SHOCK: step " BIGINT_FORMAT
                " t %g zpos %g vz %g az %g zlo %g\n",
                update->ntimestep, t, zlo, vz, paccelz, domain->boxlo[2]);
      if (logfile)
        fprintf(logfile,"SHOCK: step " BIGINT_FORMAT
                " t %g zpos %g vz %g az %g zlo %g\n",
                update->ntimestep, t, zlo, vz, paccelz, domain->boxlo[2]);
    }

  // VIRIAL PRESSURE CONTRIBUTION?

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      roughoff = 0.0;
      if (roughflag) {
        roughoff += roughdist*fabs((x[i][0] - domain->boxlo[0])/
                                   (domain->boxhi[0]-domain->boxlo[0])-0.5);
        roughoff += roughdist*fabs((x[i][1] - domain->boxlo[1])/
                                   (domain->boxhi[1]-domain->boxlo[1])-0.5);
      }
      if (zloflag && x[i][2] < zlo - roughoff) {
        x[i][2] = 2.0 * (zlo - roughoff) - x[i][2];
        v[i][2] = 2.0 * vz - v[i][2];
      }
    }
  }
  double **f = atom->f;
  int  *type = atom->type;

  double gamma1,gamma2;
  double tsqrt = sqrt(t_target);

  if (atom->mass) {
    if (tempflag) {
      for (int i = 1; i <= atom->ntypes; i++) {
        gfactor1[i] = -atom->mass[i] / t_period / force->ftm2v;
        gfactor2[i] = sqrt(atom->mass[i]) *
          sqrt(24.0*force->boltz/t_period/update->dt/force->mvv2e) /
          force->ftm2v;
      }
    }
    for (int i = 0; i < nlocal; i++) {
      // SET TEMP AHEAD OF PISTON
      if (tempflag && x[i][2] <= domain->boxlo[2] + t_extent ) {
        gamma1 = gfactor1[type[i]];
        gamma2 = gfactor2[type[i]] * tsqrt;
        f[i][0] += gamma1*v[i][0] + gamma2*(randomt->uniform()-0.5);
        f[i][1] += gamma1*v[i][1] + gamma2*(randomt->uniform()-0.5);
        f[i][2] += gamma1*(v[i][2]-vz) + gamma2*(randomt->uniform()-0.5);
      }
    }
  } else {
    double *rmass = atom->rmass;
    double boltz = force->boltz;
    double dt = update->dt;
    double mvv2e = force->mvv2e;
    double ftm2v = force->ftm2v;

    for (int i = 0; i < nlocal; i++) {
      // SET TEMP AHEAD OF PISTON
      if (tempflag && x[i][2] <= domain->boxlo[2] + t_extent ) {
        gamma1 = -rmass[i] / t_period / ftm2v;
        gamma2 = sqrt(rmass[i]) * sqrt(24.0*boltz/t_period/dt/mvv2e) / ftm2v;
        gamma2 *= tsqrt;
        f[i][0] += gamma1*v[i][0] + gamma2*(randomt->uniform()-0.5);
        f[i][1] += gamma1*v[i][1] + gamma2*(randomt->uniform()-0.5);
        f[i][2] += gamma1*v[i][2] + gamma2*(randomt->uniform()-0.5);
      }
    }
  }
}
