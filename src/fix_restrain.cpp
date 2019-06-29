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
   Contributing author: Craig Tenney (University of Notre Dame)
     support for bond and angle restraints by Andres Jaramillo-Botero (Caltech)
------------------------------------------------------------------------- */

#include <mpi.h>
#include <cmath>
#include <cstring>
#include "fix_restrain.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "domain.h"
#include "comm.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum{BOND,ANGLE,DIHEDRAL};

#define TOLERANCE 0.05
#define SMALL 0.001
#define DELTA 1

/* ---------------------------------------------------------------------- */

FixRestrain::FixRestrain(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  rstyle(NULL), mult(NULL), ids(NULL), kstart(NULL), kstop(NULL), target(NULL),
  cos_target(NULL), sin_target(NULL)
{
  if (narg < 4) error->all(FLERR,"Illegal fix restrain command");

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;
  vector_flag = 1;
  size_vector = 3;
  extvector = 1;
  respa_level_support = 1;
  ilevel_respa = 0;

  // parse args

  nrestrain = maxrestrain = 0;

  int iarg = 3;
  while (iarg < narg) {
    if (nrestrain == maxrestrain) {
      maxrestrain += DELTA;
      memory->grow(rstyle,maxrestrain,"restrain:rstyle");
      memory->grow(mult,maxrestrain,"restrain:mult");
      memory->grow(ids,maxrestrain,4,"restrain:ids");
      memory->grow(kstart,maxrestrain,"restrain:kstart");
      memory->grow(kstop,maxrestrain,"restrain:kstop");
      memory->grow(target,maxrestrain,"restrain:target");
      memory->grow(cos_target,maxrestrain,"restrain:cos_target");
      memory->grow(sin_target,maxrestrain,"restrain:sin_target");
    }

    if (strcmp(arg[iarg],"bond") == 0) {
      if (iarg+6 > narg) error->all(FLERR,"Illegal fix restrain command");
      rstyle[nrestrain] = BOND;
      ids[nrestrain][0] = force->inumeric(FLERR,arg[iarg+1]);
      ids[nrestrain][1] = force->inumeric(FLERR,arg[iarg+2]);
      kstart[nrestrain] = force->numeric(FLERR,arg[iarg+3]);
      kstop[nrestrain] = force->numeric(FLERR,arg[iarg+4]);
      target[nrestrain] = force->numeric(FLERR,arg[iarg+5]);
      iarg += 6;
    } else if (strcmp(arg[iarg],"angle") == 0) {
      if (iarg+7 > narg) error->all(FLERR,"Illegal fix restrain command");
      rstyle[nrestrain] = ANGLE;
      ids[nrestrain][0] = force->inumeric(FLERR,arg[iarg+1]);
      ids[nrestrain][1] = force->inumeric(FLERR,arg[iarg+2]);
      ids[nrestrain][2] = force->inumeric(FLERR,arg[iarg+3]);
      kstart[nrestrain] = force->numeric(FLERR,arg[iarg+4]);
      kstop[nrestrain] = force->numeric(FLERR,arg[iarg+5]);
      target[nrestrain] = force->numeric(FLERR,arg[iarg+6]);
      target[nrestrain] *= MY_PI / 180.0;
      iarg += 7;
    } else if (strcmp(arg[iarg],"dihedral") == 0) {
      if (iarg+8 > narg) error->all(FLERR,"Illegal fix restrain command");
      rstyle[nrestrain] = DIHEDRAL;
      mult[nrestrain]   = 1;
      ids[nrestrain][0] = force->inumeric(FLERR,arg[iarg+1]);
      ids[nrestrain][1] = force->inumeric(FLERR,arg[iarg+2]);
      ids[nrestrain][2] = force->inumeric(FLERR,arg[iarg+3]);
      ids[nrestrain][3] = force->inumeric(FLERR,arg[iarg+4]);
      kstart[nrestrain] = force->numeric(FLERR,arg[iarg+5]);
      kstop[nrestrain] = force->numeric(FLERR,arg[iarg+6]);
      target[nrestrain] = force->numeric(FLERR,arg[iarg+7]);
      target[nrestrain] *= MY_PI / 180.0;
      cos_target[nrestrain] = cos(target[nrestrain]);
      sin_target[nrestrain] = sin(target[nrestrain]);
      iarg += 8;
      if ((iarg < narg) && (strcmp("mult",arg[iarg]) == 0)) {
        if (iarg+1 > narg) error->all(FLERR,"Illegal fix restrain command");
        mult[nrestrain] = force->inumeric(FLERR,arg[iarg+1]);
        if (mult[nrestrain] < 0)
          error->all(FLERR,"Illegal fix restrain command");
        iarg += 2;
      }
    } else error->all(FLERR,"Illegal fix restrain command");

    nrestrain++;
  }

  // require atom map to lookup atom IDs

  if (atom->map_style == 0)
    error->all(FLERR,"Fix restrain requires an atom map, see atom_modify");
}

/* ---------------------------------------------------------------------- */

FixRestrain::~FixRestrain()
{
  memory->destroy(rstyle);
  memory->destroy(mult);
  memory->destroy(ids);
  memory->destroy(kstart);
  memory->destroy(kstop);
  memory->destroy(target);
  memory->destroy(cos_target);
  memory->destroy(sin_target);
}

/* ---------------------------------------------------------------------- */

int FixRestrain::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixRestrain::init()
{
  if (strstr(update->integrate_style,"respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixRestrain::setup(int vflag)
{
  if (strcmp(update->integrate_style,"verlet") == 0)
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixRestrain::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixRestrain::post_force(int /*vflag*/)
{
  energy = 0.0;

  ebond = 0.0;
  eangle = 0.0;
  edihed = 0.0;

  for (int m = 0; m < nrestrain; m++)
    if (rstyle[m] == BOND) restrain_bond(m);
    else if (rstyle[m] == ANGLE) restrain_angle(m);
    else if (rstyle[m] == DIHEDRAL) restrain_dihedral(m);
}

/* ---------------------------------------------------------------------- */

void FixRestrain::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixRestrain::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   apply harmonic bond restraints
---------------------------------------------------------------------- */

void FixRestrain::restrain_bond(int m)
{
  int i1,i2;
  double delx,dely,delz,fbond;
  double rsq,r,dr,rk;

  double **x = atom->x;
  double **f = atom->f;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;
  double k = kstart[m] + delta * (kstop[m] - kstart[m]);

  i1 = atom->map(ids[m][0]);
  i2 = atom->map(ids[m][1]);

  // newton_bond on: only processor owning i2 computes restraint
  // newton_bond off: only processors owning either of i1,i2 computes restraint

  if (newton_bond) {
    if (i2 == -1 || i2 >= nlocal) return;
    if (i1 == -1) {
      char str[128];
      sprintf(str,
              "Restrain atoms %d %d missing on proc %d at step " BIGINT_FORMAT,
              ids[m][0],ids[m][1],
              comm->me,update->ntimestep);
      error->one(FLERR,str);
    }
  } else {
    if ((i1 == -1 || i1 >= nlocal) && (i2 == -1 || i2 >= nlocal)) return;
    if (i1 == -1 || i2 == -1)  {
      char str[128];
      sprintf(str,
              "Restrain atoms %d %d missing on proc %d at step " BIGINT_FORMAT,
              ids[m][0],ids[m][1],
              comm->me,update->ntimestep);
      error->one(FLERR,str);
    }
  }

  delx = x[i1][0] - x[i2][0];
  dely = x[i1][1] - x[i2][1];
  delz = x[i1][2] - x[i2][2];
  domain->minimum_image(delx,dely,delz);

  rsq = delx*delx + dely*dely + delz*delz;
  r = sqrt(rsq);
  dr = r - target[m];
  rk = k * dr;

  // force & energy

  if (r > 0.0) fbond = -2.0*rk/r;
  else fbond = 0.0;

  ebond += rk*dr;
  energy += rk*dr;

  // apply force to each of 2 atoms

  if (newton_bond || i1 < nlocal) {
    f[i1][0] += delx*fbond;
    f[i1][1] += dely*fbond;
    f[i1][2] += delz*fbond;
  }

  if (newton_bond || i2 < nlocal) {
    f[i2][0] -= delx*fbond;
    f[i2][1] -= dely*fbond;
    f[i2][2] -= delz*fbond;
  }
}

/* ----------------------------------------------------------------------
   apply harmonic angle restraints
---------------------------------------------------------------------- */

void FixRestrain::restrain_angle(int m)
{
  int i1,i2,i3;
  double delx1,dely1,delz1,delx2,dely2,delz2;
  double f1[3],f3[3];
  double dtheta,tk;
  double rsq1,rsq2,r1,r2,c,s,a,a11,a12,a22;

  double **x = atom->x;
  double **f = atom->f;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;
  double k = kstart[m] + delta * (kstop[m] - kstart[m]);

  i1 = atom->map(ids[m][0]);
  i2 = atom->map(ids[m][1]);
  i3 = atom->map(ids[m][2]);

  // newton_bond on: only processor owning i2 computes restraint
  // newton_bond off: only processors owning any of i1-i3 computes restraint

  if (newton_bond) {
    if (i2 == -1 || i2 >= nlocal) return;
    if (i1 == -1 || i3 == -1) {
      char str[128];
      sprintf(str,
              "Restrain atoms %d %d %d missing on proc %d at step "
              BIGINT_FORMAT,
              ids[m][0],ids[m][1],ids[m][2],
              comm->me,update->ntimestep);
      error->one(FLERR,str);
    }
  } else {
    if ((i1 == -1 || i1 >= nlocal) && (i2 == -1 || i2 >= nlocal) &&
        (i3 == -1 || i3 >= nlocal)) return;
    if (i1 == -1 || i2 == -1 || i3 == -1) {
      char str[128];
      sprintf(str,
              "Restrain atoms %d %d %d missing on proc %d at step "
              BIGINT_FORMAT,
              ids[m][0],ids[m][1],ids[m][2],
              comm->me,update->ntimestep);
      error->one(FLERR,str);
    }
  }

  // 1st bond

  delx1 = x[i1][0] - x[i2][0];
  dely1 = x[i1][1] - x[i2][1];
  delz1 = x[i1][2] - x[i2][2];
  domain->minimum_image(delx1,dely1,delz1);

  rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1;
  r1 = sqrt(rsq1);

  // 2nd bond

  delx2 = x[i3][0] - x[i2][0];
  dely2 = x[i3][1] - x[i2][1];
  delz2 = x[i3][2] - x[i2][2];
  domain->minimum_image(delx2,dely2,delz2);

  rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2;
  r2 = sqrt(rsq2);

  // angle (cos and sin)

  c = delx1*delx2 + dely1*dely2 + delz1*delz2;
  c /= r1*r2;

  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;

  s = sqrt(1.0 - c*c);
  if (s < SMALL) s = SMALL;
  s = 1.0/s;

  // force & energy

  dtheta = acos(c) - target[m];
  tk = k * dtheta;

  eangle += tk*dtheta;
  energy += tk*dtheta;

  a = -2.0 * tk * s;
  a11 = a*c / rsq1;
  a12 = -a / (r1*r2);
  a22 = a*c / rsq2;

  f1[0] = a11*delx1 + a12*delx2;
  f1[1] = a11*dely1 + a12*dely2;
  f1[2] = a11*delz1 + a12*delz2;
  f3[0] = a22*delx2 + a12*delx1;
  f3[1] = a22*dely2 + a12*dely1;
  f3[2] = a22*delz2 + a12*delz1;

  // apply force to each of 3 atoms

  if (newton_bond || i1 < nlocal) {
    f[i1][0] += f1[0];
    f[i1][1] += f1[1];
    f[i1][2] += f1[2];
  }

  if (newton_bond || i2 < nlocal) {
    f[i2][0] -= f1[0] + f3[0];
    f[i2][1] -= f1[1] + f3[1];
    f[i2][2] -= f1[2] + f3[2];
  }

  if (newton_bond || i3 < nlocal) {
    f[i3][0] += f3[0];
    f[i3][1] += f3[1];
    f[i3][2] += f3[2];
  }
}

/* ----------------------------------------------------------------------
   apply dihedral restraints
   adopted from dihedral_charmm
---------------------------------------------------------------------- */

void FixRestrain::restrain_dihedral(int m)
{
  int i1,i2,i3,i4,i;
  double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z,vb2xm,vb2ym,vb2zm;
  double f1[3],f2[3],f3[3],f4[3];
  double ax,ay,az,bx,by,bz,rasq,rbsq,rgsq,rg,rginv,ra2inv,rb2inv,rabinv;
  double df,df1,ddf1,fg,hg,fga,hgb,gaa,gbb;
  double dtfx,dtfy,dtfz,dtgx,dtgy,dtgz,dthx,dthy,dthz;
  double c,s,p,sx2,sy2,sz2;

  double **x = atom->x;
  double **f = atom->f;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;
  double k = kstart[m] + delta * (kstop[m] - kstart[m]);

  i1 = atom->map(ids[m][0]);
  i2 = atom->map(ids[m][1]);
  i3 = atom->map(ids[m][2]);
  i4 = atom->map(ids[m][3]);

  // newton_bond on: only processor owning i2 computes restraint
  // newton_bond off: only processors owning any of i1-i4 computes restraint

  if (newton_bond) {
    if (i2 == -1 || i2 >= nlocal) return;
    if (i1 == -1 || i3 == -1 || i4 == -1) {
      char str[128];
      sprintf(str,
              "Restrain atoms %d %d %d %d missing on proc %d at step "
              BIGINT_FORMAT,
              ids[m][0],ids[m][1],ids[m][2],ids[m][3],
              comm->me,update->ntimestep);
      error->one(FLERR,str);
    }
  } else {
    if ((i1 == -1 || i1 >= nlocal) && (i2 == -1 || i2 >= nlocal) &&
        (i3 == -1 || i3 >= nlocal) && (i4 == -1 || i3 >= nlocal)) return;
    if (i1 == -1 || i2 == -1 || i3 == -1 || i4 == -1) {
      char str[128];
      sprintf(str,
              "Restrain atoms %d %d %d %d missing on proc %d at step "
              BIGINT_FORMAT,
              ids[m][0],ids[m][1],ids[m][2],ids[m][3],
              comm->me,update->ntimestep);
      error->one(FLERR,str);
    }
  }

  // 1st bond

  vb1x = x[i1][0] - x[i2][0];
  vb1y = x[i1][1] - x[i2][1];
  vb1z = x[i1][2] - x[i2][2];
  domain->minimum_image(vb1x,vb1y,vb1z);

  // 2nd bond

  vb2x = x[i3][0] - x[i2][0];
  vb2y = x[i3][1] - x[i2][1];
  vb2z = x[i3][2] - x[i2][2];
  domain->minimum_image(vb2x,vb2y,vb2z);

  vb2xm = -vb2x;
  vb2ym = -vb2y;
  vb2zm = -vb2z;
  domain->minimum_image(vb2xm,vb2ym,vb2zm);

  // 3rd bond

  vb3x = x[i4][0] - x[i3][0];
  vb3y = x[i4][1] - x[i3][1];
  vb3z = x[i4][2] - x[i3][2];
  domain->minimum_image(vb3x,vb3y,vb3z);

  ax = vb1y*vb2zm - vb1z*vb2ym;
  ay = vb1z*vb2xm - vb1x*vb2zm;
  az = vb1x*vb2ym - vb1y*vb2xm;
  bx = vb3y*vb2zm - vb3z*vb2ym;
  by = vb3z*vb2xm - vb3x*vb2zm;
  bz = vb3x*vb2ym - vb3y*vb2xm;

  rasq = ax*ax + ay*ay + az*az;
  rbsq = bx*bx + by*by + bz*bz;
  rgsq = vb2xm*vb2xm + vb2ym*vb2ym + vb2zm*vb2zm;
  rg = sqrt(rgsq);

  rginv = ra2inv = rb2inv = 0.0;
  if (rg > 0) rginv = 1.0/rg;
  if (rasq > 0) ra2inv = 1.0/rasq;
  if (rbsq > 0) rb2inv = 1.0/rbsq;
  rabinv = sqrt(ra2inv*rb2inv);

  c = (ax*bx + ay*by + az*bz)*rabinv;
  s = rg*rabinv*(ax*vb3x + ay*vb3y + az*vb3z);

  // error check

  if (c > 1.0 + TOLERANCE || c < (-1.0 - TOLERANCE)) {
    int me;
    MPI_Comm_rank(world,&me);
    if (screen) {
      char str[128];
      sprintf(str,"Restrain problem: %d " BIGINT_FORMAT " "
              TAGINT_FORMAT " " TAGINT_FORMAT " "
              TAGINT_FORMAT " " TAGINT_FORMAT,
              me,update->ntimestep,
              atom->tag[i1],atom->tag[i2],atom->tag[i3],atom->tag[i4]);
      error->warning(FLERR,str);
      fprintf(screen,"  1st atom: %d %g %g %g\n",
              me,x[i1][0],x[i1][1],x[i1][2]);
      fprintf(screen,"  2nd atom: %d %g %g %g\n",
              me,x[i2][0],x[i2][1],x[i2][2]);
      fprintf(screen,"  3rd atom: %d %g %g %g\n",
              me,x[i3][0],x[i3][1],x[i3][2]);
      fprintf(screen,"  4th atom: %d %g %g %g\n",
              me,x[i4][0],x[i4][1],x[i4][2]);
    }
  }

  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;

  p = 1.0;
  df1 = 0.0;

  for (i = 0; i < mult[m]; i++) {
    ddf1 = p*c - df1*s;
    df1 = p*s + df1*c;
    p = ddf1;
  }

  p = p*cos_target[m] + df1*sin_target[m];
  df1 = df1*cos_target[m] - ddf1*sin_target[m];
  df1 *= -mult[m];
  p += 1.0;

  edihed += k * p;
  energy += k * p;

  fg = vb1x*vb2xm + vb1y*vb2ym + vb1z*vb2zm;
  hg = vb3x*vb2xm + vb3y*vb2ym + vb3z*vb2zm;
  fga = fg*ra2inv*rginv;
  hgb = hg*rb2inv*rginv;
  gaa = -ra2inv*rg;
  gbb = rb2inv*rg;

  dtfx = gaa*ax;
  dtfy = gaa*ay;
  dtfz = gaa*az;
  dtgx = fga*ax - hgb*bx;
  dtgy = fga*ay - hgb*by;
  dtgz = fga*az - hgb*bz;
  dthx = gbb*bx;
  dthy = gbb*by;
  dthz = gbb*bz;

  df = -k * df1;

  sx2 = df*dtgx;
  sy2 = df*dtgy;
  sz2 = df*dtgz;

  f1[0] = df*dtfx;
  f1[1] = df*dtfy;
  f1[2] = df*dtfz;

  f2[0] = sx2 - f1[0];
  f2[1] = sy2 - f1[1];
  f2[2] = sz2 - f1[2];

  f4[0] = df*dthx;
  f4[1] = df*dthy;
  f4[2] = df*dthz;

  f3[0] = -sx2 - f4[0];
  f3[1] = -sy2 - f4[1];
  f3[2] = -sz2 - f4[2];

  // apply force to each of 4 atoms

  if (newton_bond || i1 < nlocal) {
    f[i1][0] += f1[0];
    f[i1][1] += f1[1];
    f[i1][2] += f1[2];
  }

  if (newton_bond || i2 < nlocal) {
    f[i2][0] += f2[0];
    f[i2][1] += f2[1];
    f[i2][2] += f2[2];
  }

  if (newton_bond || i3 < nlocal) {
    f[i3][0] += f3[0];
    f[i3][1] += f3[1];
    f[i3][2] += f3[2];
  }

  if (newton_bond || i4 < nlocal) {
    f[i4][0] += f4[0];
    f[i4][1] += f4[1];
    f[i4][2] += f4[2];
  }
}

/* ----------------------------------------------------------------------
   potential energy of added force
------------------------------------------------------------------------- */

double FixRestrain::compute_scalar()
{
  MPI_Allreduce(&energy,&energy_all,1,MPI_DOUBLE,MPI_SUM,world);
  return energy_all;
}

/* ----------------------------------------------------------------------
  return individual energy contributions
------------------------------------------------------------------------- */

double FixRestrain::compute_vector(int n)
{
  if (n == 0) {
    MPI_Allreduce(&ebond,&ebond_all,1,MPI_DOUBLE,MPI_SUM,world);
    return ebond_all;
  } else if (n == 1) {
    MPI_Allreduce(&eangle,&eangle_all,1,MPI_DOUBLE,MPI_SUM,world);
    return eangle_all;
  } else if (n == 2) {
    MPI_Allreduce(&edihed,&edihed_all,1,MPI_DOUBLE,MPI_SUM,world);
    return edihed_all;
  } else {
    return 0.0;
  }
}
