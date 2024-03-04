// clang-format off
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
   Contributing author: Axel Kohlmeyer (UPenn)
   based on fix spring by: Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include "fix_smd.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "respa.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

enum { SMD_NONE=0,
       SMD_TETHER=1<<0, SMD_COUPLE=1<<1,
       SMD_CVEL=1<<2, SMD_CFOR=1<<3,
       SMD_AUTOX=1<<4, SMD_AUTOY=1<<5, SMD_AUTOZ=1<<6};

static constexpr double SMALL = 0.001;

/* ---------------------------------------------------------------------- */

FixSMD::FixSMD(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  styleflag = SMD_NONE;
  k_smd = f_smd = v_smd = -1.0;
  xflag = yflag = zflag = 1;
  xc = yc = zc = 0.0;
  xn = yn = zn = 1.0;
  pmf = r_old = r_now = r0 = 0.0;

  restart_global = 1;
  vector_flag = 1;
  size_vector = 7;
  global_freq = 1;
  extvector = 1;
  respa_level_support = 1;
  ilevel_respa = 0;
  virial_global_flag = virial_peratom_flag = 1;

  int argoffs=3;
  if (strcmp(arg[argoffs],"cvel") == 0) {
    if (narg < argoffs+3) error->all(FLERR,"Illegal fix smd command");
    styleflag |= SMD_CVEL;
    k_smd = utils::numeric(FLERR,arg[argoffs+1],false,lmp);
    v_smd = utils::numeric(FLERR,arg[argoffs+2],false,lmp); // to be multiplied by update->dt when used.
    argoffs += 3;
  } else if (strcmp(arg[argoffs],"cfor") == 0) {
    if (narg < argoffs+2) error->all(FLERR,"Illegal fix smd command");
    styleflag |= SMD_CFOR;
    f_smd = utils::numeric(FLERR,arg[argoffs+1],false,lmp);
    argoffs += 2;
  } else error->all(FLERR,"Illegal fix smd command");

  if (strcmp(arg[argoffs],"tether") == 0) {
    if (narg < argoffs+5) error->all(FLERR,"Illegal fix smd command");
    styleflag |= SMD_TETHER;
    if (strcmp(arg[argoffs+1],"NULL") == 0) xflag = 0;
    else xc = utils::numeric(FLERR,arg[argoffs+1],false,lmp);
    if (strcmp(arg[argoffs+2],"NULL") == 0) yflag = 0;
    else yc = utils::numeric(FLERR,arg[argoffs+2],false,lmp);
    if (strcmp(arg[argoffs+3],"NULL") == 0) zflag = 0;
    else zc = utils::numeric(FLERR,arg[argoffs+3],false,lmp);
    r0 = utils::numeric(FLERR,arg[argoffs+4],false,lmp);
    if (r0 < 0) error->all(FLERR,"R0 < 0 for fix smd command");
    argoffs += 5;
  } else if (strcmp(arg[argoffs],"couple") == 0) {
    if (narg < argoffs+6) error->all(FLERR,"Illegal fix smd command");
    styleflag |= SMD_COUPLE;
    igroup2 = group->find(arg[argoffs+1]);
    if (igroup2 == -1)
      error->all(FLERR,"Could not find fix smd couple group ID");
    if (igroup2 == igroup)
      error->all(FLERR,"Two groups cannot be the same in fix smd couple");
    group2bit = group->bitmask[igroup2];

    if (strcmp(arg[argoffs+2],"NULL") == 0) xflag = 0;
    else if (strcmp(arg[argoffs+2],"auto") == 0) styleflag |= SMD_AUTOX;
    else xc = utils::numeric(FLERR,arg[argoffs+2],false,lmp);
    if (strcmp(arg[argoffs+3],"NULL") == 0) yflag = 0;
    else if (strcmp(arg[argoffs+3],"auto") == 0) styleflag |= SMD_AUTOY;
    else yc = utils::numeric(FLERR,arg[argoffs+3],false,lmp);
    if (strcmp(arg[argoffs+4],"NULL") == 0) zflag = 0;
    else if (strcmp(arg[argoffs+4],"auto") == 0) styleflag |= SMD_AUTOZ;
    else zc = utils::numeric(FLERR,arg[argoffs+4],false,lmp);

    r0 = utils::numeric(FLERR,arg[argoffs+5],false,lmp);
    if (r0 < 0) error->all(FLERR,"R0 < 0 for fix smd command");
    argoffs +=6;
  } else error->all(FLERR,"Illegal fix smd command");

  force_flag = 0;
  ftotal[0] = ftotal[1] = ftotal[2] = 0.0;
}

/* ---------------------------------------------------------------------- */

int FixSMD::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSMD::init()
{
  double xcm[3], xcm2[3];
  masstotal = group->mass(igroup);
  group->xcm(igroup,masstotal,xcm);

  double dx,dy,dz;
  if (styleflag & SMD_TETHER) {
    dx = xc - xcm[0];
    dy = yc - xcm[1];
    dz = zc - xcm[2];
  } else {     /* SMD_COUPLE */
    masstotal2 = group->mass(igroup2);
    group->xcm(igroup2,masstotal2,xcm2);
    if (styleflag & SMD_AUTOX) dx = xcm2[0] - xcm[0];
    else dx = xc;
    if (styleflag & SMD_AUTOY) dy = xcm2[1] - xcm[1];
    else dy = yc;
    if (styleflag & SMD_AUTOZ) dz = xcm2[2] - xcm[2];
    else dz = zc;
  }

  if (!xflag) dx = 0.0;
  if (!yflag) dy = 0.0;
  if (!zflag) dz = 0.0;
  r_old = sqrt(dx*dx + dy*dy + dz*dz);
  if (r_old > SMALL) {
    xn = dx/r_old;
    yn = dy/r_old;
    zn = dz/r_old;
  }

  if (utils::strmatch(update->integrate_style,"^respa")) {
    ilevel_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixSMD::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style,"^verlet"))
    post_force(vflag);
  else {
    (dynamic_cast<Respa *>(update->integrate))->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    (dynamic_cast<Respa *>(update->integrate))->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixSMD::post_force(int vflag)
{
  // virial setup

  v_init(vflag);

  if (styleflag & SMD_TETHER) smd_tether();
  else smd_couple();

  if (styleflag & SMD_CVEL) {
    if (utils::strmatch(update->integrate_style,"^verlet"))
      r_old += v_smd * update->dt;
    else
      r_old += v_smd * (dynamic_cast<Respa *>(update->integrate))->step[ilevel_respa];
  }
}

/* ---------------------------------------------------------------------- */

void FixSMD::smd_tether()
{
  double xcm[3];
  group->xcm(igroup,masstotal,xcm);

  double dt = update->dt;
  if (utils::strmatch(update->integrate_style,"^respa"))
    dt = (dynamic_cast<Respa *>(update->integrate))->step[ilevel_respa];

  // fx,fy,fz = components of k * (r-r0)

  double dx,dy,dz,fx,fy,fz,r,dr;

  dx = xcm[0] - xc;
  dy = xcm[1] - yc;
  dz = xcm[2] - zc;
  r_now = sqrt(dx*dx + dy*dy + dz*dz);

  if (!xflag) dx = 0.0;
  if (!yflag) dy = 0.0;
  if (!zflag) dz = 0.0;
  r = sqrt(dx*dx + dy*dy + dz*dz);
  if (styleflag & SMD_CVEL) {
    if (r > SMALL) {
      dr = r - r0 - r_old;
      fx = k_smd*dx*dr/r;
      fy = k_smd*dy*dr/r;
      fz = k_smd*dz*dr/r;
      pmf += (fx*xn + fy*yn + fz*zn) * v_smd * dt;
    } else {
      fx = 0.0;
      fy = 0.0;
      fz = 0.0;
    }
  } else {
    r_old = r;
    fx = f_smd*dx/r;
    fy = f_smd*dy/r;
    fz = f_smd*dz/r;
  }

  // apply restoring force to atoms in group
  // f = -k*(r-r0)*mass/masstotal

  double **x = atom->x;
  double **f = atom->f;
  imageint *image = atom->image;
  int *mask = atom->mask;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  double massfrac;
  double unwrap[3],v[6];
  int nlocal = atom->nlocal;

  ftotal[0] = ftotal[1] = ftotal[2] = 0.0;
  force_flag = 0;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        massfrac = rmass[i]/masstotal;
        f[i][0] -= fx*massfrac;
        f[i][1] -= fy*massfrac;
        f[i][2] -= fz*massfrac;
        ftotal[0] -= fx*massfrac;
        ftotal[1] -= fy*massfrac;
        ftotal[2] -= fz*massfrac;
        if (evflag) {
          domain->unmap(x[i],image[i],unwrap);
          v[0] = -fx*massfrac*unwrap[0];
          v[1] = -fy*massfrac*unwrap[1];
          v[2] = -fz*massfrac*unwrap[2];
          v[3] = -fx*massfrac*unwrap[1];
          v[4] = -fx*massfrac*unwrap[2];
          v[5] = -fy*massfrac*unwrap[2];
          v_tally(i,v);
        }
      }
  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        massfrac = mass[type[i]]/masstotal;
        f[i][0] -= fx*massfrac;
        f[i][1] -= fy*massfrac;
        f[i][2] -= fz*massfrac;
        ftotal[0] -= fx*massfrac;
        ftotal[1] -= fy*massfrac;
        ftotal[2] -= fz*massfrac;
        if (evflag) {
          domain->unmap(x[i],image[i],unwrap);
          v[0] = -fx*massfrac*unwrap[0];
          v[1] = -fy*massfrac*unwrap[1];
          v[2] = -fz*massfrac*unwrap[2];
          v[3] = -fx*massfrac*unwrap[1];
          v[4] = -fx*massfrac*unwrap[2];
          v[5] = -fy*massfrac*unwrap[2];
          v_tally(i,v);
        }
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixSMD::smd_couple()
{
  double xcm[3],xcm2[3];
  group->xcm(igroup,masstotal,xcm);
  group->xcm(igroup2,masstotal2,xcm2);

  double dt = update->dt;
  if (utils::strmatch(update->integrate_style,"^respa"))
    dt = (dynamic_cast<Respa *>(update->integrate))->step[ilevel_respa];

  // renormalize direction of spring
  double dx,dy,dz,r,dr;
  if (styleflag & SMD_AUTOX) dx = xcm2[0] - xcm[0];
  else dx = xn*r_old;
  if (styleflag & SMD_AUTOY) dy = xcm2[1] - xcm[1];
  else dy = yn*r_old;
  if (styleflag & SMD_AUTOZ) dz = xcm2[2] - xcm[2];
  else dz = zn*r_old;
  if (!xflag) dx = 0.0;
  if (!yflag) dy = 0.0;
  if (!zflag) dz = 0.0;
  r = sqrt(dx*dx + dy*dy + dz*dz);
  if (r > SMALL) {
    xn = dx/r; yn = dy/r; zn = dz/r;
  }

  double fx,fy,fz;
  if (styleflag & SMD_CVEL) {
    dx = xcm2[0] - xcm[0];
    dy = xcm2[1] - xcm[1];
    dz = xcm2[2] - xcm[2];
    r_now = sqrt(dx*dx + dy*dy + dz*dz);

    dx -= xn*r_old;
    dy -= yn*r_old;
    dz -= zn*r_old;

    if (!xflag) dx = 0.0;
    if (!yflag) dy = 0.0;
    if (!zflag) dz = 0.0;
    r = sqrt(dx*dx + dy*dy + dz*dz);
    dr = r - r0;

    if (r > SMALL) {
      double fsign;
      fsign  = (v_smd<0.0) ? -1.0 : 1.0;

      fx = k_smd*dx*dr/r;
      fy = k_smd*dy*dr/r;
      fz = k_smd*dz*dr/r;
      pmf += (fx*xn + fy*yn + fz*zn) * fsign * v_smd * dt;
    } else {
      fx = 0.0;
      fy = 0.0;
      fz = 0.0;
    }
  } else {
    dx = xcm2[0] - xcm[0];
    dy = xcm2[1] - xcm[1];
    dz = xcm2[2] - xcm[2];
    r_now = sqrt(dx*dx + dy*dy + dz*dz);
    r_old = r;

    fx = f_smd*xn;
    fy = f_smd*yn;
    fz = f_smd*zn;
  }

  // apply restoring force to atoms in group
  // f = -k*(r-r0)*mass/masstotal

  double **f = atom->f;
  int *mask = atom->mask;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  ftotal[0] = ftotal[1] = ftotal[2] = 0.0;
  force_flag = 0;

  double massfrac;
  if (rmass) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        massfrac = rmass[i]/masstotal;
        f[i][0] += fx*massfrac;
        f[i][1] += fy*massfrac;
        f[i][2] += fz*massfrac;
        ftotal[0] += fx*massfrac;
        ftotal[1] += fy*massfrac;
        ftotal[2] += fz*massfrac;
      }
      if (mask[i] & group2bit) {
        massfrac = rmass[i]/masstotal2;
        f[i][0] -= fx*massfrac;
        f[i][1] -= fy*massfrac;
        f[i][2] -= fz*massfrac;
      }
    }
  } else {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        massfrac = mass[type[i]]/masstotal;
        f[i][0] += fx*massfrac;
        f[i][1] += fy*massfrac;
        f[i][2] += fz*massfrac;
        ftotal[0] += fx*massfrac;
        ftotal[1] += fy*massfrac;
        ftotal[2] += fz*massfrac;
      }
      if (mask[i] & group2bit) {
        massfrac = mass[type[i]]/masstotal2;
        f[i][0] -= fx*massfrac;
        f[i][1] -= fy*massfrac;
        f[i][2] -= fz*massfrac;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixSMD::write_restart(FILE *fp)
{
  static constexpr int RESTART_ITEMS = 5;
  double buf[RESTART_ITEMS], fsign;

  if (comm->me == 0) {
    // make sure we project the force into the direction of the pulling.
    fsign  = (v_smd<0.0) ? -1.0 : 1.0;
    buf[0] = r_old;
    buf[1] = xn*fsign;
    buf[2] = yn*fsign;
    buf[3] = zn*fsign;
    buf[4] = pmf;
    int size = RESTART_ITEMS*sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(&buf[0],sizeof(double),RESTART_ITEMS,fp);
  }
}

/* ---------------------------------------------------------------------- */

void FixSMD::restart(char *buf)
{
  auto list = (double *)buf;
  r_old = list[0];
  xn=list[1];
  yn=list[2];
  zn=list[3];
  pmf=list[4];
}

/* ---------------------------------------------------------------------- */

void FixSMD::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ----------------------------------------------------------------------
   return components of total smd force on fix group
------------------------------------------------------------------------- */

double FixSMD::compute_vector(int n)
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(ftotal,ftotal_all,3,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
    if (styleflag & SMD_CVEL) {
      ftotal_all[3]=ftotal_all[0]*xn+ftotal_all[1]*yn+ftotal_all[2]*zn;
      ftotal_all[4]=r_old;
    } else {
      ftotal_all[3]=f_smd;
      ftotal_all[4]=r_old;
    }
    ftotal_all[5]=r_now;
    ftotal_all[6]=pmf;
  }
  return ftotal_all[n];
}
