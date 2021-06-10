// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Ravi Agrawal (Northwestern U)
------------------------------------------------------------------------- */

#include "fix_indent.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "input.h"
#include "lattice.h"
#include "modify.h"
#include "respa.h"
#include "update.h"
#include "variable.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,SPHERE,CYLINDER,PLANE};
enum{INSIDE,OUTSIDE};

/* ---------------------------------------------------------------------- */

FixIndent::FixIndent(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  xstr(nullptr), ystr(nullptr), zstr(nullptr), rstr(nullptr), pstr(nullptr)
{
  if (narg < 4) error->all(FLERR,"Illegal fix indent command");

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 3;
  energy_global_flag = 1;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;
  respa_level_support = 1;
  ilevel_respa = 0;

  k = utils::numeric(FLERR,arg[3],false,lmp);
  k3 = k/3.0;

  // read options from end of input line

  options(narg-4,&arg[4]);

  // setup scaling

  double xscale,yscale,zscale;
  if (scaleflag) {
    xscale = domain->lattice->xlattice;
    yscale = domain->lattice->ylattice;
    zscale = domain->lattice->zlattice;
  }
  else xscale = yscale = zscale = 1.0;

  // apply scaling factors to geometry

  if (istyle == SPHERE || istyle == CYLINDER) {
    if (!xstr) xvalue *= xscale;
    if (!ystr) yvalue *= yscale;
    if (!zstr) zvalue *= zscale;
    if (!rstr) rvalue *= xscale;
  } else if (istyle == PLANE) {
    if (cdim == 0 && !pstr) pvalue *= xscale;
    else if (cdim == 1 && !pstr) pvalue *= yscale;
    else if (cdim == 2 && !pstr) pvalue *= zscale;
  } else error->all(FLERR,"Illegal fix indent command");

  varflag = 0;
  if (xstr || ystr || zstr || rstr || pstr) varflag = 1;

  indenter_flag = 0;
  indenter[0] = indenter[1] = indenter[2] = indenter[3] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixIndent::~FixIndent()
{
  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
  delete [] rstr;
  delete [] pstr;
}

/* ---------------------------------------------------------------------- */

int FixIndent::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixIndent::init()
{
  if (xstr) {
    xvar = input->variable->find(xstr);
    if (xvar < 0)
      error->all(FLERR,"Variable name for fix indent does not exist");
    if (!input->variable->equalstyle(xvar))
      error->all(FLERR,"Variable for fix indent is invalid style");
  }
  if (ystr) {
    yvar = input->variable->find(ystr);
    if (yvar < 0)
      error->all(FLERR,"Variable name for fix indent does not exist");
    if (!input->variable->equalstyle(yvar))
      error->all(FLERR,"Variable for fix indent is not equal style");
  }
  if (zstr) {
    zvar = input->variable->find(zstr);
    if (zvar < 0)
      error->all(FLERR,"Variable name for fix indent does not exist");
    if (!input->variable->equalstyle(zvar))
      error->all(FLERR,"Variable for fix indent is not equal style");
  }
  if (rstr) {
    rvar = input->variable->find(rstr);
    if (rvar < 0)
      error->all(FLERR,"Variable name for fix indent does not exist");
    if (!input->variable->equalstyle(rvar))
      error->all(FLERR,"Variable for fix indent is not equal style");
  }
  if (pstr) {
    pvar = input->variable->find(pstr);
    if (pvar < 0)
      error->all(FLERR,"Variable name for fix indent does not exist");
    if (!input->variable->equalstyle(pvar))
      error->all(FLERR,"Variable for fix indent is not equal style");
  }

  if (utils::strmatch(update->integrate_style,"^respa")) {
    ilevel_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixIndent::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style,"^verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixIndent::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixIndent::post_force(int /*vflag*/)
{
  // indenter values, 0 = energy, 1-3 = force components
  // wrap variable evaluations with clear/add

  if (varflag) modify->clearstep_compute();

  indenter_flag = 0;
  indenter[0] = indenter[1] = indenter[2] = indenter[3] = 0.0;

  // spherical indenter

  if (istyle == SPHERE) {

    // ctr = current indenter center
    // remap into periodic box

    double ctr[3];
    if (xstr) ctr[0] = input->variable->compute_equal(xvar);
    else ctr[0] = xvalue;
    if (ystr) ctr[1] = input->variable->compute_equal(yvar);
    else ctr[1] = yvalue;
    if (zstr) ctr[2] = input->variable->compute_equal(zvar);
    else ctr[2] = zvalue;
    domain->remap(ctr);

    double radius;
    if (rstr) radius = input->variable->compute_equal(rvar);
    else radius = rvalue;

    double **x = atom->x;
    double **f = atom->f;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    double delx,dely,delz,r,dr,fmag,fx,fy,fz;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        delx = x[i][0] - ctr[0];
        dely = x[i][1] - ctr[1];
        delz = x[i][2] - ctr[2];
        domain->minimum_image(delx,dely,delz);
        r = sqrt(delx*delx + dely*dely + delz*delz);
        if (side == OUTSIDE) {
          dr = r - radius;
          fmag = k*dr*dr;
        } else {
          dr = radius - r;
          fmag = -k*dr*dr;
        }
        if (dr >= 0.0) continue;
        fx = delx*fmag/r;
        fy = dely*fmag/r;
        fz = delz*fmag/r;
        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;
        indenter[0] -= k3 * dr*dr*dr;
        indenter[1] -= fx;
        indenter[2] -= fy;
        indenter[3] -= fz;
      }

  // cylindrical indenter

  } else if (istyle == CYLINDER) {

    // ctr = current indenter axis
    // remap into periodic box
    // 3rd coord is just near box for remap(), since isn't used

    double ctr[3];
    if (cdim == 0) {
      ctr[0] = domain->boxlo[0];
      if (ystr) ctr[1] = input->variable->compute_equal(yvar);
      else ctr[1] = yvalue;
      if (zstr) ctr[2] = input->variable->compute_equal(zvar);
      else ctr[2] = zvalue;
    } else if (cdim == 1) {
      if (xstr) ctr[0] = input->variable->compute_equal(xvar);
      else ctr[0] = xvalue;
      ctr[1] = domain->boxlo[1];
      if (zstr) ctr[2] = input->variable->compute_equal(zvar);
      else ctr[2] = zvalue;
    } else {
      if (xstr) ctr[0] = input->variable->compute_equal(xvar);
      else ctr[0] = xvalue;
      if (ystr) ctr[1] = input->variable->compute_equal(yvar);
      else ctr[1] = yvalue;
      ctr[2] = domain->boxlo[2];
    }
    domain->remap(ctr);

    double radius;
    if (rstr) radius = input->variable->compute_equal(rvar);
    else radius = rvalue;

    double **x = atom->x;
    double **f = atom->f;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    double delx,dely,delz,r,dr,fmag,fx,fy,fz;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (cdim == 0) {
          delx = 0;
          dely = x[i][1] - ctr[1];
          delz = x[i][2] - ctr[2];
        } else if (cdim == 1) {
          delx = x[i][0] - ctr[0];
          dely = 0;
          delz = x[i][2] - ctr[2];
        } else {
          delx = x[i][0] - ctr[0];
          dely = x[i][1] - ctr[1];
          delz = 0;
        }
        domain->minimum_image(delx,dely,delz);
        r = sqrt(delx*delx + dely*dely + delz*delz);
        if (side == OUTSIDE) {
          dr = r - radius;
          fmag = k*dr*dr;
        } else {
          dr = radius - r;
          fmag = -k*dr*dr;
        }
        if (dr >= 0.0) continue;
        fx = delx*fmag/r;
        fy = dely*fmag/r;
        fz = delz*fmag/r;
        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;
        indenter[0] -= k3 * dr*dr*dr;
        indenter[1] -= fx;
        indenter[2] -= fy;
        indenter[3] -= fz;
      }

  // planar indenter

  } else {

    // plane = current plane position

    double plane;
    if (pstr) plane = input->variable->compute_equal(pvar);
    else plane = pvalue;

    double **x = atom->x;
    double **f = atom->f;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    double dr,fatom;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dr = planeside * (plane - x[i][cdim]);
        if (dr >= 0.0) continue;
        fatom = -planeside * k*dr*dr;
        f[i][cdim] += fatom;
        indenter[0] -= k3 * dr*dr*dr;
        indenter[cdim+1] -= fatom;
      }
  }

  if (varflag) modify->addstep_compute(update->ntimestep + 1);
}

/* ---------------------------------------------------------------------- */

void FixIndent::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixIndent::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   energy of indenter interaction
------------------------------------------------------------------------- */

double FixIndent::compute_scalar()
{
  // only sum across procs one time

  if (indenter_flag == 0) {
    MPI_Allreduce(indenter,indenter_all,4,MPI_DOUBLE,MPI_SUM,world);
    indenter_flag = 1;
  }
  return indenter_all[0];
}

/* ----------------------------------------------------------------------
   components of force on indenter
------------------------------------------------------------------------- */

double FixIndent::compute_vector(int n)
{
  // only sum across procs one time

  if (indenter_flag == 0) {
    MPI_Allreduce(indenter,indenter_all,4,MPI_DOUBLE,MPI_SUM,world);
    indenter_flag = 1;
  }
  return indenter_all[n+1];
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of input line
------------------------------------------------------------------------- */

void FixIndent::options(int narg, char **arg)
{
  if (narg < 0) error->all(FLERR,"Illegal fix indent command");

  istyle = NONE;
  xstr = ystr = zstr = rstr = pstr = nullptr;
  xvalue = yvalue = zvalue = rvalue = pvalue = 0.0;
  scaleflag = 1;
  side = OUTSIDE;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"sphere") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal fix indent command");

      if (utils::strmatch(arg[iarg+1],"^v_")) {
        xstr = utils::strdup(arg[iarg+1]+2);
      } else xvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (utils::strmatch(arg[iarg+2],"^v_")) {
        ystr = utils::strdup(arg[iarg+2]+2);
      } else yvalue = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      if (utils::strmatch(arg[iarg+3],"^v_")) {
        zstr = utils::strdup(arg[iarg+3]+2);
      } else zvalue = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      if (utils::strmatch(arg[iarg+4],"^v_")) {
        rstr = utils::strdup(arg[iarg+4]+2);
      } else rvalue = utils::numeric(FLERR,arg[iarg+4],false,lmp);

      istyle = SPHERE;
      iarg += 5;

    } else if (strcmp(arg[iarg],"cylinder") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal fix indent command");

      if (strcmp(arg[iarg+1],"x") == 0) {
        cdim = 0;
        if (utils::strmatch(arg[iarg+2],"^v_")) {
          ystr = utils::strdup(arg[iarg+2]+2);
        } else yvalue = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        if (utils::strmatch(arg[iarg+3],"^v_")) {
          zstr = utils::strdup(arg[iarg+3]+2);
        } else zvalue = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      } else if (strcmp(arg[iarg+1],"y") == 0) {
        cdim = 1;
        if (utils::strmatch(arg[iarg+2],"^v_")) {
          xstr = utils::strdup(arg[iarg+2]+2);
        } else xvalue = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        if (utils::strmatch(arg[iarg+3],"^v_")) {
          zstr = utils::strdup(arg[iarg+3]+2);
        } else zvalue = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      } else if (strcmp(arg[iarg+1],"z") == 0) {
        cdim = 2;
        if (utils::strmatch(arg[iarg+2],"^v_")) {
          xstr = utils::strdup(arg[iarg+2]+2);
        } else xvalue = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        if (utils::strmatch(arg[iarg+3],"^v_")) {
          ystr = utils::strdup(arg[iarg+3]+2);
        } else yvalue = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      } else error->all(FLERR,"Illegal fix indent command");

      if (utils::strmatch(arg[iarg+4],"^v_")) {
        rstr = utils::strdup(arg[iarg+4]+2);
      } else rvalue = utils::numeric(FLERR,arg[iarg+4],false,lmp);

      istyle = CYLINDER;
      iarg += 5;

    } else if (strcmp(arg[iarg],"plane") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix indent command");
      if (strcmp(arg[iarg+1],"x") == 0) cdim = 0;
      else if (strcmp(arg[iarg+1],"y") == 0) cdim = 1;
      else if (strcmp(arg[iarg+1],"z") == 0) cdim = 2;
      else error->all(FLERR,"Illegal fix indent command");

      if (utils::strmatch(arg[iarg+2],"^v_")) {
        pstr = utils::strdup(arg[iarg+2]+2);
      } else pvalue = utils::numeric(FLERR,arg[iarg+2],false,lmp);

      if (strcmp(arg[iarg+3],"lo") == 0) planeside = -1;
      else if (strcmp(arg[iarg+3],"hi") == 0) planeside = 1;
      else error->all(FLERR,"Illegal fix indent command");
      istyle = PLANE;
      iarg += 4;

    } else if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix indent command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all(FLERR,"Illegal fix indent command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"side") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix indent command");
      if (strcmp(arg[iarg+1],"in") == 0) side = INSIDE;
      else if (strcmp(arg[iarg+1],"out") == 0) side = OUTSIDE;
      else error->all(FLERR,"Illegal fix indent command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix indent command");
  }
}
