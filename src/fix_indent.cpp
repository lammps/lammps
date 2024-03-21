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
   Contributing author: Ravi Agrawal (Northwestern U)
------------------------------------------------------------------------- */

#include "fix_indent.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "input.h"
#include "lattice.h"
#include "math_extra.h"
#include "modify.h"
#include "respa.h"
#include "update.h"
#include "variable.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

enum { NONE, SPHERE, CYLINDER, PLANE, CONE };
enum { INSIDE, OUTSIDE };

/* ---------------------------------------------------------------------- */

FixIndent::FixIndent(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), xstr(nullptr), ystr(nullptr), zstr(nullptr), rstr(nullptr), pstr(nullptr),
    rlostr(nullptr), rhistr(nullptr), lostr(nullptr), histr(nullptr)
{
  if (narg < 4) utils::missing_cmd_args(FLERR, "fix indent", error);

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 3;
  energy_global_flag = 1;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;
  respa_level_support = 1;
  ilevel_respa = 0;

  k = utils::numeric(FLERR, arg[3], false, lmp);
  if (k < 0.0) error->all(FLERR, "Illegal fix indent force constant: {}", k);
  k3 = k / 3.0;

  // read geometry of indenter and optional args

  int iarg = geometry(narg - 4, &arg[4]) + 4;
  options(narg - iarg, &arg[iarg]);

  // setup scaling

  const double xscale{scaleflag ? domain->lattice->xlattice : 1.0};
  const double yscale{scaleflag ? domain->lattice->ylattice : 1.0};
  const double zscale{scaleflag ? domain->lattice->zlattice : 1.0};

  // apply scaling factors to geometry

  if (istyle == SPHERE || istyle == CYLINDER) {
    if (!xstr) xvalue *= xscale;
    if (!ystr) yvalue *= yscale;
    if (!zstr) zvalue *= zscale;
    if (!rstr) rvalue *= xscale;

  } else if (istyle == CONE) {
    if (!xstr) xvalue *= xscale;
    if (!ystr) yvalue *= yscale;
    if (!zstr) zvalue *= zscale;

    double scaling_factor = 1.0;
    switch (cdim) {
      case 0:
        scaling_factor = xscale;
        break;
      case 1:
        scaling_factor = yscale;
        break;
      case 2:
        scaling_factor = zscale;
        break;
    }

    if (!rlostr) rlovalue *= scaling_factor;
    if (!rhistr) rhivalue *= scaling_factor;
    if (!lostr) lovalue *= scaling_factor;
    if (!histr) hivalue *= scaling_factor;

  } else if (istyle == PLANE) {
    if (cdim == 0 && !pstr)
      pvalue *= xscale;
    else if (cdim == 1 && !pstr)
      pvalue *= yscale;
    else if (cdim == 2 && !pstr)
      pvalue *= zscale;

  } else
    error->all(FLERR, "Unknown fix indent keyword: {}", istyle);

  varflag = 0;
  if (xstr || ystr || zstr || rstr || pstr || rlostr || rhistr || lostr || histr) varflag = 1;

  indenter_flag = 0;
  indenter[0] = indenter[1] = indenter[2] = indenter[3] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixIndent::~FixIndent()
{
  delete[] xstr;
  delete[] ystr;
  delete[] zstr;
  delete[] rstr;
  delete[] pstr;
  delete[] rlostr;
  delete[] rhistr;
  delete[] lostr;
  delete[] histr;
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
    if (xvar < 0) error->all(FLERR, "Variable {} for fix indent does not exist", xstr);
    if (!input->variable->equalstyle(xvar))
      error->all(FLERR, "Variable {} for fix indent is invalid style", xstr);
  }
  if (ystr) {
    yvar = input->variable->find(ystr);
    if (yvar < 0) error->all(FLERR, "Variable {} for fix indent does not exist", ystr);
    if (!input->variable->equalstyle(yvar))
      error->all(FLERR, "Variable {} for fix indent is invalid style", ystr);
  }
  if (zstr) {
    zvar = input->variable->find(zstr);
    if (zvar < 0) error->all(FLERR, "Variable {} for fix indent does not exist", zstr);
    if (!input->variable->equalstyle(zvar))
      error->all(FLERR, "Variable {} for fix indent is invalid style", zstr);
  }
  if (rstr) {
    rvar = input->variable->find(rstr);
    if (rvar < 0) error->all(FLERR, "Variable {} for fix indent does not exist", rstr);
    if (!input->variable->equalstyle(rvar))
      error->all(FLERR, "Variable {} for fix indent is invalid style", rstr);
  }
  if (pstr) {
    pvar = input->variable->find(pstr);
    if (pvar < 0) error->all(FLERR, "Variable {} for fix indent does not exist", pstr);
    if (!input->variable->equalstyle(pvar))
      error->all(FLERR, "Variable {} for fix indent is invalid style", pstr);
  }
  if (rlostr) {
    rlovar = input->variable->find(rlostr);
    if (rlovar < 0) error->all(FLERR, "Variable {} for fix indent does not exist", rlostr);
    if (!input->variable->equalstyle(rlovar))
      error->all(FLERR, "Variable {} for fix indent is invalid style", rlostr);
  }
  if (rhistr) {
    rhivar = input->variable->find(rhistr);
    if (rhivar < 0) error->all(FLERR, "Variable {} for fix indent does not exist", rhistr);
    if (!input->variable->equalstyle(rhivar))
      error->all(FLERR, "Variable {} for fix indent is invalid style", rhistr);
  }
  if (lostr) {
    lovar = input->variable->find(lostr);
    if (lovar < 0) error->all(FLERR, "Variable {} for fix indent does not exist", lostr);
    if (!input->variable->equalstyle(lovar))
      error->all(FLERR, "Variable {} for fix indent is invalid style", lostr);
  }
  if (histr) {
    hivar = input->variable->find(histr);
    if (hivar < 0) error->all(FLERR, "Variable {} for fix indent does not exist", histr);
    if (!input->variable->equalstyle(hivar))
      error->all(FLERR, "Variable {} for fix indent is invalid style", histr);
  }

  if (utils::strmatch(update->integrate_style, "^respa")) {
    ilevel_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels - 1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level, ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixIndent::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style, "^verlet"))
    post_force(vflag);
  else {
    (dynamic_cast<Respa *>(update->integrate))->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag, ilevel_respa, 0);
    (dynamic_cast<Respa *>(update->integrate))->copy_f_flevel(ilevel_respa);
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

  // ctr = current indenter centerz

  double ctr[3] = {xvalue, yvalue, zvalue};
  if (xstr) ctr[0] = input->variable->compute_equal(xvar);
  if (ystr) ctr[1] = input->variable->compute_equal(yvar);
  if (zstr) ctr[2] = input->variable->compute_equal(zvar);

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double delx, dely, delz, r, dr, fmag, fx, fy, fz;

  // spherical indenter

  if (istyle == SPHERE) {

    // remap indenter center into periodic box

    domain->remap(ctr);

    double radius = rstr ? input->variable->compute_equal(rvar) : rvalue;
    if (radius < 0.0) error->all(FLERR, "Illegal fix indent sphere radius: {}", radius);

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        delx = x[i][0] - ctr[0];
        dely = x[i][1] - ctr[1];
        delz = x[i][2] - ctr[2];
        domain->minimum_image(delx, dely, delz);
        r = sqrt(delx * delx + dely * dely + delz * delz);
        if (side == OUTSIDE) {
          dr = r - radius;
          fmag = k * dr * dr;
        } else {
          dr = radius - r;
          fmag = -k * dr * dr;
        }
        if (dr >= 0.0) continue;
        fx = delx * fmag / r;
        fy = dely * fmag / r;
        fz = delz * fmag / r;
        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;
        indenter[0] -= k3 * dr * dr * dr;
        indenter[1] -= fx;
        indenter[2] -= fy;
        indenter[3] -= fz;
      }

    // cylindrical indenter

  } else if (istyle == CYLINDER) {

    // ctr = current indenter axis
    // remap into periodic box
    // 3rd coord is just near box for remap(), since isn't used

    ctr[cdim] = domain->boxlo[cdim];
    domain->remap(ctr);

    double radius{rstr ? input->variable->compute_equal(rvar) : rvalue};
    if (radius < 0.0) error->all(FLERR, "Illegal fix indent cylinder radius: {}", radius);

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        double del[3] = {x[i][0] - ctr[0], x[i][1] - ctr[1], x[i][2] - ctr[2]};
        del[cdim] = 0;
        domain->minimum_image(del[0], del[1], del[2]);
        r = sqrt(del[0] * del[0] + del[1] * del[1] + del[2] * del[2]);
        if (side == OUTSIDE) {
          dr = r - radius;
          fmag = k * dr * dr;
        } else {
          dr = radius - r;
          fmag = -k * dr * dr;
        }
        if (dr >= 0.0) continue;
        fx = del[0] * fmag / r;
        fy = del[1] * fmag / r;
        fz = del[2] * fmag / r;
        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;
        indenter[0] -= k3 * dr * dr * dr;
        indenter[1] -= fx;
        indenter[2] -= fy;
        indenter[3] -= fz;
      }

    // conical indenter

  } else if (istyle == CONE) {

    double radiuslo{rlostr ? input->variable->compute_equal(rlovar) : rlovalue};
    if (radiuslo < 0.0) error->all(FLERR, "Illegal fix indent cone lower radius: {}", radiuslo);
    double radiushi{rhistr ? input->variable->compute_equal(rhivar) : rhivalue};
    if (radiushi < 0.0) error->all(FLERR, "Illegal fix indent cone high radius: {}", radiushi);

    double initial_lo{lostr ? input->variable->compute_equal(lovar) : lovalue};
    double initial_hi{histr ? input->variable->compute_equal(hivar) : hivalue};

    ctr[cdim] = 0.5 * (initial_hi + initial_lo);

    domain->remap(ctr);

    double hi = ctr[cdim] + 0.5 * (initial_hi - initial_lo);
    double lo = ctr[cdim] - 0.5 * (initial_hi - initial_lo);

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {

        delx = x[i][0] - ctr[0];
        dely = x[i][1] - ctr[1];
        delz = x[i][2] - ctr[2];
        domain->minimum_image(delx, dely, delz);

        double x0[3] = {delx + ctr[0], dely + ctr[1], delz + ctr[2]};
        r = sqrt(delx * delx + dely * dely + delz * delz);

        // check if particle is inside or outside the cone

        bool point_inside_cone = PointInsideCone(cdim, ctr, lo, hi, radiuslo, radiushi, x0);

        if (side == INSIDE && point_inside_cone) continue;
        if (side == OUTSIDE && !point_inside_cone) continue;

        // find the distance between the point and the cone

        if (point_inside_cone) {
          DistanceInteriorPoint(cdim, ctr, lo, hi, radiuslo, radiushi, x0[0], x0[1], x0[2]);
        } else {
          DistanceExteriorPoint(cdim, ctr, lo, hi, radiuslo, radiushi, x0[0], x0[1], x0[2]);
        }

        // compute the force from the center of the cone
        // this is different from how it is done in fix wall/region

        dr = sqrt(x0[0] * x0[0] + x0[1] * x0[1] + x0[2] * x0[2]);

        int force_sign = {point_inside_cone ? 1 : -1};
        fmag = force_sign * k * dr * dr;

        fx = delx * fmag / r;
        fy = dely * fmag / r;
        fz = delz * fmag / r;
        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;
        indenter[0] -= k3 * dr * dr * dr;
        indenter[1] -= fx;
        indenter[2] -= fy;
        indenter[3] -= fz;
      }
    }

    // planar indenter

  } else {

    // plane = current plane position

    double plane{pstr ? input->variable->compute_equal(pvar) : pvalue};

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dr = planeside * (plane - x[i][cdim]);
        if (dr >= 0.0) continue;
        fmag = -planeside * k * dr * dr;
        f[i][cdim] += fmag;
        indenter[0] -= k3 * dr * dr * dr;
        indenter[cdim + 1] -= fmag;
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
    MPI_Allreduce(indenter, indenter_all, 4, MPI_DOUBLE, MPI_SUM, world);
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
    MPI_Allreduce(indenter, indenter_all, 4, MPI_DOUBLE, MPI_SUM, world);
    indenter_flag = 1;
  }
  return indenter_all[n + 1];
}

/* ----------------------------------------------------------------------
   parse input args for geometry of indenter
------------------------------------------------------------------------- */

int FixIndent::geometry(int narg, char **arg)
{
  if (narg < 0) utils::missing_cmd_args(FLERR, "fix indent", error);

  istyle = NONE;
  xstr = ystr = zstr = rstr = pstr = nullptr;
  xvalue = yvalue = zvalue = rvalue = pvalue = 0.0;

  // sphere

  if (strcmp(arg[0], "sphere") == 0) {
    if (istyle != NONE) error->all(FLERR, "Fix indent requires a single geometry keyword");
    if (5 > narg) utils::missing_cmd_args(FLERR, "fix indent sphere", error);

    if (utils::strmatch(arg[1], "^v_")) {
      xstr = utils::strdup(arg[1] + 2);
    } else
      xvalue = utils::numeric(FLERR, arg[1], false, lmp);
    if (utils::strmatch(arg[2], "^v_")) {
      ystr = utils::strdup(arg[2] + 2);
    } else
      yvalue = utils::numeric(FLERR, arg[2], false, lmp);
    if (utils::strmatch(arg[3], "^v_")) {
      zstr = utils::strdup(arg[3] + 2);
    } else
      zvalue = utils::numeric(FLERR, arg[3], false, lmp);
    if (utils::strmatch(arg[4], "^v_")) {
      rstr = utils::strdup(arg[4] + 2);
    } else
      rvalue = utils::numeric(FLERR, arg[4], false, lmp);

    istyle = SPHERE;
    return 5;
  }

  // cylinder

  if (strcmp(arg[0], "cylinder") == 0) {
    if (istyle != NONE) error->all(FLERR, "Fix indent requires a single geometry keyword");
    if (5 > narg) utils::missing_cmd_args(FLERR, "fix indent cylinder", error);

    if (strcmp(arg[1], "x") == 0) {
      cdim = 0;
      if (utils::strmatch(arg[2], "^v_")) {
        ystr = utils::strdup(arg[2] + 2);
      } else
        yvalue = utils::numeric(FLERR, arg[2], false, lmp);
      if (utils::strmatch(arg[3], "^v_")) {
        zstr = utils::strdup(arg[3] + 2);
      } else
        zvalue = utils::numeric(FLERR, arg[3], false, lmp);
    } else if (strcmp(arg[1], "y") == 0) {
      cdim = 1;
      if (utils::strmatch(arg[2], "^v_")) {
        xstr = utils::strdup(arg[2] + 2);
      } else
        xvalue = utils::numeric(FLERR, arg[2], false, lmp);
      if (utils::strmatch(arg[3], "^v_")) {
        zstr = utils::strdup(arg[3] + 2);
      } else
        zvalue = utils::numeric(FLERR, arg[3], false, lmp);
    } else if (strcmp(arg[1], "z") == 0) {
      cdim = 2;
      if (utils::strmatch(arg[2], "^v_")) {
        xstr = utils::strdup(arg[2] + 2);
      } else
        xvalue = utils::numeric(FLERR, arg[2], false, lmp);
      if (utils::strmatch(arg[3], "^v_")) {
        ystr = utils::strdup(arg[3] + 2);
      } else
        yvalue = utils::numeric(FLERR, arg[3], false, lmp);
    } else
      error->all(FLERR, "Unknown fix indent cylinder argument: {}", arg[1]);

    if (utils::strmatch(arg[4], "^v_")) {
      rstr = utils::strdup(arg[4] + 2);
    } else
      rvalue = utils::numeric(FLERR, arg[4], false, lmp);

    istyle = CYLINDER;
    return 5;
  }

  // cone

  if (strcmp(arg[0], "cone") == 0) {
    if (istyle != NONE) error->all(FLERR, "Fix indent requires a single geometry keyword");
    if (8 > narg) utils::missing_cmd_args(FLERR, "fix indent cone", error);

    if (strcmp(arg[1], "x") == 0) {
      cdim = 0;
      if (utils::strmatch(arg[2], "^v_")) {
        ystr = utils::strdup(arg[2] + 2);
      } else
        yvalue = utils::numeric(FLERR, arg[2], false, lmp);
      if (utils::strmatch(arg[3], "^v_")) {
        zstr = utils::strdup(arg[3] + 2);
      } else
        zvalue = utils::numeric(FLERR, arg[3], false, lmp);

    } else if (strcmp(arg[1], "y") == 0) {
      cdim = 1;
      if (utils::strmatch(arg[2], "^v_")) {
        xstr = utils::strdup(arg[2] + 2);
      } else
        xvalue = utils::numeric(FLERR, arg[2], false, lmp);
      if (utils::strmatch(arg[3], "^v_")) {
        zstr = utils::strdup(arg[3] + 2);
      } else
        zvalue = utils::numeric(FLERR, arg[3], false, lmp);

    } else if (strcmp(arg[1], "z") == 0) {
      cdim = 2;
      if (utils::strmatch(arg[2], "^v_")) {
        xstr = utils::strdup(arg[2] + 2);
      } else
        xvalue = utils::numeric(FLERR, arg[2], false, lmp);
      if (utils::strmatch(arg[3], "^v_")) {
        ystr = utils::strdup(arg[3] + 2);
      } else
        yvalue = utils::numeric(FLERR, arg[3], false, lmp);

    } else
      error->all(FLERR, "Unknown fix indent cone argument: {}", arg[1]);

    if (utils::strmatch(arg[4], "^v_")) {
      rlostr = utils::strdup(arg[4] + 2);
    } else
      rlovalue = utils::numeric(FLERR, arg[4], false, lmp);
    if (utils::strmatch(arg[5], "^v_")) {
      rhistr = utils::strdup(arg[5] + 2);
    } else
      rhivalue = utils::numeric(FLERR, arg[5], false, lmp);
    if (utils::strmatch(arg[6], "^v_")) {
      lostr = utils::strdup(arg[6] + 2);
    } else
      lovalue = utils::numeric(FLERR, arg[6], false, lmp);
    if (utils::strmatch(arg[7], "^v_")) {
      histr = utils::strdup(arg[7] + 2);
    } else
      hivalue = utils::numeric(FLERR, arg[7], false, lmp);

    istyle = CONE;
    return 8;
  }

  // plane

  if (strcmp(arg[0], "plane") == 0) {
    if (istyle != NONE) error->all(FLERR, "Fix indent requires a single geometry keyword");
    if (4 > narg) utils::missing_cmd_args(FLERR, "fix indent plane", error);
    if (strcmp(arg[1], "x") == 0)
      cdim = 0;
    else if (strcmp(arg[1], "y") == 0)
      cdim = 1;
    else if (strcmp(arg[1], "z") == 0)
      cdim = 2;
    else
      error->all(FLERR, "Unknown fix indent plane argument: {}", arg[1]);

    if (utils::strmatch(arg[2], "^v_")) {
      pstr = utils::strdup(arg[2] + 2);
    } else
      pvalue = utils::numeric(FLERR, arg[2], false, lmp);

    if (strcmp(arg[3], "lo") == 0)
      planeside = -1;
    else if (strcmp(arg[3], "hi") == 0)
      planeside = 1;
    else
      error->all(FLERR, "Unknown fix indent plane argument: {}", arg[3]);
    istyle = PLANE;
    return 4;
  }

  // invalid istyle arg

  error->all(FLERR, "Unknown fix indent argument: {}", arg[0]);

  return 0;
}

/* ----------------------------------------------------------------------
   parse optional input args
------------------------------------------------------------------------- */

void FixIndent::options(int narg, char **arg)
{
  scaleflag = 1;
  side = OUTSIDE;

  int iarg = 0;

  while (iarg < narg) {
    if (strcmp(arg[iarg], "units") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix indent units", error);
      if (strcmp(arg[iarg + 1], "box") == 0)
        scaleflag = 0;
      else if (strcmp(arg[iarg + 1], "lattice") == 0)
        scaleflag = 1;
      else
        error->all(FLERR, "Unknown fix indent units argument: {}", arg[iarg + 1]);
      iarg += 2;

    } else if (strcmp(arg[iarg], "side") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix indent side", error);
      if (strcmp(arg[iarg + 1], "in") == 0)
        side = INSIDE;
      else if (strcmp(arg[iarg + 1], "out") == 0)
        side = OUTSIDE;
      else
        error->all(FLERR, "Unknown fix indent side argument: {}", arg[iarg + 1]);
      iarg += 2;

    } else
      error->all(FLERR, "Unknown fix indent argument: {}", arg[iarg]);
  }
}

/* ----------------------------------------------------------------------
   determines if a point is inside (true) or outside (false) of a cone
------------------------------------------------------------------------- */

bool FixIndent::PointInsideCone(int dir, double *center, double lo, double hi, double rlo,
                                double rhi, double *x)
{
  if ((x[dir] > hi) || (x[dir] < lo)) return false;

  double del[3] = {x[0] - center[0], x[1] - center[1], x[2] - center[2]};
  del[dir] = 0.0;

  double dist = sqrt(del[0] * del[0] + del[1] * del[1] + del[2] * del[2]);
  double currentradius = rlo + (x[dir] - lo) * (rhi - rlo) / (hi - lo);

  if (dist > currentradius) return false;

  return true;
}

/* ----------------------------------------------------------------------
   distance between an exterior point and a cone
------------------------------------------------------------------------- */

void FixIndent::DistanceExteriorPoint(int dir, double *center, double lo, double hi, double rlo,
                                      double rhi, double &x, double &y, double &z)
{
  double xp[3], nearest[3], corner1[3], corner2[3];
  double point[3] = {x, y, z};
  double del[3] = {x - center[0], y - center[1], z - center[2]};

  del[dir] = 0.0;
  double r = sqrt(del[0] * del[0] + del[1] * del[1] + del[2] * del[2]);

  corner1[0] = center[0] + del[0] * rlo / r;
  corner1[1] = center[1] + del[1] * rlo / r;
  corner1[2] = center[2] + del[2] * rlo / r;
  corner1[dir] = lo;

  corner2[0] = center[0] + del[0] * rhi / r;
  corner2[1] = center[1] + del[1] * rhi / r;
  corner2[2] = center[2] + del[2] * rhi / r;
  corner2[dir] = hi;

  double corner3[3] = {center[0], center[1], center[2]};
  corner3[dir] = lo;

  double corner4[3] = {center[0], center[1], center[2]};
  corner4[dir] = hi;

  // initialize distance to a big number

  double distsq = 1.0e20;

  // check the first triangle

  point_on_line_segment(corner1, corner2, point, xp);
  distsq = closest(point, xp, nearest, distsq);

  // check the second triangle

  point_on_line_segment(corner1, corner3, point, xp);
  distsq = closest(point, xp, nearest, distsq);

  // check the third triangle

  point_on_line_segment(corner2, corner4, point, xp);
  distsq = closest(point, xp, nearest, distsq);

  x -= nearest[0];
  y -= nearest[1];
  z -= nearest[2];

  return;
}

/* ----------------------------------------------------------------------
   distance between an interior point and a cone
------------------------------------------------------------------------- */

void FixIndent::DistanceInteriorPoint(int dir, double *center, double lo, double hi, double rlo,
                                      double rhi, double &x, double &y, double &z)
{
  double r, dist_disk, dist_surf;
  double surflo[3], surfhi[3], xs[3];
  double initial_point[3] = {x, y, z};
  double point[3] = {0.0, 0.0, 0.0};

  // initial check with the two disks

  if ((initial_point[dir] - lo) < (hi - initial_point[dir])) {
    dist_disk = (initial_point[dir] - lo) * (initial_point[dir] - lo);
    point[dir] = initial_point[dir] - lo;
  } else {
    dist_disk = (hi - initial_point[dir]) * (hi - initial_point[dir]);
    point[dir] = initial_point[dir] - hi;
  }

  // check with the points in the conical surface

  double del[3] = {x - center[0], y - center[1], z - center[2]};
  del[dir] = 0.0;
  r = sqrt(del[0] * del[0] + del[1] * del[1] + del[2] * del[2]);

  surflo[0] = center[0] + del[0] * rlo / r;
  surflo[1] = center[1] + del[1] * rlo / r;
  surflo[2] = center[2] + del[2] * rlo / r;
  surflo[dir] = lo;

  surfhi[0] = center[0] + del[0] * rhi / r;
  surfhi[1] = center[1] + del[1] * rhi / r;
  surfhi[2] = center[2] + del[2] * rhi / r;
  surfhi[dir] = hi;

  point_on_line_segment(surflo, surfhi, initial_point, xs);

  double dx[3] = {initial_point[0] - xs[0], initial_point[1] - xs[1], initial_point[2] - xs[2]};
  dist_surf = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
  if (dist_surf < dist_disk) {
    x = dx[0];
    y = dx[1];
    z = dx[2];
  } else {
    x = point[0];
    y = point[1];
    z = point[2];
  }

  return;
}

/* ----------------------------------------------------------------------
   helper function extracted from region.cpp
------------------------------------------------------------------------- */

void FixIndent::point_on_line_segment(double *a, double *b, double *c, double *d)
{
  double ba[3], ca[3];

  MathExtra::sub3(b, a, ba);
  MathExtra::sub3(c, a, ca);
  double t = MathExtra::dot3(ca, ba) / MathExtra::dot3(ba, ba);
  if (t <= 0.0) {
    d[0] = a[0];
    d[1] = a[1];
    d[2] = a[2];
  } else if (t >= 1.0) {
    d[0] = b[0];
    d[1] = b[1];
    d[2] = b[2];
  } else {
    d[0] = a[0] + t * ba[0];
    d[1] = a[1] + t * ba[1];
    d[2] = a[2] + t * ba[2];
  }
}

/* ----------------------------------------------------------------------
   helper function extracted from region_cone.cpp
------------------------------------------------------------------------- */

double FixIndent::closest(double *x, double *near, double *nearest, double dsq)
{
  double dx = x[0] - near[0];
  double dy = x[1] - near[1];
  double dz = x[2] - near[2];
  double rsq = dx * dx + dy * dy + dz * dz;
  if (rsq >= dsq) return dsq;

  nearest[0] = near[0];
  nearest[1] = near[1];
  nearest[2] = near[2];
  return rsq;
}
