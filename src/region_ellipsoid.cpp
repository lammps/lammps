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

#include "region_ellipsoid.h"

#include "domain.h"
#include "error.h"
#include "input.h"
#include "variable.h"

#include <cmath>
#include <limits>

using namespace LAMMPS_NS;

enum { CONSTANT, VARIABLE };

static double GetRoot2D(double r0, double z0, double z1, double g);
static double GetRoot3D(double r0, double r1, double z0, double z1, double z2, double g);

static double DistancePointEllipse(double e0, double e1, double y0, double y1, double &x0,
                                   double &x1);
static double DistancePointEllipsoid(double e0, double e1, double e2, double y0, double y1,
                                     double y2, double &x0, double &x1, double &x2);

static constexpr int maxIterations =
    std::numeric_limits<double>::digits - std::numeric_limits<double>::min_exponent;
static constexpr double EPSILON = std::numeric_limits<double>::epsilon() * 2.0;

/* ---------------------------------------------------------------------- */

RegEllipsoid::RegEllipsoid(LAMMPS *lmp, int narg, char **arg) :
    Region(lmp, narg, arg), xstr(nullptr), ystr(nullptr), zstr(nullptr), astr(nullptr),
    bstr(nullptr), cstr(nullptr)
{
  options(narg - 8, &arg[8]);

  if (utils::strmatch(arg[2], "^v_")) {
    xstr = utils::strdup(arg[2] + 2);
    xc = 0.0;
    xstyle = VARIABLE;
    varshape = 1;
  } else {
    xc = xscale * utils::numeric(FLERR, arg[2], false, lmp);
    xstyle = CONSTANT;
  }

  if (utils::strmatch(arg[3], "^v_")) {
    ystr = utils::strdup(arg[3] + 2);
    yc = 0.0;
    ystyle = VARIABLE;
    varshape = 1;
  } else {
    yc = yscale * utils::numeric(FLERR, arg[3], false, lmp);
    ystyle = CONSTANT;
  }

  if (utils::strmatch(arg[4], "^v_")) {
    zstr = utils::strdup(arg[4] + 2);
    zc = 0.0;
    zstyle = VARIABLE;
    varshape = 1;
  } else {
    zc = zscale * utils::numeric(FLERR, arg[4], false, lmp);
    zstyle = CONSTANT;
  }

  if (utils::strmatch(arg[5], "^v_")) {
    astr = utils::strdup(arg[5] + 2);
    a = 0.0;
    astyle = VARIABLE;
    varshape = 1;
  } else {
    a = xscale * utils::numeric(FLERR, arg[5], false, lmp);
    astyle = CONSTANT;
  }

  if (utils::strmatch(arg[6], "^v_")) {
    bstr = utils::strdup(arg[6] + 2);
    b = 0.0;
    bstyle = VARIABLE;
    varshape = 1;
  } else {
    b = yscale * utils::numeric(FLERR, arg[6], false, lmp);
    bstyle = CONSTANT;
  }

  if (utils::strmatch(arg[7], "^v_")) {
    cstr = utils::strdup(arg[7] + 2);
    c = 0.0;
    cstyle = VARIABLE;
    varshape = 1;
  } else {
    c = zscale * utils::numeric(FLERR, arg[7], false, lmp);
    cstyle = CONSTANT;
  }

  if (varshape) {
    variable_check();
    shape_update();
  }

  // error check

  if (a < 0.0) error->all(FLERR, "Illegal region ellipsoid a: {}", a);
  if (b < 0.0) error->all(FLERR, "Illegal region ellipsoid b: {}", b);
  if (c < 0.0) error->all(FLERR, "Illegal region ellipsoid c: {}", c);

  // extent of ellipsoid
  // for variable axes, uses initial axes and origin for variable center

  if (interior) {
    bboxflag = 1;
    extent_xlo = xc - a;
    extent_xhi = xc + a;
    extent_ylo = yc - b;
    extent_yhi = yc + b;
    extent_zlo = zc - c;
    extent_zhi = zc + c;
  } else
    bboxflag = 0;

  cmax = 1;
  contact = new Contact[cmax];
  tmax = 1;
}

/* ---------------------------------------------------------------------- */

RegEllipsoid::~RegEllipsoid()
{
  delete[] xstr;
  delete[] ystr;
  delete[] zstr;
  delete[] astr;
  delete[] bstr;
  delete[] cstr;
  delete[] contact;
}

/* ---------------------------------------------------------------------- */

void RegEllipsoid::init()
{
  Region::init();
  if (varshape) variable_check();
}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside or on surface
   inside = 0 if x,y,z is outside and not on surface
------------------------------------------------------------------------- */

int RegEllipsoid::inside(double x, double y, double z)
{
  if (domain->dimension == 3) {
    double delx = b * c * (x - xc);
    double dely = a * c * (y - yc);
    double delz = a * b * (z - zc);
    double r = delx * delx + dely * dely + delz * delz;
    double rc = a * a * b * b * c * c;
    if (r <= rc) return 1;
  } else {
    double delx = b * (x - xc);
    double dely = a * (y - yc);
    double r = delx * delx + dely * dely;
    double rc = a * a * b * b;
    if (r <= rc) return 1;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   one contact if 0 <= x < cutoff from inner surface of ellipsoid
   no contact if outside (possible if called from union/intersect)
   delxyz = vector from nearest point on ellipsoid to x
   special case: no contact if x is at center of ellipsoid
------------------------------------------------------------------------- */

int RegEllipsoid::surface_interior(double *x, double cutoff)
{

  if (domain->dimension == 3) {
    double delx = b * c * (x[0] - xc);
    double dely = a * c * (x[1] - yc);
    double delz = a * b * (x[2] - zc);
    double r = delx * delx + dely * dely + delz * delz;
    double rc = a * a * b * b * c * c;
    if (r > rc || r == 0.0) return 0;

    double a_r = a - cutoff;
    double b_r = b - cutoff;
    double c_r = c - cutoff;
    double delx_r = b_r * c_r * (x[0] - xc);
    double dely_r = a_r * c_r * (x[1] - yc);
    double delz_r = a_r * b_r * (x[2] - zc);
    double r_r = delx_r * delx_r + dely_r * dely_r + delz_r * delz_r;
    double rc_r = a_r * a_r * b_r * b_r * c_r * c_r;

    if (r_r > rc_r) {
      // sort the values
      double axes[3] = {a, b, c};
      double coords[3] = {fabs(x[0] - xc), fabs(x[1] - yc), fabs(x[2] - zc)};

      int min, max;
      if (axes[0] < axes[1]) {
        min = 0;
        max = 1;
      } else {
        min = 1;
        max = 0;
      }
      if (axes[min] > axes[2]) { min = 2; }
      if (axes[max] < axes[2]) { max = 2; }
      int mid = 3 - min - max;
      int sorting[3] = {min, mid, max};

      double x0[3];
      contact[0].r = DistancePointEllipsoid(
          axes[sorting[2]], axes[sorting[1]], axes[sorting[0]], coords[sorting[2]],
          coords[sorting[1]], coords[sorting[0]], x0[sorting[2]], x0[sorting[1]], x0[sorting[0]]);
      contact[0].delx = x[0] - (copysign(x0[0], x[0] - xc) + xc);
      contact[0].dely = x[1] - (copysign(x0[1], x[1] - yc) + yc);
      contact[0].delz = x[2] - (copysign(x0[2], x[2] - zc) + zc);
      //      contact[0].radius = -radius;
      contact[0].iwall = 0;
      contact[0].varflag = 1;
      return 1;
    }
    return 0;
  } else {
    double delx = b * (x[0] - xc);
    double dely = a * (x[1] - yc);
    double r = delx * delx + dely * dely;
    double rc = a * a * b * b;
    if (r > rc || r == 0.0) return 0;

    double a_r = a - cutoff;
    double b_r = b - cutoff;
    double delx_r = b_r * (x[0] - xc);
    double dely_r = a_r * (x[1] - yc);
    double r_r = delx_r * delx_r + dely_r * dely_r;
    double rc_r = a_r * a_r * b_r * b_r;

    if (r_r > rc_r) {
      double x0, x1;
      if (a >= b) {
        contact[0].r = DistancePointEllipse(a, b, fabs(x[0] - xc), fabs(x[1] - yc), x0, x1);
        contact[0].delx = x[0] - (copysign(x0, x[0] - xc) + xc);
        contact[0].dely = x[1] - (copysign(x1, x[1] - yc) + yc);
      } else {
        contact[0].r = DistancePointEllipse(b, a, fabs(x[1] - yc), fabs(x[0] - xc), x0, x1);
        contact[0].delx = x[0] - (copysign(x1, x[0] - xc) + xc);
        contact[0].dely = x[1] - (copysign(x0, x[1] - yc) + yc);
      }
      contact[0].delz = 0;
      //     contact[0].radius = -radius;
      contact[0].iwall = 0;
      contact[0].varflag = 1;
      return 1;
    }
    return 0;
  }
}

/* ----------------------------------------------------------------------
   one contact if 0 <= x < cutoff from outer surface of ellipsoid
   no contact if inside (possible if called from union/intersect)
   delxyz = vector from nearest point on ellipsoid to x
------------------------------------------------------------------------- */

int RegEllipsoid::surface_exterior(double *x, double cutoff)
{
  if (domain->dimension == 3) {
    double delx = b * c * (x[0] - xc);
    double dely = a * c * (x[1] - yc);
    double delz = a * b * (x[2] - zc);
    double r = delx * delx + dely * dely + delz * delz;
    double rc = a * a * b * b * c * c;
    if (r < rc) return 0;

    double a_r = a + cutoff;
    double b_r = b + cutoff;
    double c_r = c + cutoff;
    double delx_r = b_r * c_r * (x[0] - xc);
    double dely_r = a_r * c_r * (x[1] - yc);
    double delz_r = a_r * b_r * (x[2] - zc);
    double r_r = delx_r * delx_r + dely_r * dely_r + delz_r * delz_r;
    double rc_r = a_r * a_r * b_r * b_r * c_r * c_r;

    if (r_r < rc_r) {
      // sort the values
      double axes[3] = {a, b, c};
      double coords[3] = {fabs(x[0] - xc), fabs(x[1] - yc), fabs(x[2] - zc)};

      int min, max;
      if (axes[0] < axes[1]) {
        min = 0;
        max = 1;
      } else {
        min = 1;
        max = 0;
      }
      if (axes[min] > axes[2]) { min = 2; }
      if (axes[max] < axes[2]) { max = 2; }
      int mid = 3 - min - max;
      int sorting[3] = {min, mid, max};

      double x0[3];
      contact[0].r = DistancePointEllipsoid(
          axes[sorting[2]], axes[sorting[1]], axes[sorting[0]], coords[sorting[2]],
          coords[sorting[1]], coords[sorting[0]], x0[sorting[2]], x0[sorting[1]], x0[sorting[0]]);
      contact[0].delx = x[0] - (copysign(x0[0], x[0] - xc) + xc);
      contact[0].dely = x[1] - (copysign(x0[1], x[1] - yc) + yc);
      contact[0].delz = x[2] - (copysign(x0[2], x[2] - zc) + zc);
      //      contact[0].radius = radius;
      contact[0].iwall = 0;
      contact[0].varflag = 1;
      return 1;
    }
    return 0;
  } else {
    double delx = b * c * (x[0] - xc);
    double dely = a * c * (x[1] - yc);
    double r = delx * delx + dely * dely;
    double rc = a * a * b * b;
    if (r < rc) return 0;

    double a_r = a + cutoff;
    double b_r = b + cutoff;
    double delx_r = b_r * (x[0] - xc);
    double dely_r = a_r * (x[1] - yc);
    double r_r = delx_r * delx_r + dely_r * dely_r;
    double rc_r = a_r * a_r * b_r * b_r;

    if (r_r < rc_r) {
      double x0, x1;
      if (a >= b) {
        contact[0].r = DistancePointEllipse(a, b, fabs(x[0] - xc), fabs(x[1] - yc), x0, x1);
        contact[0].delx = x[0] - (copysign(x0, x[0] - xc) + xc);
        contact[0].dely = x[1] - (copysign(x1, x[1] - yc) + yc);
      } else {
        contact[0].r = DistancePointEllipse(b, a, fabs(x[1] - yc), fabs(x[0] - xc), x0, x1);
        contact[0].delx = x[0] - (copysign(x1, x[0] - xc) + xc);
        contact[0].dely = x[1] - (copysign(x0, x[1] - yc) + yc);
      }
      contact[0].delz = 0;
      //    contact[0].radius = radius;
      contact[0].iwall = 0;
      contact[0].varflag = 1;
      return 1;
    }
    return 0;
  }
}

/* ----------------------------------------------------------------------
   change region shape via variable evaluation
------------------------------------------------------------------------- */

void RegEllipsoid::shape_update()
{
  if (xstyle == VARIABLE) xc = xscale * input->variable->compute_equal(xvar);
  if (ystyle == VARIABLE) yc = yscale * input->variable->compute_equal(yvar);
  if (zstyle == VARIABLE) zc = zscale * input->variable->compute_equal(zvar);

  if (astyle == VARIABLE) {
    a = xscale * input->variable->compute_equal(avar);
    if (a < 0.0) error->one(FLERR, "Variable evaluation in region gave bad value");
  }
  if (bstyle == VARIABLE) {
    b = yscale * input->variable->compute_equal(bvar);
    if (b < 0.0) error->one(FLERR, "Variable evaluation in region gave bad value");
  }
  if (cstyle == VARIABLE) {
    c = zscale * input->variable->compute_equal(cvar);
    if (c < 0.0) error->one(FLERR, "Variable evaluation in region gave bad value");
  }
}

/* ----------------------------------------------------------------------
   error check on existence of variable
------------------------------------------------------------------------- */

void RegEllipsoid::variable_check()
{
  if (xstyle == VARIABLE) {
    xvar = input->variable->find(xstr);
    if (xvar < 0) error->all(FLERR, "Variable name {} for region ellipsoid does not exist", xstr);
    if (!input->variable->equalstyle(xvar))
      error->all(FLERR, "Variable {} for region ellipsoid is invalid style", xstr);
  }

  if (ystyle == VARIABLE) {
    yvar = input->variable->find(ystr);
    if (yvar < 0) error->all(FLERR, "Variable name {} for region ellipsoid does not exist", ystr);
    if (!input->variable->equalstyle(yvar))
      error->all(FLERR, "Variable {} for region ellipsoid is invalid style", ystr);
  }

  if (zstyle == VARIABLE) {
    zvar = input->variable->find(zstr);
    if (zvar < 0) error->all(FLERR, "Variable name {} for region ellipsoid does not exist", zstr);
    if (!input->variable->equalstyle(zvar))
      error->all(FLERR, "Variable {} for region ellipsoid is invalid style", zstr);
  }

  if (astyle == VARIABLE) {
    avar = input->variable->find(astr);
    if (avar < 0) error->all(FLERR, "Variable name {} for region ellipsoid does not exist", astr);
    if (!input->variable->equalstyle(avar))
      error->all(FLERR, "Variable {} for region ellipsoid is invalid style", astr);
  }

  if (bstyle == VARIABLE) {
    bvar = input->variable->find(bstr);
    if (bvar < 0) error->all(FLERR, "Variable name {} for region ellipsoid does not exist", bstr);
    if (!input->variable->equalstyle(bvar))
      error->all(FLERR, "Variable {} for region ellipsoid is invalid style", bstr);
  }

  if (cstyle == VARIABLE) {
    cvar = input->variable->find(cstr);
    if (cvar < 0) error->all(FLERR, "Variable name {} for region ellipsoid does not exist", cstr);
    if (!input->variable->equalstyle(cvar))
      error->all(FLERR, "Variable {} for region ellipsoid is invalid style", cstr);
  }
}

// ------------------------------------------------------------------
// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2021
// Based on https://www.geometrictools.com/Documentation/DistancePointEllipseEllipsoid.pdf
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
// https://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// Version: 4.0.2021.08.01
// ------------------------------------------------------------------

/* ----------------------------------------------------------------------
   static helper functions for the 2D case
------------------------------------------------------------------------- */

double GetRoot2D(double r0, double z0, double z1, double g)
{
  const double n0 = r0 * z0;
  double s0 = z1 - 1.0;
  double s1 = (g < 0.0 ? 0.0 : sqrt(n0 * n0 + z1 * z1) - 1.0);
  double s = 0.0;
  for (int i = 0; i < maxIterations; ++i) {
    s = (s0 + s1) / 2.0;
    if (s == s0 || s == s1) { break; }
    const double ratio0 = n0 / (s + r0);
    const double ratio1 = z1 / (s + 1.0);
    g = ratio0 * ratio0 + ratio1 * ratio1 - 1.0;
    if ((g > 0.0) && (g > EPSILON)) {
      s0 = s;
    } else if ((g < 0.0) && (g < -EPSILON)) {
      s1 = s;
    } else {
      break;
    }
  }
  return s;
}

double DistancePointEllipse(double e0, double e1, double y0, double y1, double &x0, double &x1)
{
  double distance;
  if (y1 > 0.0) {
    if (y0 > 0.0) {
      double z0 = y0 / e0;
      double z1 = y1 / e1;
      double g = z0 * z0 + z1 * z1 - 1.0;
      if (g != 0.0) {
        double r0 = (e0 * e0) / (e1 * e1);
        double sbar = GetRoot2D(r0, z0, z1, g);
        x0 = r0 * y0 / (sbar + r0);
        x1 = y1 / (sbar + 1.0);
        distance = sqrt((x0 - y0) * (x0 - y0) + (x1 - y1) * (x1 - y1));
      } else {
        x0 = y0;
        x1 = y1;
        distance = 0.0;
      }
    } else {
      x0 = 0.0;
      x1 = e1;
      distance = fabs(y1 - e1);
    }
  } else {
    double numer0 = e0 * y0;
    double denom0 = e0 * e0 - e1 * e1;
    if (numer0 < denom0) {
      double xde0 = numer0 / denom0;
      x0 = e0 * xde0;
      x1 = e1 * sqrt(1 - xde0 * xde0);
      distance = sqrt((x0 - y0) * (x0 - y0) + x1 * x1);
    } else {
      x0 = e0;
      x1 = 0.0;
      distance = fabs(y0 - e0);
    }
  }
  return distance;
}

/* ----------------------------------------------------------------------
   static helper functions for the 3D case
------------------------------------------------------------------------- */

double GetRoot3D(double r0, double r1, double z0, double z1, double z2, double g)
{
  const double n0 = r0 * z0;
  const double n1 = r1 * z1;
  double s0 = z2 - 1.0;
  double s1 = (g < 0.0 ? 0.0 : sqrt(n0 * n0 + n1 * n1 + z2 * z2) - 1.0);
  double s = 0.0;
  for (int i = 0; i < maxIterations; ++i) {
    s = (s0 + s1) / 2.0;
    if (s == s0 || s == s1) { break; }
    const double ratio0 = n0 / (s + r0);
    const double ratio1 = n1 / (s + r1);
    const double ratio2 = z2 / (s + 1.0);
    g = ratio0 * ratio0 + ratio1 * ratio1 + ratio2 * ratio2 - 1.0;
    if ((g > 0.0) && (g > EPSILON)) {
      s0 = s;
    } else if ((g < 0.0) && (g < -EPSILON)) {
      s1 = s;
    } else {
      break;
    }
  }
  return s;
}

double DistancePointEllipsoid(double e0, double e1, double e2, double y0, double y1, double y2,
                              double &x0, double &x1, double &x2)
{
  double distance;
  if (y2 > 0.0) {
    if (y1 > 0.0) {
      if (y0 > 0.0) {
        double z0 = y0 / e0;
        double z1 = y1 / e1;
        double z2 = y2 / e2;
        double g = z0 * z0 + z1 * z1 + z2 * z2 - 1.0;
        if (g != 0.0) {
          double r0 = e0 * e0 / (e2 * e2);
          double r1 = e1 * e1 / (e2 * e2);
          double sbar = GetRoot3D(r0, r1, z0, z1, z2, g);
          x0 = r0 * y0 / (sbar + r0);
          x1 = r1 * y1 / (sbar + r1);
          x2 = y2 / (sbar + 1.0);
          distance = sqrt((x0 - y0) * (x0 - y0) + (x1 - y1) * (x1 - y1) + (x2 - y2) * (x2 - y2));
        } else {
          x0 = y0;
          x1 = y1;
          x2 = y2;
          distance = 0.0;
        }
      } else {
        x0 = 0.0;
        distance = DistancePointEllipse(e1, e2, y1, y2, x1, x2);
      }
    } else {
      if (y0 > 0.0) {
        x1 = 0.0;
        distance = DistancePointEllipse(e0, e2, y0, y2, x0, x2);
      } else {
        x0 = 0.0;
        x1 = 0.0;
        x2 = e2;
        distance = fabs(y2 - e2);
      }
    }
  } else {
    double denom0 = e0 * e0 - e2 * e2;
    double denom1 = e1 * e1 - e2 * e2;
    double numer0 = e0 * y0;
    double numer1 = e1 * y1;
    bool computed = false;
    if (numer0 < denom0 && numer1 < denom1) {
      double xde0 = numer0 / denom0;
      double xde1 = numer1 / denom1;
      double xde0sqr = xde0 * xde0;
      double xde1sqr = xde1 * xde1;
      double discr = 1.0 - xde0sqr - xde1sqr;
      if (discr > 0.0) {
        x0 = e0 * xde0;
        x1 = e1 * xde1;
        x2 = e2 * sqrt(discr);
        distance = sqrt((x0 - y0) * (x0 - y0) + (x1 - y1) * (x1 - y1) + x2 * x2);
        computed = true;
      }
    }
    if (!computed) {
      x2 = 0.0;
      distance = DistancePointEllipse(e0, e1, y0, y1, x0, x1);
    }
  }
  return distance;
}
