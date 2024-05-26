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

#include "region_sphere.h"

#include "error.h"
#include "input.h"
#include "update.h"
#include "variable.h"

#include <cmath>

using namespace LAMMPS_NS;

enum { CONSTANT, VARIABLE };

/* ---------------------------------------------------------------------- */

RegSphere::RegSphere(LAMMPS *lmp, int narg, char **arg) :
    Region(lmp, narg, arg), xstr(nullptr), ystr(nullptr), zstr(nullptr), rstr(nullptr)
{
  options(narg - 6, &arg[6]);

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
    rstr = utils::strdup(arg[5] + 2);
    radius = 0.0;
    rstyle = VARIABLE;
    varshape = 1;
  } else {
    radius = xscale * utils::numeric(FLERR, arg[5], false, lmp);
    rstyle = CONSTANT;
  }

  if (varshape) {
    variable_check();
    RegSphere::shape_update();
  }

  // error check

  if (radius < 0.0) error->all(FLERR, "Illegal region sphere radius: {}", radius);

  // extent of sphere
  // for variable radius, uses initial radius and origin for variable center

  if (interior) {
    bboxflag = 1;
    extent_xlo = xc - radius;
    extent_xhi = xc + radius;
    extent_ylo = yc - radius;
    extent_yhi = yc + radius;
    extent_zlo = zc - radius;
    extent_zhi = zc + radius;
  } else
    bboxflag = 0;

  cmax = 1;
  contact = new Contact[cmax];
  tmax = 1;
}

/* ---------------------------------------------------------------------- */

RegSphere::~RegSphere()
{
  delete[] xstr;
  delete[] ystr;
  delete[] zstr;
  delete[] rstr;
  delete[] contact;
}

/* ---------------------------------------------------------------------- */

void RegSphere::init()
{
  Region::init();
  if (varshape) variable_check();
}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside or on surface
   inside = 0 if x,y,z is outside and not on surface
------------------------------------------------------------------------- */

int RegSphere::inside(double x, double y, double z)
{
  double delx = x - xc;
  double dely = y - yc;
  double delz = z - zc;
  double r = sqrt(delx * delx + dely * dely + delz * delz);

  if (r <= radius) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   one contact if 0 <= x < cutoff from inner surface of sphere
   no contact if outside (possible if called from union/intersect)
   delxyz = vector from nearest point on sphere to x
   special case: no contact if x is at center of sphere
------------------------------------------------------------------------- */

int RegSphere::surface_interior(double *x, double cutoff)
{
  double delx = x[0] - xc;
  double dely = x[1] - yc;
  double delz = x[2] - zc;
  double r = sqrt(delx * delx + dely * dely + delz * delz);
  if (r > radius || r == 0.0) return 0;

  double delta = radius - r;
  if (delta < cutoff) {
    contact[0].r = delta;
    contact[0].delx = delx * (1.0 - radius / r);
    contact[0].dely = dely * (1.0 - radius / r);
    contact[0].delz = delz * (1.0 - radius / r);
    contact[0].radius = -radius;
    contact[0].iwall = 0;
    contact[0].varflag = 1;
    return 1;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   one contact if 0 <= x < cutoff from outer surface of sphere
   no contact if inside (possible if called from union/intersect)
   delxyz = vector from nearest point on sphere to x
------------------------------------------------------------------------- */

int RegSphere::surface_exterior(double *x, double cutoff)
{
  double delx = x[0] - xc;
  double dely = x[1] - yc;
  double delz = x[2] - zc;
  double r = sqrt(delx * delx + dely * dely + delz * delz);
  if (r < radius) return 0;

  double delta = r - radius;
  if (delta < cutoff) {
    contact[0].r = delta;
    contact[0].delx = delx * (1.0 - radius / r);
    contact[0].dely = dely * (1.0 - radius / r);
    contact[0].delz = delz * (1.0 - radius / r);
    contact[0].radius = radius;
    contact[0].iwall = 0;
    contact[0].varflag = 1;
    return 1;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   change region shape via variable evaluation
------------------------------------------------------------------------- */

void RegSphere::shape_update()
{
  if (xstyle == VARIABLE) xc = xscale * input->variable->compute_equal(xvar);

  if (ystyle == VARIABLE) yc = yscale * input->variable->compute_equal(yvar);

  if (zstyle == VARIABLE) zc = zscale * input->variable->compute_equal(zvar);

  if (rstyle == VARIABLE) {
    radius = xscale * input->variable->compute_equal(rvar);
    if (radius < 0.0) error->one(FLERR, "Variable evaluation in region gave bad value");
  }
}

/* ----------------------------------------------------------------------
   error check on existence of variable
------------------------------------------------------------------------- */

void RegSphere::variable_check()
{
  if (xstyle == VARIABLE) {
    xvar = input->variable->find(xstr);
    if (xvar < 0) error->all(FLERR, "Variable {} for region sphere does not exist", xstr);
    if (!input->variable->equalstyle(xvar))
      error->all(FLERR, "Variable {} for region sphere is invalid style", xstr);
  }

  if (ystyle == VARIABLE) {
    yvar = input->variable->find(ystr);
    if (yvar < 0) error->all(FLERR, "Variable {} for region sphere does not exist", ystr);
    if (!input->variable->equalstyle(yvar))
      error->all(FLERR, "Variable {} for region sphere is invalid style", ystr);
  }

  if (zstyle == VARIABLE) {
    zvar = input->variable->find(zstr);
    if (zvar < 0) error->all(FLERR, "Variable {} for region sphere does not exist", zstr);
    if (!input->variable->equalstyle(zvar))
      error->all(FLERR, "Variable {} for region sphere is invalid style", zstr);
  }

  if (rstyle == VARIABLE) {
    rvar = input->variable->find(rstr);
    if (rvar < 0) error->all(FLERR, "Variable {} for region sphere does not exist", rstr);
    if (!input->variable->equalstyle(rvar))
      error->all(FLERR, "Variable {} for region sphere is invalid style", rstr);
  }
}

/* ----------------------------------------------------------------------
   Set values needed to calculate velocity due to shape changes.
   These values do not depend on the contact, so this function is
   called once per timestep by fix/wall/gran/region.

------------------------------------------------------------------------- */

void RegSphere::set_velocity_shape()
{
  xcenter[0] = xc;
  xcenter[1] = yc;
  xcenter[2] = zc;
  forward_transform(xcenter[0], xcenter[1], xcenter[2]);
  if (update->ntimestep > 0)
    rprev = prev[4];
  else
    rprev = radius;
  prev[4] = radius;
}

/* ----------------------------------------------------------------------
   add velocity due to shape change to wall velocity
------------------------------------------------------------------------- */

void RegSphere::velocity_contact_shape(double *vwall, double *xc)
{
  double delx, dely, delz;    // Displacement of contact point in x,y,z

  delx = (xc[0] - xcenter[0]) * (1 - rprev / radius);
  dely = (xc[1] - xcenter[1]) * (1 - rprev / radius);
  delz = (xc[2] - xcenter[2]) * (1 - rprev / radius);

  vwall[0] += delx / update->dt;
  vwall[1] += dely / update->dt;
  vwall[2] += delz / update->dt;
}
