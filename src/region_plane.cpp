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

#include "region_plane.h"

#include "error.h"
#include "input.h"
#include "update.h"
#include "variable.h"

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

RegPlane::RegPlane(LAMMPS *lmp, int narg, char **arg) : Region(lmp, narg, arg),
    xstr(nullptr), ystr(nullptr), zstr(nullptr)
{
  xvar = yvar = zvar = 0.0;

  options(narg - 8, &arg[8]);

  if (utils::strmatch(arg[2], "^v_")) {
    xstr = utils::strdup(arg[2] + 2);
    xp = 0.0;
    xstyle = VARIABLE;
    varshape = 1;
  } else {
    xp = xscale * utils::numeric(FLERR, arg[2], false, lmp);
    xstyle = CONSTANT;
  }

  if (utils::strmatch(arg[3], "^v_")) {
    ystr = utils::strdup(arg[3] + 2);
    yp = 0.0;
    ystyle = VARIABLE;
    varshape = 1;
  } else {
    yp = yscale * utils::numeric(FLERR, arg[3], false, lmp);
    ystyle = CONSTANT;
  }

  if (utils::strmatch(arg[4], "^v_")) {
    zstr = utils::strdup(arg[4] + 2);
    zp = 0.0;
    zstyle = VARIABLE;
    varshape = 1;
  } else {
    zp = zscale * utils::numeric(FLERR, arg[4], false, lmp);
    zstyle = CONSTANT;
  }

  if (varshape) {
    variable_check();
    RegPlane::shape_update();
  }

  normal[0] = xscale * utils::numeric(FLERR, arg[5], false, lmp);
  normal[1] = yscale * utils::numeric(FLERR, arg[6], false, lmp);
  normal[2] = zscale * utils::numeric(FLERR, arg[7], false, lmp);

  // enforce unit normal

  double rsq = normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2];
  if (rsq == 0.0)
    error->all(FLERR, "Illegal region plane normal vector: {} {} {}", normal[0], normal[1],
               normal[2]);
  normal[0] /= sqrt(rsq);
  normal[1] /= sqrt(rsq);
  normal[2] /= sqrt(rsq);

  // plane has no bounding box

  bboxflag = 0;
  cmax = 1;
  contact = new Contact[cmax];
  tmax = 1;
}

/* ---------------------------------------------------------------------- */

RegPlane::~RegPlane()
{
  delete[] xstr;
  delete[] ystr;
  delete[] zstr;
  delete[] contact;
}

/* ---------------------------------------------------------------------- */

void RegPlane::init()
{
  Region::init();
  if (varshape) variable_check();
}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is on normal side of plane or on plane
   inside = 0 if x,y,z is on non-normal side of plane and not on plane
   x,y,z is inside if (x-xp) dot normal >= 0
------------------------------------------------------------------------- */

int RegPlane::inside(double x, double y, double z)
{
  double dot = (x - xp) * normal[0] + (y - yp) * normal[1] + (z - zp) * normal[2];

  if (dot >= 0.0) return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   one contact if 0 <= x < cutoff from normal side of plane
   no contact if on other side (possible if called from union/intersect)
   delxyz = vector from nearest projected point on plane to x
------------------------------------------------------------------------- */

int RegPlane::surface_interior(double *x, double cutoff)
{
  double dot = (x[0] - xp) * normal[0] + (x[1] - yp) * normal[1] + (x[2] - zp) * normal[2];
  if (dot < cutoff && dot >= 0.0) {
    contact[0].r = dot;
    contact[0].delx = dot * normal[0];
    contact[0].dely = dot * normal[1];
    contact[0].delz = dot * normal[2];
    contact[0].radius = 0;
    contact[0].iwall = 0;
    return 1;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   one contact if 0 <= x < cutoff from non-normal side of plane
   no contact if on other side (possible if called from union/intersect)
   delxyz = vector from nearest projected point on plane to x
------------------------------------------------------------------------- */

int RegPlane::surface_exterior(double *x, double cutoff)
{
  double dot = (x[0] - xp) * normal[0] + (x[1] - yp) * normal[1] + (x[2] - zp) * normal[2];
  dot = -dot;
  if (dot < cutoff && dot >= 0.0) {
    contact[0].r = dot;
    contact[0].delx = -dot * normal[0];
    contact[0].dely = -dot * normal[1];
    contact[0].delz = -dot * normal[2];
    contact[0].radius = 0;
    contact[0].iwall = 0;
    return 1;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   change region shape via variable evaluation
------------------------------------------------------------------------- */

void RegPlane::shape_update()
{
  if (xstyle == VARIABLE) xp = xscale * input->variable->compute_equal(xvar);

  if (ystyle == VARIABLE) yp = yscale * input->variable->compute_equal(yvar);

  if (zstyle == VARIABLE) zp = zscale * input->variable->compute_equal(zvar);
}

/* ----------------------------------------------------------------------
   error check on existence of variable
------------------------------------------------------------------------- */

void RegPlane::variable_check()
{
  if (xstyle == VARIABLE) {
    xvar = input->variable->find(xstr);
    if (xvar < 0) error->all(FLERR, "Variable {} for region plane does not exist", xstr);
    if (!input->variable->equalstyle(xvar))
      error->all(FLERR, "Variable {} for region plane is invalid style", xstr);
  }

  if (ystyle == VARIABLE) {
    yvar = input->variable->find(ystr);
    if (yvar < 0) error->all(FLERR, "Variable {} for region plane does not exist", ystr);
    if (!input->variable->equalstyle(yvar))
      error->all(FLERR, "Variable {} for region plane is invalid style", ystr);
  }

  if (zstyle == VARIABLE) {
    zvar = input->variable->find(zstr);
    if (zvar < 0) error->all(FLERR, "Variable {} for region plane does not exist", zstr);
    if (!input->variable->equalstyle(zvar))
      error->all(FLERR, "Variable {} for region plane is invalid style", zstr);
  }
}

