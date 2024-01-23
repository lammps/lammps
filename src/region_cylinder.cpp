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

#include "region_cylinder.h"

#include "domain.h"
#include "error.h"
#include "input.h"
#include "update.h"
#include "variable.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

static constexpr double BIG = 1.0e20;

enum { CONSTANT, VARIABLE };

/* ---------------------------------------------------------------------- */

RegCylinder::RegCylinder(LAMMPS *lmp, int narg, char **arg) :
    Region(lmp, narg, arg), c1str(nullptr), c2str(nullptr), rstr(nullptr)
{
  options(narg - 8, &arg[8]);

  // check open face settings

  if (openflag)
    for (int i = 3; i < 6; i++)
      if (open_faces[i]) error->all(FLERR, "Illegal region cylinder open face: {}", i + 1);

  if (strcmp(arg[2], "x") != 0 && strcmp(arg[2], "y") != 0 && strcmp(arg[2], "z") != 0)
    error->all(FLERR, "Illegal region cylinder axis: {}", arg[2]);
  axis = arg[2][0];

  if (axis == 'x') {
    if (utils::strmatch(arg[3], "^v_")) {
      c1str = utils::strdup(arg[3] + 2);
      c1 = 0.0;
      c1style = VARIABLE;
      varshape = 1;
    } else {
      c1 = yscale * utils::numeric(FLERR, arg[3], false, lmp);
      c1style = CONSTANT;
    }
    if (utils::strmatch(arg[4], "^v_")) {
      c2str = utils::strdup(arg[4] + 2);
      c2 = 0.0;
      c2style = VARIABLE;
      varshape = 1;
    } else {
      c2 = zscale * utils::numeric(FLERR, arg[4], false, lmp);
      c2style = CONSTANT;
    }
  } else if (axis == 'y') {
    if (utils::strmatch(arg[3], "^v_")) {
      c1str = utils::strdup(arg[3] + 2);
      c1 = 0.0;
      c1style = VARIABLE;
      varshape = 1;
    } else {
      c1 = xscale * utils::numeric(FLERR, arg[3], false, lmp);
      c1style = CONSTANT;
    }
    if (utils::strmatch(arg[4], "^v_")) {
      c2str = utils::strdup(arg[4] + 2);
      c2 = 0.0;
      c2style = VARIABLE;
      varshape = 1;
    } else {
      c2 = zscale * utils::numeric(FLERR, arg[4], false, lmp);
      c2style = CONSTANT;
    }
  } else if (axis == 'z') {
    if (utils::strmatch(arg[3], "^v_")) {
      c1str = utils::strdup(arg[3] + 2);
      c1 = 0.0;
      c1style = VARIABLE;
      varshape = 1;
    } else {
      c1 = xscale * utils::numeric(FLERR, arg[3], false, lmp);
      c1style = CONSTANT;
    }
    if (utils::strmatch(arg[4], "^v_")) {
      c2str = utils::strdup(arg[4] + 2);
      c2 = 0.0;
      c2style = VARIABLE;
      varshape = 1;
    } else {
      c2 = yscale * utils::numeric(FLERR, arg[4], false, lmp);
      c2style = CONSTANT;
    }
  }

  if (utils::strmatch(arg[5], "^v_")) {
    rstr = utils::strdup(arg[5] + 2);
    radius = 0.0;
    rstyle = VARIABLE;
    varshape = 1;
  } else {
    radius = utils::numeric(FLERR, arg[5], false, lmp);
    if (axis == 'x')
      radius *= yscale;
    else
      radius *= xscale;
    rstyle = CONSTANT;
  }

  if (varshape) {
    variable_check();
    RegCylinder::shape_update();
  }

  if (strcmp(arg[6], "INF") == 0 || strcmp(arg[6], "EDGE") == 0) {
    if (domain->box_exist == 0)
      error->all(FLERR, "Cannot use region INF or EDGE when box does not exist");
    if (axis == 'x') {
      if (strcmp(arg[6], "INF") == 0)
        lo = -BIG;
      else if (domain->triclinic == 0)
        lo = domain->boxlo[0];
      else
        lo = domain->boxlo_bound[0];
    }
    if (axis == 'y') {
      if (strcmp(arg[6], "INF") == 0)
        lo = -BIG;
      else if (domain->triclinic == 0)
        lo = domain->boxlo[1];
      else
        lo = domain->boxlo_bound[1];
    }
    if (axis == 'z') {
      if (strcmp(arg[6], "INF") == 0)
        lo = -BIG;
      else if (domain->triclinic == 0)
        lo = domain->boxlo[2];
      else
        lo = domain->boxlo_bound[2];
    }
  } else {
    if (axis == 'x') lo = xscale * utils::numeric(FLERR, arg[6], false, lmp);
    if (axis == 'y') lo = yscale * utils::numeric(FLERR, arg[6], false, lmp);
    if (axis == 'z') lo = zscale * utils::numeric(FLERR, arg[6], false, lmp);
  }

  if (strcmp(arg[7], "INF") == 0 || strcmp(arg[7], "EDGE") == 0) {
    if (domain->box_exist == 0)
      error->all(FLERR, "Cannot use region INF or EDGE when box does not exist");
    if (axis == 'x') {
      if (strcmp(arg[7], "INF") == 0)
        hi = BIG;
      else if (domain->triclinic == 0)
        hi = domain->boxhi[0];
      else
        hi = domain->boxhi_bound[0];
    }
    if (axis == 'y') {
      if (strcmp(arg[7], "INF") == 0)
        hi = BIG;
      else if (domain->triclinic == 0)
        hi = domain->boxhi[1];
      else
        hi = domain->boxhi_bound[1];
    }
    if (axis == 'z') {
      if (strcmp(arg[7], "INF") == 0)
        hi = BIG;
      else if (domain->triclinic == 0)
        hi = domain->boxhi[2];
      else
        hi = domain->boxhi_bound[2];
    }
  } else {
    if (axis == 'x') hi = xscale * utils::numeric(FLERR, arg[7], false, lmp);
    if (axis == 'y') hi = yscale * utils::numeric(FLERR, arg[7], false, lmp);
    if (axis == 'z') hi = zscale * utils::numeric(FLERR, arg[7], false, lmp);
  }

  // error check

  if (radius <= 0.0) error->all(FLERR, "Illegal radius {} in region cylinder command", radius);

  // extent of cylinder
  // for variable radius, uses initial radius

  if (interior) {
    bboxflag = 1;
    if (axis == 'x') {
      extent_xlo = lo;
      extent_xhi = hi;
      extent_ylo = c1 - radius;
      extent_yhi = c1 + radius;
      extent_zlo = c2 - radius;
      extent_zhi = c2 + radius;
    }
    if (axis == 'y') {
      extent_xlo = c1 - radius;
      extent_xhi = c1 + radius;
      extent_ylo = lo;
      extent_yhi = hi;
      extent_zlo = c2 - radius;
      extent_zhi = c2 + radius;
    }
    if (axis == 'z') {
      extent_xlo = c1 - radius;
      extent_xhi = c1 + radius;
      extent_ylo = c2 - radius;
      extent_yhi = c2 + radius;
      extent_zlo = lo;
      extent_zhi = hi;
    }
  } else
    bboxflag = 0;

  // particle could be close to cylinder surface and 2 ends
  // particle can only touch surface and 1 end

  cmax = 3;
  contact = new Contact[cmax];
  if (interior)
    tmax = 2;
  else
    tmax = 1;
}

/* ---------------------------------------------------------------------- */

RegCylinder::~RegCylinder()
{
  delete[] c1str;
  delete[] c2str;
  delete[] rstr;
  delete[] contact;
}

/* ---------------------------------------------------------------------- */

void RegCylinder::init()
{
  Region::init();
  if (varshape) variable_check();
}

/* ----------------------------------------------------------------------
   inside = 1 if x,y,z is inside or on surface
   inside = 0 if x,y,z is outside and not on surface
------------------------------------------------------------------------- */

int RegCylinder::inside(double x, double y, double z)
{
  double del1, del2, dist;
  int inside;

  if (axis == 'x') {
    del1 = y - c1;
    del2 = z - c2;
    dist = sqrt(del1 * del1 + del2 * del2);
    if (dist <= radius && x >= lo && x <= hi)
      inside = 1;
    else
      inside = 0;
  } else if (axis == 'y') {
    del1 = x - c1;
    del2 = z - c2;
    dist = sqrt(del1 * del1 + del2 * del2);
    if (dist <= radius && y >= lo && y <= hi)
      inside = 1;
    else
      inside = 0;
  } else {
    del1 = x - c1;
    del2 = y - c2;
    dist = sqrt(del1 * del1 + del2 * del2);
    if (dist <= radius && z >= lo && z <= hi)
      inside = 1;
    else
      inside = 0;
  }

  return inside;
}

/* ----------------------------------------------------------------------
   contact if 0 <= x < cutoff from one or more inner surfaces of cylinder
   can be one contact for each of 3 cylinder surfaces
   no contact if outside (possible if called from union/intersect)
   delxyz = vector from nearest point on cylinder to x
   special case: no contact with curved surf if x is on center axis
------------------------------------------------------------------------- */

int RegCylinder::surface_interior(double *x, double cutoff)
{
  double del1, del2, r, delta;

  int n = 0;

  if (axis == 'x') {
    del1 = x[1] - c1;
    del2 = x[2] - c2;
    r = sqrt(del1 * del1 + del2 * del2);

    // x is exterior to cylinder

    if (r > radius || x[0] < lo || x[0] > hi) return 0;

    // x is interior to cylinder or on its surface

    delta = radius - r;
    if (delta < cutoff && r > 0.0 && !open_faces[2]) {
      contact[n].r = delta;
      contact[n].delx = 0.0;
      contact[n].dely = del1 * (1.0 - radius / r);
      contact[n].delz = del2 * (1.0 - radius / r);
      contact[n].radius = -2.0 * radius;
      contact[n].iwall = 2;
      contact[n].varflag = 1;
      n++;
    }
    delta = x[0] - lo;
    if (delta < cutoff && !open_faces[0]) {
      contact[n].r = delta;
      contact[n].delx = delta;
      contact[n].dely = contact[n].delz = 0.0;
      contact[n].radius = 0;
      contact[n].iwall = 0;
      contact[n].varflag = 0;
      n++;
    }
    delta = hi - x[0];
    if (delta < cutoff && !open_faces[1]) {
      contact[n].r = delta;
      contact[n].delx = -delta;
      contact[n].dely = contact[n].delz = 0.0;
      contact[n].radius = 0;
      contact[n].iwall = 1;
      contact[n].varflag = 0;
      n++;
    }

  } else if (axis == 'y') {
    del1 = x[0] - c1;
    del2 = x[2] - c2;
    r = sqrt(del1 * del1 + del2 * del2);

    // y is exterior to cylinder

    if (r > radius || x[1] < lo || x[1] > hi) return 0;

    // y is interior to cylinder or on its surface

    delta = radius - r;
    if (delta < cutoff && r > 0.0 && !open_faces[2]) {
      contact[n].r = delta;
      contact[n].delx = del1 * (1.0 - radius / r);
      contact[n].dely = 0.0;
      contact[n].delz = del2 * (1.0 - radius / r);
      contact[n].radius = -2.0 * radius;
      contact[n].iwall = 2;
      contact[n].varflag = 1;
      n++;
    }
    delta = x[1] - lo;
    if (delta < cutoff && !open_faces[0]) {
      contact[n].r = delta;
      contact[n].dely = delta;
      contact[n].delx = contact[n].delz = 0.0;
      contact[n].radius = 0;
      contact[n].iwall = 0;
      contact[n].varflag = 0;
      n++;
    }
    delta = hi - x[1];
    if (delta < cutoff && !open_faces[1]) {
      contact[n].r = delta;
      contact[n].dely = -delta;
      contact[n].delx = contact[n].delz = 0.0;
      contact[n].radius = 0;
      contact[n].iwall = 1;
      contact[n].varflag = 0;
      n++;
    }

  } else {
    del1 = x[0] - c1;
    del2 = x[1] - c2;
    r = sqrt(del1 * del1 + del2 * del2);

    // z is exterior to cylinder

    if (r > radius || x[2] < lo || x[2] > hi) return 0;

    // z is interior to cylinder or on its surface

    delta = radius - r;
    if (delta < cutoff && r > 0.0 && !open_faces[2]) {
      contact[n].r = delta;
      contact[n].delx = del1 * (1.0 - radius / r);
      contact[n].dely = del2 * (1.0 - radius / r);
      contact[n].delz = 0.0;
      contact[n].radius = -2.0 * radius;
      contact[n].iwall = 2;
      contact[n].varflag = 1;
      n++;
    }
    delta = x[2] - lo;
    if (delta < cutoff && !open_faces[0]) {
      contact[n].r = delta;
      contact[n].delz = delta;
      contact[n].delx = contact[n].dely = 0.0;
      contact[n].radius = 0;
      contact[n].iwall = 0;
      contact[n].varflag = 0;
      n++;
    }
    delta = hi - x[2];
    if (delta < cutoff && !open_faces[1]) {
      contact[n].r = delta;
      contact[n].delz = -delta;
      contact[n].delx = contact[n].dely = 0.0;
      contact[n].radius = 0;
      contact[n].iwall = 1;
      contact[n].varflag = 0;
      n++;
    }
  }

  return n;
}

/* ----------------------------------------------------------------------
   one contact if 0 <= x < cutoff from outer surface of cylinder
   no contact if inside (possible if called from union/intersect)
   delxyz = vector from nearest point on cylinder to x
------------------------------------------------------------------------- */

int RegCylinder::surface_exterior(double *x, double cutoff)
{
  double del1, del2, r;
  double xp, yp, zp;
  double dx, dr, dr2, d2, d2prev;

  // radius of curvature for granular
  // 0 for flat surfaces (infinite case), 2*radius for curved portion

  double crad = 0.0;
  int varflag = 0;

  if (axis == 'x') {
    del1 = x[1] - c1;
    del2 = x[2] - c2;
    r = sqrt(del1 * del1 + del2 * del2);

    // x is far enough from cylinder that there is no contact
    // x is interior to cylinder

    if (r >= radius + cutoff || x[0] <= lo - cutoff || x[0] >= hi + cutoff) return 0;
    if (r < radius && x[0] > lo && x[0] < hi) return 0;

    // x is exterior to cylinder or on its surface
    // xp,yp,zp = point on surface of cylinder that x is closest to
    //            could be edge of cylinder
    // do not add contact point if r >= cutoff

    d2prev = BIG;
    if (!openflag) {
      if (r > radius) {
        yp = c1 + del1 * radius / r;
        zp = c2 + del2 * radius / r;
        crad = 2.0 * radius;
        varflag = 1;
      } else {
        yp = x[1];
        zp = x[2];
      }
      if (x[0] < lo)
        xp = lo;
      else if (x[0] > hi)
        xp = hi;
      else
        xp = x[0];

    } else {

      // closest point on curved surface

      dr = r - radius;
      dr2 = dr * dr;
      if (!open_faces[2]) {
        yp = c1 + del1 * radius / r;
        zp = c2 + del2 * radius / r;
        if (x[0] < lo) {
          dx = lo - x[0];
          xp = lo;
        } else if (x[0] > hi) {
          dx = x[0] - hi;
          xp = hi;
        } else {
          dx = 0;
          xp = x[0];
        }
        d2 = d2prev = dr2 + dx * dx;
      }

      // closest point on bottom cap

      if (!open_faces[0]) {
        dx = lo - x[0];
        if (r < radius)
          d2 = dx * dx;
        else
          d2 = dr2 + dx * dx;
        if (d2 < d2prev) {
          xp = lo;
          if (r < radius) {
            yp = x[1];
            zp = x[2];
          }
          d2prev = d2;
        }
      }

      // closest point on top cap

      if (!open_faces[1]) {
        dx = hi - x[0];
        if (r < radius)
          d2 = dx * dx;
        else
          d2 = dr2 + dx * dx;
        if (d2 < d2prev) {
          xp = hi;
          if (r < radius) {
            yp = x[1];
            zp = x[2];
          }
        }
      }
    }

    add_contact(0, x, xp, yp, zp);
    contact[0].radius = crad;
    contact[0].varflag = varflag;
    contact[0].iwall = 0;
    if (contact[0].r < cutoff) return 1;
    return 0;

  } else if (axis == 'y') {
    del1 = x[0] - c1;
    del2 = x[2] - c2;
    r = sqrt(del1 * del1 + del2 * del2);

    // y is far enough from cylinder that there is no contact
    // y is interior to cylinder

    if (r >= radius + cutoff || x[1] <= lo - cutoff || x[1] >= hi + cutoff) return 0;
    if (r < radius && x[1] > lo && x[1] < hi) return 0;

    // y is exterior to cylinder or on its surface
    // xp,yp,zp = point on surface of cylinder that x is closest to
    //            could be edge of cylinder
    // do not add contact point if r >= cutoff

    d2prev = BIG;
    if (!openflag) {
      if (r > radius) {
        xp = c1 + del1 * radius / r;
        zp = c2 + del2 * radius / r;
        crad = 2.0 * radius;
        varflag = 1;
      } else {
        xp = x[0];
        zp = x[2];
      }
      if (x[1] < lo)
        yp = lo;
      else if (x[1] > hi)
        yp = hi;
      else
        yp = x[1];

    } else {

      // closest point on curved surface

      dr = r - radius;
      dr2 = dr * dr;
      if (!open_faces[2]) {
        xp = c1 + del1 * radius / r;
        zp = c2 + del2 * radius / r;
        if (x[1] < lo) {
          dx = lo - x[1];
          yp = lo;
        } else if (x[1] > hi) {
          dx = x[1] - hi;
          yp = hi;
        } else {
          dx = 0;
          yp = x[1];
        }
        d2 = d2prev = dr2 + dx * dx;
      }

      // closest point on bottom cap

      if (!open_faces[0]) {
        dx = lo - x[1];
        if (r < radius)
          d2 = dx * dx;
        else
          d2 = dr2 + dx * dx;
        if (d2 < d2prev) {
          yp = lo;
          if (r < radius) {
            xp = x[0];
            zp = x[2];
          }
          d2prev = d2;
        }
      }

      // closest point on top cap

      if (!open_faces[1]) {
        dx = hi - x[1];
        if (r < radius)
          d2 = dx * dx;
        else
          d2 = dr2 + dx * dx;
        if (d2 < d2prev) {
          yp = hi;
          if (r < radius) {
            xp = x[0];
            zp = x[2];
          }
        }
      }
    }

    add_contact(0, x, xp, yp, zp);
    contact[0].radius = crad;
    contact[0].varflag = varflag;
    contact[0].iwall = 0;
    if (contact[0].r < cutoff) return 1;
    return 0;

  } else {
    del1 = x[0] - c1;
    del2 = x[1] - c2;
    r = sqrt(del1 * del1 + del2 * del2);

    // z is far enough from cylinder that there is no contact
    // z is interior to cylinder

    if (r >= radius + cutoff || x[2] <= lo - cutoff || x[2] >= hi + cutoff) return 0;
    if (r < radius && x[2] > lo && x[2] < hi) return 0;

    // z is exterior to cylinder or on its surface
    // xp,yp,zp = point on surface of cylinder that x is closest to
    //            could be edge of cylinder
    // do not add contact point if r >= cutoff

    d2prev = BIG;
    if (!openflag) {
      if (r > radius) {
        xp = c1 + del1 * radius / r;
        yp = c2 + del2 * radius / r;
        crad = 2.0 * radius;
        varflag = 1;
      } else {
        xp = x[0];
        yp = x[1];
      }
      if (x[2] < lo)
        zp = lo;
      else if (x[2] > hi)
        zp = hi;
      else
        zp = x[2];

    } else {

      // closest point on curved surface

      dr = r - radius;
      dr2 = dr * dr;
      if (!open_faces[2]) {
        xp = c1 + del1 * radius / r;
        yp = c2 + del2 * radius / r;
        if (x[2] < lo) {
          dx = lo - x[2];
          zp = lo;
        } else if (x[2] > hi) {
          dx = x[2] - hi;
          zp = hi;
        } else {
          dx = 0;
          zp = x[2];
        }
        d2prev = dr2 + dx * dx;
      }

      // closest point on bottom cap

      if (!open_faces[0]) {
        dx = lo - x[2];
        if (r < radius)
          d2 = dx * dx;
        else
          d2 = dr2 + dx * dx;
        if (d2 < d2prev) {
          zp = lo;
          if (r < radius) {
            xp = x[0];
            yp = x[1];
          }
          d2prev = d2;
        }
      }

      // closest point on top cap

      if (!open_faces[1]) {
        dx = hi - x[2];
        if (r < radius)
          d2 = dx * dx;
        else
          d2 = dr2 + dx * dx;
        if (d2 < d2prev) {
          zp = hi;
          if (r < radius) {
            xp = x[0];
            yp = x[1];
          }
        }
      }
    }

    add_contact(0, x, xp, yp, zp);
    contact[0].radius = crad;
    contact[0].varflag = varflag;
    contact[0].iwall = 0;
    if (contact[0].r < cutoff) return 1;
    return 0;
  }
}

/* ----------------------------------------------------------------------
   change region shape via variable evaluation
------------------------------------------------------------------------- */

void RegCylinder::shape_update()
{
  if (c1style == VARIABLE) c1 = input->variable->compute_equal(c1var);
  if (c2style == VARIABLE) c2 = input->variable->compute_equal(c2var);
  if (rstyle == VARIABLE) {
    radius = input->variable->compute_equal(rvar);
    if (radius < 0.0) error->one(FLERR, "Variable evaluation in region gave bad value");
  }

  if (axis == 'x') {
    if (c1style == VARIABLE) c1 *= yscale;
    if (c2style == VARIABLE) c2 *= zscale;
    if (rstyle == VARIABLE) radius *= yscale;
  } else if (axis == 'y') {
    if (c1style == VARIABLE) c1 *= xscale;
    if (c2style == VARIABLE) c2 *= zscale;
    if (rstyle == VARIABLE) radius *= xscale;
  } else {    // axis == 'z'
    if (c1style == VARIABLE) c1 *= xscale;
    if (c2style == VARIABLE) c2 *= yscale;
    if (rstyle == VARIABLE) radius *= xscale;
  }
}

/* ----------------------------------------------------------------------
   error check on existence of variable
------------------------------------------------------------------------- */

void RegCylinder::variable_check()
{
  if (c1style == VARIABLE) {
    c1var = input->variable->find(c1str);
    if (c1var < 0) error->all(FLERR, "Variable {} for region cylinder does not exist", c1str);
    if (!input->variable->equalstyle(c1var))
      error->all(FLERR, "Variable {} for region cylinder is invalid style", c1str);
  }

  if (c2style == VARIABLE) {
    c2var = input->variable->find(c2str);
    if (c2var < 0) error->all(FLERR, "Variable {} for region cylinder does not exist", c2str);
    if (!input->variable->equalstyle(c2var))
      error->all(FLERR, "Variable {} for region cylinder is invalid style", c2str);
  }

  if (rstyle == VARIABLE) {
    rvar = input->variable->find(rstr);
    if (rvar < 0) error->all(FLERR, "Variable {} for region cylinder does not exist", rstr);
    if (!input->variable->equalstyle(rvar))
      error->all(FLERR, "Variable {} for region cylinder is invalid style", rstr);
  }
}

/* ----------------------------------------------------------------------
   Set values needed to calculate velocity due to shape changes.
   These values do not depend on the contact, so this function is
   called once per timestep by fix/wall/gran/region.

------------------------------------------------------------------------- */

void RegCylinder::set_velocity_shape()
{
  if (axis == 'x') {
    xcenter[0] = 0;
    xcenter[1] = c1;
    xcenter[2] = c2;
  } else if (axis == 'y') {
    xcenter[0] = c1;
    xcenter[1] = 0;
    xcenter[2] = c2;
  } else {
    xcenter[0] = c1;
    xcenter[1] = c2;
    xcenter[2] = 0;
  }
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

void RegCylinder::velocity_contact_shape(double *vwall, double *xc)
{
  double delx, dely, delz;    // Displacement of contact point in x,y,z
  if (axis == 'x') {
    delx = 0;
    dely = (xc[1] - xcenter[1]) * (1 - rprev / radius);
    delz = (xc[2] - xcenter[2]) * (1 - rprev / radius);
  } else if (axis == 'y') {
    delx = (xc[0] - xcenter[0]) * (1 - rprev / radius);
    dely = 0;
    delz = (xc[2] - xcenter[2]) * (1 - rprev / radius);
  } else {
    delx = (xc[0] - xcenter[0]) * (1 - rprev / radius);
    dely = (xc[1] - xcenter[1]) * (1 - rprev / radius);
    delz = 0;
  }
  vwall[0] += delx / update->dt;
  vwall[1] += dely / update->dt;
  vwall[2] += delz / update->dt;
}
