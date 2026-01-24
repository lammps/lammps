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

#include "fix_graphics_objects.h"

#include "comm.h"
#include "domain.h"
#include "error.h"
#include "graphics.h"
#include "input.h"
#include "lattice.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

#include <cstring>
#include <utility>

using namespace LAMMPS_NS;
using namespace FixConst;

enum { SPHERE, CYLINDER, ARROW, CONE, PROGBAR };
enum { X = 0, Y, Z };

/* ---------------------------------------------------------------------- */

#define PARSE_VARIABLE(value, name, index)      \
  if (strstr(arg[index], "v_") == arg[index]) { \
    varflag = 1;                                \
    delete[] name;                              \
    name = utils::strdup(arg[index] + 2);       \
  } else                                        \
    value = utils::numeric(FLERR, arg[index], false, lmp)

FixGraphicsObjects::FixGraphicsObjects(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), imgobjs(nullptr), imgparms(nullptr)
{
  if (narg < 4) utils::missing_cmd_args(FLERR, "fix graphics/objects", error);

  // parse mandatory arg

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  if (nevery <= 0) error->all(FLERR, 3, "Illegal fix graphics/objects nevery value {}", nevery);
  global_freq = nevery;
  dynamic_group_allow = 1;

  numobjs = 0;
  varflag = 0;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "sphere") == 0) {
      if (iarg + 6 > narg) utils::missing_cmd_args(FLERR, "fix graphics/objects sphere", error);
      // clang-format off
      SphereItem sphere{SPHERE, 1, {0.0, 0.0, 0.0}, 0.0, nullptr, nullptr, nullptr, nullptr,
                        -1, -1, -1, -1};
      // clang-format on
      sphere.type = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      PARSE_VARIABLE(sphere.pos[X], sphere.xstr, iarg + 2);
      PARSE_VARIABLE(sphere.pos[Y], sphere.ystr, iarg + 3);
      PARSE_VARIABLE(sphere.pos[Z], sphere.zstr, iarg + 4);
      PARSE_VARIABLE(sphere.diameter, sphere.dstr, iarg + 5);
      sphere.diameter *= 2.0;
      items.emplace_back(sphere);
      ++numobjs;
      iarg += 6;
    } else if (strcmp(arg[iarg], "cylinder") == 0) {
      if (iarg + 9 > narg) utils::missing_cmd_args(FLERR, "fix graphics/objects cylinder", error);
      // clang-format off
      CylinderItem cylinder{CYLINDER, 1, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0.0,
                            nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                            -1, -1, -1, -1, -1, -1, -1};
      // clang-format on
      cylinder.type = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      PARSE_VARIABLE(cylinder.pos1[X], cylinder.x1str, iarg + 2);
      PARSE_VARIABLE(cylinder.pos1[Y], cylinder.y1str, iarg + 3);
      PARSE_VARIABLE(cylinder.pos1[Z], cylinder.z1str, iarg + 4);
      PARSE_VARIABLE(cylinder.pos2[X], cylinder.x2str, iarg + 5);
      PARSE_VARIABLE(cylinder.pos2[Y], cylinder.y2str, iarg + 6);
      PARSE_VARIABLE(cylinder.pos2[Z], cylinder.z2str, iarg + 7);
      PARSE_VARIABLE(cylinder.diameter, cylinder.dstr, iarg + 8);
      cylinder.diameter *= 2.0;
      items.emplace_back(cylinder);
      ++numobjs;
      iarg += 9;
    } else if (strcmp(arg[iarg], "arrow") == 0) {
      if (iarg + 10 > narg) utils::missing_cmd_args(FLERR, "fix graphics/objects arrow", error);
      // clang-format off
      ArrowItem arrow{ARROW, 1, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0.0, 0.1,
                            nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                            -1, -1, -1, -1, -1, -1, -1};
      // clang-format on
      arrow.type = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      PARSE_VARIABLE(arrow.bot[X], arrow.x1str, iarg + 2);
      PARSE_VARIABLE(arrow.bot[Y], arrow.y1str, iarg + 3);
      PARSE_VARIABLE(arrow.bot[Z], arrow.z1str, iarg + 4);
      PARSE_VARIABLE(arrow.tip[X], arrow.x2str, iarg + 5);
      PARSE_VARIABLE(arrow.tip[Y], arrow.y2str, iarg + 6);
      PARSE_VARIABLE(arrow.tip[Z], arrow.z2str, iarg + 7);
      PARSE_VARIABLE(arrow.diameter, arrow.dstr, iarg + 8);
      arrow.diameter *= 2.0;
      arrow.ratio = utils::numeric(FLERR, arg[iarg + 9], false, lmp);
      if ((arrow.ratio < 0.1) || (arrow.ratio > 0.5))
        error->all(FLERR, iarg + 9, "Arrow tip ratio must be between 0.1 and 0.5");
      items.emplace_back(arrow);
      ++numobjs;
      iarg += 10;
    } else if (strcmp(arg[iarg], "cone") == 0) {
      if (iarg + 11 > narg) utils::missing_cmd_args(FLERR, "fix graphics/objects cone", error);
      // clang-format off
      ConeItem cone{CONE, 1, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0.0, 0.0,
                    nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                    -1, -1, -1, -1, -1, -1, -1, -1, Graphics::CONE_ALL};
      // clang-format on
      cone.type = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      PARSE_VARIABLE(cone.bot[X], cone.x1str, iarg + 2);
      PARSE_VARIABLE(cone.bot[Y], cone.y1str, iarg + 3);
      PARSE_VARIABLE(cone.bot[Z], cone.z1str, iarg + 4);
      PARSE_VARIABLE(cone.top[X], cone.x2str, iarg + 5);
      PARSE_VARIABLE(cone.top[Y], cone.y2str, iarg + 6);
      PARSE_VARIABLE(cone.top[Z], cone.z2str, iarg + 7);
      PARSE_VARIABLE(cone.botdiam, cone.d1str, iarg + 8);
      PARSE_VARIABLE(cone.topdiam, cone.d2str, iarg + 9);
      cone.botdiam *= 2.0;
      cone.topdiam *= 2.0;
      cone.sides = utils::inumeric(FLERR, arg[iarg + 10], false, lmp);
      if ((cone.sides < 0) || (cone.sides > 7))
        error->all(FLERR, iarg + 10, "Cone sides value must be between 0 and 7");
      items.emplace_back(cone);
      ++numobjs;
      iarg += 11;
    } else if (strcmp(arg[iarg], "progbar") == 0) {
      if (iarg + 11 > narg) utils::missing_cmd_args(FLERR, "fix graphics/objects progbar", error);
      // clang-format off
      ProgbarItem progbar{PROGBAR, 1, 2, Y, 0, {0.0, 0.0, 0.0}, 0.0, 0.0, 0.0, nullptr, -1};
      // clang-format on
      progbar.type1 = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      progbar.type2 = utils::inumeric(FLERR, arg[iarg + 2], false, lmp);
      if (strcmp(arg[iarg + 3], "x") == 0) {
        progbar.dim = X;
      } else if (strcmp(arg[iarg + 3], "y") == 0) {
        progbar.dim = Y;
      } else if (strcmp(arg[iarg + 3], "z") == 0) {
        progbar.dim = Z;
      } else {
        error->all(FLERR, iarg + 3, "Unsupported progress bar dimension string {}", arg[iarg + 3]);
      }
      progbar.pos[X] = utils::numeric(FLERR, arg[iarg + 4], false, lmp);
      progbar.pos[Y] = utils::numeric(FLERR, arg[iarg + 5], false, lmp);
      progbar.pos[Z] = utils::numeric(FLERR, arg[iarg + 6], false, lmp);
      progbar.length = utils::numeric(FLERR, arg[iarg + 7], false, lmp);
      if ((progbar.length <= 0.0) || (progbar.length > 2.0 * domain->prd[progbar.dim]))
        error->all(FLERR, iarg + 7, "Illegal progress bar length {}", arg[iarg + 7]);
      progbar.diameter = 2.0 * utils::numeric(FLERR, arg[iarg + 8], false, lmp);
      PARSE_VARIABLE(progbar.progress, progbar.pstr, iarg + 9);
      progbar.tics = utils::inumeric(FLERR, arg[iarg + 10], false, lmp);
      if ((progbar.tics < 0) || (progbar.tics > 20))
        error->all(FLERR, iarg + 10, "Unsupported number of progress bar tics {}", arg[iarg + 10]);
      items.emplace_back(progbar);
      numobjs += 2 + progbar.tics;
      iarg += 11;
    } else {
      error->all(FLERR, iarg, "Unknown fix graphics/objects keyword {}", arg[iarg]);
    }
  }
  memory->create(imgobjs, numobjs, "fix_graphics/objects:imgobjs");
  memory->create(imgparms, numobjs, 10, "fix_graphics/objects:imgparms");
}
#undef PARSE_VARIABLE
/* ---------------------------------------------------------------------- */

FixGraphicsObjects::~FixGraphicsObjects()
{
  for (auto &gi : items) {
    switch (gi.style) {
      case SPHERE:
        delete[] gi.sphere.xstr;
        delete[] gi.sphere.ystr;
        delete[] gi.sphere.zstr;
        delete[] gi.sphere.dstr;
        break;
      case CYLINDER:
        delete[] gi.cylinder.x1str;
        delete[] gi.cylinder.y1str;
        delete[] gi.cylinder.z1str;
        delete[] gi.cylinder.x2str;
        delete[] gi.cylinder.y2str;
        delete[] gi.cylinder.z2str;
        delete[] gi.cylinder.dstr;
        break;
      case ARROW:
        delete[] gi.arrow.x1str;
        delete[] gi.arrow.y1str;
        delete[] gi.arrow.z1str;
        delete[] gi.arrow.x2str;
        delete[] gi.arrow.y2str;
        delete[] gi.arrow.z2str;
        delete[] gi.arrow.dstr;
        break;
      case CONE:
        delete[] gi.cone.x1str;
        delete[] gi.cone.y1str;
        delete[] gi.cone.z1str;
        delete[] gi.cone.x2str;
        delete[] gi.cone.y2str;
        delete[] gi.cone.z2str;
        delete[] gi.cone.d1str;
        delete[] gi.cone.d2str;
        break;
      case PROGBAR:
        delete[] gi.progbar.pstr;
        break;
      default:;    // do nothing
        break;
    }
  }

  memory->destroy(imgobjs);
  memory->destroy(imgparms);
}

/* ---------------------------------------------------------------------- */

int FixGraphicsObjects::setmask()
{
  return END_OF_STEP;
}

/* ---------------------------------------------------------------------- */

#define CHECK_VARIABLE(index, name)                                                     \
  if (name) {                                                                           \
    int ivar = input->variable->find(name);                                             \
    if (ivar < 0)                                                                       \
      error->all(FLERR, Error::NOLASTLINE,                                              \
                 "Variable name {} for fix graphics/objects does not exist", name);     \
    if (input->variable->equalstyle(ivar) == 0)                                         \
      error->all(FLERR, Error::NOLASTLINE,                                              \
                 "Fix graphics/objects variable {} is not equal-style variable", name); \
    index = ivar;                                                                       \
  }

void FixGraphicsObjects::init()
{
  int n = 0;
  for (auto &gi : items) {
    if (gi.style == SPHERE) {
      imgobjs[n] = Graphics::SPHERE;
      imgparms[n][0] = gi.sphere.type;
      CHECK_VARIABLE(gi.sphere.xvar, gi.sphere.xstr);
      CHECK_VARIABLE(gi.sphere.yvar, gi.sphere.ystr);
      CHECK_VARIABLE(gi.sphere.zvar, gi.sphere.zstr);
      CHECK_VARIABLE(gi.sphere.dvar, gi.sphere.dstr);
      ++n;
    } else if (gi.style == CYLINDER) {
      imgobjs[n] = Graphics::CYLINDER;
      imgparms[n][0] = gi.cylinder.type;
      CHECK_VARIABLE(gi.cylinder.x1var, gi.cylinder.x1str);
      CHECK_VARIABLE(gi.cylinder.y1var, gi.cylinder.y1str);
      CHECK_VARIABLE(gi.cylinder.z1var, gi.cylinder.z1str);
      CHECK_VARIABLE(gi.cylinder.x2var, gi.cylinder.x1str);
      CHECK_VARIABLE(gi.cylinder.y2var, gi.cylinder.y1str);
      CHECK_VARIABLE(gi.cylinder.z2var, gi.cylinder.z1str);
      CHECK_VARIABLE(gi.cylinder.dvar, gi.cylinder.dstr);
      ++n;
    } else if (gi.style == ARROW) {
      imgobjs[n] = Graphics::ARROW;
      imgparms[n][0] = gi.arrow.type;
      CHECK_VARIABLE(gi.arrow.x1var, gi.arrow.x1str);
      CHECK_VARIABLE(gi.arrow.y1var, gi.arrow.y1str);
      CHECK_VARIABLE(gi.arrow.z1var, gi.arrow.z1str);
      CHECK_VARIABLE(gi.arrow.x2var, gi.arrow.x2str);
      CHECK_VARIABLE(gi.arrow.y2var, gi.arrow.y2str);
      CHECK_VARIABLE(gi.arrow.z2var, gi.arrow.z2str);
      CHECK_VARIABLE(gi.arrow.dvar, gi.arrow.dstr);
      imgparms[n][9] = gi.arrow.ratio;
      ++n;
    } else if (gi.style == CONE) {
      imgobjs[n] = Graphics::CONE;
      imgparms[n][0] = gi.cone.type;
      CHECK_VARIABLE(gi.cone.x1var, gi.cone.x1str);
      CHECK_VARIABLE(gi.cone.y1var, gi.cone.y1str);
      CHECK_VARIABLE(gi.cone.z1var, gi.cone.z1str);
      CHECK_VARIABLE(gi.cone.x2var, gi.cone.x2str);
      CHECK_VARIABLE(gi.cone.y2var, gi.cone.y2str);
      CHECK_VARIABLE(gi.cone.z2var, gi.cone.z2str);
      CHECK_VARIABLE(gi.cone.d1var, gi.cone.d1str);
      CHECK_VARIABLE(gi.cone.d2var, gi.cone.d2str);
      ++n;
    } else if (gi.style == PROGBAR) {
      imgobjs[n] = Graphics::CYLINDER;
      imgparms[n][0] = gi.progbar.type1;
      imgparms[n][1] = gi.progbar.pos[X];
      imgparms[n][2] = gi.progbar.pos[Y];
      imgparms[n][3] = gi.progbar.pos[Z];
      imgparms[n][4] = gi.progbar.pos[X];
      imgparms[n][5] = gi.progbar.pos[Y];
      imgparms[n][6] = gi.progbar.pos[Z];
      imgparms[n][7] = gi.progbar.diameter;
      switch (gi.progbar.dim) {
        case X:
          imgparms[n][1] -= 0.5 * gi.progbar.length;
          imgparms[n][4] += 0.5 * gi.progbar.length;
          break;
        case Y:
          imgparms[n][2] -= 0.5 * gi.progbar.length;
          imgparms[n][5] += 0.5 * gi.progbar.length;
          break;
        case Z:
          imgparms[n][3] -= 0.5 * gi.progbar.length;
          imgparms[n][6] += 0.5 * gi.progbar.length;
          break;
        default:;    // do nothing
      }
      ++n;
      imgobjs[n] = Graphics::CYLINDER;
      imgparms[n][0] = gi.progbar.type2;
      imgparms[n][1] = gi.progbar.pos[X];
      imgparms[n][2] = gi.progbar.pos[Y];
      imgparms[n][3] = gi.progbar.pos[Z];
      imgparms[n][4] = gi.progbar.pos[X];
      imgparms[n][5] = gi.progbar.pos[Y];
      imgparms[n][6] = gi.progbar.pos[Z];
      imgparms[n][7] = 0.75 * gi.progbar.diameter;
      switch (gi.progbar.dim) {
        case X:
          imgparms[n][1] -= 0.5 * gi.progbar.length;
          imgparms[n][4] -= 0.5 * gi.progbar.length;
          imgparms[n][3] += 0.2 * gi.progbar.diameter;
          imgparms[n][6] += 0.2 * gi.progbar.diameter;
          break;
        case Y:
          imgparms[n][2] -= 0.5 * gi.progbar.length;
          imgparms[n][5] -= 0.5 * gi.progbar.length;
          imgparms[n][1] += 0.15 * gi.progbar.diameter;
          imgparms[n][4] += 0.15 * gi.progbar.diameter;
          break;
        case Z:
          imgparms[n][3] -= 0.5 * gi.progbar.length;
          imgparms[n][6] -= 0.5 * gi.progbar.length;
          imgparms[n][1] += 0.15 * gi.progbar.diameter;
          imgparms[n][4] += 0.15 * gi.progbar.diameter;
          break;
        default:
          break;
      }
      ++n;
      double delta = gi.progbar.length / (double) (gi.progbar.tics - 1);
      double lo = gi.progbar.pos[gi.progbar.dim] - 0.5 * gi.progbar.length;
      for (int i = 0; i < gi.progbar.tics; ++i) {
        imgobjs[n] = Graphics::CYLINDER;
        imgparms[n][0] = gi.progbar.type1;
        imgparms[n][1] = gi.progbar.pos[X];
        imgparms[n][2] = gi.progbar.pos[Y];
        imgparms[n][3] = gi.progbar.pos[Z];
        imgparms[n][4] = gi.progbar.pos[X];
        imgparms[n][5] = gi.progbar.pos[Y];
        imgparms[n][6] = gi.progbar.pos[Z];
        imgparms[n][7] = 1.1 * gi.progbar.diameter;
        switch (gi.progbar.dim) {
          case X:
            imgparms[n][1] = lo + delta * i - 0.05 * delta;
            imgparms[n][4] = lo + delta * i + 0.05 * delta;
            break;
          case Y:
            imgparms[n][2] = lo + delta * i - 0.05 * delta;
            imgparms[n][5] = lo + delta * i + 0.05 * delta;
            break;
          case Z:
            imgparms[n][3] = lo + delta * i - 0.05 * delta;
            imgparms[n][6] = lo + delta * i + 0.05 * delta;
            break;
        }
        ++n;
      }
      CHECK_VARIABLE(gi.progbar.pvar, gi.progbar.pstr);
    }
  }
  end_of_step();
}

#undef CHECK_VARIABLE
/* ---------------------------------------------------------------------- */

void FixGraphicsObjects::end_of_step()
{
  // evaluate variable if necessary, wrap with clear/add

  if (varflag) modify->clearstep_compute();

  int n = 0;
  for (auto &gi : items) {
    if (gi.style == SPHERE) {
      if (gi.sphere.xstr) gi.sphere.pos[X] = input->variable->compute_equal(gi.sphere.xvar);
      if (gi.sphere.ystr) gi.sphere.pos[Y] = input->variable->compute_equal(gi.sphere.yvar);
      if (gi.sphere.zstr) gi.sphere.pos[Z] = input->variable->compute_equal(gi.sphere.zvar);
      if (gi.sphere.dstr) gi.sphere.diameter = 2.0 * input->variable->compute_equal(gi.sphere.dvar);
      imgparms[n][1] = gi.sphere.pos[X];
      imgparms[n][2] = gi.sphere.pos[Y];
      imgparms[n][3] = gi.sphere.pos[Z];
      imgparms[n][4] = gi.sphere.diameter;
      ++n;
    } else if (gi.style == CYLINDER) {
      if (gi.cylinder.x1str)
        gi.cylinder.pos1[X] = input->variable->compute_equal(gi.cylinder.x1var);
      if (gi.cylinder.y1str)
        gi.cylinder.pos1[Y] = input->variable->compute_equal(gi.cylinder.y1var);
      if (gi.cylinder.z1str)
        gi.cylinder.pos1[Z] = input->variable->compute_equal(gi.cylinder.z1var);
      if (gi.cylinder.x2str)
        gi.cylinder.pos2[X] = input->variable->compute_equal(gi.cylinder.x2var);
      if (gi.cylinder.y2str)
        gi.cylinder.pos2[Y] = input->variable->compute_equal(gi.cylinder.y2var);
      if (gi.cylinder.z2str)
        gi.cylinder.pos2[Z] = input->variable->compute_equal(gi.cylinder.z2var);
      if (gi.cylinder.dstr)
        gi.cylinder.diameter = 2.0 * input->variable->compute_equal(gi.cylinder.dvar);
      imgparms[n][1] = gi.cylinder.pos1[X];
      imgparms[n][2] = gi.cylinder.pos1[Y];
      imgparms[n][3] = gi.cylinder.pos1[Z];
      imgparms[n][4] = gi.cylinder.pos2[X];
      imgparms[n][5] = gi.cylinder.pos2[Y];
      imgparms[n][6] = gi.cylinder.pos2[Z];
      imgparms[n][7] = gi.cylinder.diameter;
      ++n;
    } else if (gi.style == ARROW) {
      if (gi.arrow.x1str) gi.arrow.bot[X] = input->variable->compute_equal(gi.arrow.x1var);
      if (gi.arrow.y1str) gi.arrow.bot[Y] = input->variable->compute_equal(gi.arrow.y1var);
      if (gi.arrow.z1str) gi.arrow.bot[Z] = input->variable->compute_equal(gi.arrow.z1var);
      if (gi.arrow.x2str) gi.arrow.tip[X] = input->variable->compute_equal(gi.arrow.x2var);
      if (gi.arrow.y2str) gi.arrow.tip[Y] = input->variable->compute_equal(gi.arrow.y2var);
      if (gi.arrow.z2str) gi.arrow.tip[Z] = input->variable->compute_equal(gi.arrow.z2var);
      if (gi.arrow.dstr) gi.arrow.diameter = 2.0 * input->variable->compute_equal(gi.arrow.dvar);

      double mid[3], vec[3];
      MathExtra::add3(gi.arrow.tip, gi.arrow.bot, vec);
      MathExtra::scale3(0.5, vec, mid);
      MathExtra::sub3(gi.arrow.tip, gi.arrow.bot, vec);
      imgparms[n][1] = mid[X];
      imgparms[n][2] = mid[Y];
      imgparms[n][3] = mid[Z];
      imgparms[n][7] = MathExtra::len3(vec);
      MathExtra::norm3(vec);
      imgparms[n][4] = vec[X];
      imgparms[n][5] = vec[Y];
      imgparms[n][6] = vec[Z];
      imgparms[n][8] = gi.arrow.diameter;
      ++n;
    } else if (gi.style == CONE) {
      if (gi.cone.x1str) gi.cone.bot[X] = input->variable->compute_equal(gi.cone.x1var);
      if (gi.cone.y1str) gi.cone.bot[Y] = input->variable->compute_equal(gi.cone.y1var);
      if (gi.cone.z1str) gi.cone.bot[Z] = input->variable->compute_equal(gi.cone.z1var);
      if (gi.cone.x2str) gi.cone.top[X] = input->variable->compute_equal(gi.cone.x2var);
      if (gi.cone.y2str) gi.cone.top[Y] = input->variable->compute_equal(gi.cone.y2var);
      if (gi.cone.z2str) gi.cone.top[Z] = input->variable->compute_equal(gi.cone.z2var);
      if (gi.cone.d1str) gi.cone.botdiam = 2.0 * input->variable->compute_equal(gi.cone.d1var);
      if (gi.cone.d2str) gi.cone.topdiam = 2.0 * input->variable->compute_equal(gi.cone.d2var);

      imgparms[n][1] = gi.cone.bot[X];
      imgparms[n][2] = gi.cone.bot[Y];
      imgparms[n][3] = gi.cone.bot[Z];
      imgparms[n][4] = gi.cone.top[X];
      imgparms[n][5] = gi.cone.top[Y];
      imgparms[n][6] = gi.cone.top[Z];
      imgparms[n][7] = gi.cone.botdiam;
      imgparms[n][8] = gi.cone.topdiam;
      imgparms[n][9] = gi.cone.sides;
      ++n;
    } else if (gi.style == PROGBAR) {
      ++n;
      if (gi.progbar.pstr) gi.progbar.progress = input->variable->compute_equal(gi.progbar.pvar);
      // bracket into (0.0;1.0] rather than throwing an error for just a viz item
      gi.progbar.progress = std::max(std::min(gi.progbar.progress, 1.0), 1.0e-10);
      switch (gi.progbar.dim) {
        case X:
          imgparms[n][1] = gi.progbar.pos[X] + (gi.progbar.progress - 0.5) * gi.progbar.length;
          break;
        case Y:
          imgparms[n][2] = gi.progbar.pos[Y] + (gi.progbar.progress - 0.5) * gi.progbar.length;
          break;
        case Z:
          imgparms[n][3] = gi.progbar.pos[Z] + (gi.progbar.progress - 0.5) * gi.progbar.length;
          break;
        default:
          break;
      }
      ++n;
      n += gi.progbar.tics;
    }
  }
  if (varflag) modify->addstep_compute((update->ntimestep / nevery) * nevery + nevery);
}

/* ----------------------------------------------------------------------
   provide graphics information to dump image
------------------------------------------------------------------------- */

int FixGraphicsObjects::image(int *&objs, double **&parms)
{
  objs = imgobjs;
  parms = imgparms;
  return numobjs;
}
