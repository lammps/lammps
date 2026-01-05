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

#include "fix_graphics.h"

#include "comm.h"
#include "domain.h"
#include "dump_image.h"
#include "error.h"
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

enum { SPHERE, CYLINDER, ARROW, PROGBAR };
enum { X = 0, Y, Z };

/* ---------------------------------------------------------------------- */

FixGraphics::FixGraphics(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), imgobjs(nullptr), imgparms(nullptr)
{
  if (narg < 4) utils::missing_cmd_args(FLERR, "fix graphics", error);

  // parse mandatory arg

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  if (nevery <= 0) error->all(FLERR, 3, "Illegal fix graphics nevery value {}", nevery);
  global_freq = nevery;
  dynamic_group_allow = 1;

  numobjs = 0;
  varflag = 0;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "sphere") == 0) {
      if (iarg + 6 > narg) utils::missing_cmd_args(FLERR, "fix graphics sphere", error);
      // clang-format off
      SphereItem sphere{SPHERE, 1, {0.0, 0.0, 0.0}, 0.0, nullptr, nullptr, nullptr, nullptr,
                        -1, -1, -1, -1};
      // clang-format on
      sphere.type = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (strstr(arg[iarg + 2], "v_") == arg[iarg + 2]) {
        varflag = 1;
        sphere.xstr = utils::strdup(arg[iarg + 2] + 2);
      } else
        sphere.pos[X] = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      if (strstr(arg[iarg + 3], "v_") == arg[iarg + 3]) {
        varflag = 1;
        sphere.ystr = utils::strdup(arg[iarg + 3] + 2);
      } else
        sphere.pos[Y] = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      if (strstr(arg[iarg + 4], "v_") == arg[iarg + 4]) {
        varflag = 1;
        sphere.zstr = utils::strdup(arg[iarg + 4] + 2);
      } else
        sphere.pos[Z] = utils::numeric(FLERR, arg[iarg + 4], false, lmp);
      if (strstr(arg[iarg + 5], "v_") == arg[iarg + 5]) {
        varflag = 1;
        sphere.dstr = utils::strdup(arg[iarg + 5] + 2);
      } else
        sphere.diameter = 2.0 * utils::numeric(FLERR, arg[iarg + 5], false, lmp);
      items.emplace_back(std::move(sphere));
      ++numobjs;
      iarg += 6;
    } else if (strcmp(arg[iarg], "cylinder") == 0) {
      if (iarg + 9 > narg) utils::missing_cmd_args(FLERR, "fix graphics cylinder", error);
      // clang-format off
      CylinderItem cylinder{CYLINDER, 1, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0.0,
                            nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                            -1, -1, -1, -1, -1, -1, -1};
      // clang-format on
      cylinder.type = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (strstr(arg[iarg + 2], "v_") == arg[iarg + 2]) {
        varflag = 1;
        cylinder.x1str = utils::strdup(arg[iarg + 2] + 2);
      } else
        cylinder.pos1[X] = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      if (strstr(arg[iarg + 3], "v_") == arg[iarg + 3]) {
        varflag = 1;
        cylinder.y1str = utils::strdup(arg[iarg + 3] + 2);
      } else
        cylinder.pos1[Y] = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      if (strstr(arg[iarg + 4], "v_") == arg[iarg + 4]) {
        varflag = 1;
        cylinder.z1str = utils::strdup(arg[iarg + 4] + 2);
      } else
        cylinder.pos1[Z] = utils::numeric(FLERR, arg[iarg + 4], false, lmp);
      if (strstr(arg[iarg + 5], "v_") == arg[iarg + 5]) {
        varflag = 1;
        cylinder.x2str = utils::strdup(arg[iarg + 5] + 2);
      } else
        cylinder.pos2[X] = utils::numeric(FLERR, arg[iarg + 5], false, lmp);
      if (strstr(arg[iarg + 6], "v_") == arg[iarg + 6]) {
        varflag = 1;
        cylinder.y2str = utils::strdup(arg[iarg + 6] + 2);
      } else
        cylinder.pos2[Y] = utils::numeric(FLERR, arg[iarg + 6], false, lmp);
      if (strstr(arg[iarg + 7], "v_") == arg[iarg + 7]) {
        varflag = 1;
        cylinder.z2str = utils::strdup(arg[iarg + 7] + 2);
      } else
        cylinder.pos2[Z] = utils::numeric(FLERR, arg[iarg + 7], false, lmp);
      if (strstr(arg[iarg + 8], "v_") == arg[iarg + 8]) {
        varflag = 1;
        cylinder.dstr = utils::strdup(arg[iarg + 8] + 2);
      } else
        cylinder.diameter = 2.0 * utils::numeric(FLERR, arg[iarg + 8], false, lmp);
      items.emplace_back(std::move(cylinder));
      ++numobjs;
      iarg += 9;
    } else if (strcmp(arg[iarg], "arrow") == 0) {
      if (iarg + 10 > narg) utils::missing_cmd_args(FLERR, "fix graphics arrow", error);
      // clang-format off
      ArrowItem arrow{ARROW, 1, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 0.0, 0.1,
                            nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                            -1, -1, -1, -1, -1, -1, -1};
      // clang-format on
      arrow.type = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (strstr(arg[iarg + 2], "v_") == arg[iarg + 2]) {
        varflag = 1;
        arrow.x1str = utils::strdup(arg[iarg + 2] + 2);
      } else
        arrow.tip[X] = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      if (strstr(arg[iarg + 3], "v_") == arg[iarg + 3]) {
        varflag = 1;
        arrow.y1str = utils::strdup(arg[iarg + 3] + 2);
      } else
        arrow.tip[Y] = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      if (strstr(arg[iarg + 4], "v_") == arg[iarg + 4]) {
        varflag = 1;
        arrow.z1str = utils::strdup(arg[iarg + 4] + 2);
      } else
        arrow.tip[Z] = utils::numeric(FLERR, arg[iarg + 4], false, lmp);
      if (strstr(arg[iarg + 5], "v_") == arg[iarg + 5]) {
        varflag = 1;
        arrow.x2str = utils::strdup(arg[iarg + 5] + 2);
      } else
        arrow.bot[X] = utils::numeric(FLERR, arg[iarg + 5], false, lmp);
      if (strstr(arg[iarg + 6], "v_") == arg[iarg + 6]) {
        varflag = 1;
        arrow.y2str = utils::strdup(arg[iarg + 6] + 2);
      } else
        arrow.bot[Y] = utils::numeric(FLERR, arg[iarg + 6], false, lmp);
      if (strstr(arg[iarg + 7], "v_") == arg[iarg + 7]) {
        varflag = 1;
        arrow.z2str = utils::strdup(arg[iarg + 7] + 2);
      } else
        arrow.bot[Z] = utils::numeric(FLERR, arg[iarg + 7], false, lmp);
      if (strstr(arg[iarg + 8], "v_") == arg[iarg + 8]) {
        varflag = 1;
        arrow.dstr = utils::strdup(arg[iarg + 8] + 2);
      } else
        arrow.diameter = 2.0 * utils::numeric(FLERR, arg[iarg + 8], false, lmp);
      arrow.ratio = utils::numeric(FLERR, arg[iarg + 9], false, lmp);
      if ((arrow.ratio < 0.1) || (arrow.ratio > 0.5))
        error->all(FLERR, iarg + 9, "Arrow tip ratio must be between 0.1 and 0.5");
      items.emplace_back(std::move(arrow));
      numobjs += 2;
      iarg += 10;
    } else if (strcmp(arg[iarg], "progbar") == 0) {
      if (iarg + 11 > narg) utils::missing_cmd_args(FLERR, "fix graphics progbar", error);
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
      if (strstr(arg[iarg + 9], "v_") == arg[iarg + 9]) {
        varflag = 1;
        progbar.pstr = utils::strdup(arg[iarg + 9] + 2);
      } else {
        progbar.progress = utils::numeric(FLERR, arg[iarg + 9], false, lmp);
      }
      progbar.tics = utils::inumeric(FLERR, arg[iarg + 10], false, lmp);
      if ((progbar.tics < 0) || (progbar.tics > 20))
        error->all(FLERR, iarg + 10, "Unsupported number of progress bar tics {}", arg[iarg + 10]);
      items.emplace_back(std::move(progbar));
      numobjs += 2 + progbar.tics;
      iarg += 11;
    } else {
      error->all(FLERR, iarg, "Unknown fix graphics keyword {}", arg[iarg]);
    }
  }
  memory->create(imgobjs, numobjs, "fix_graphics:imgobjs");
  memory->create(imgparms, numobjs, 8, "fix_graphics:imgparms");
}

/* ---------------------------------------------------------------------- */

FixGraphics::~FixGraphics()
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

int FixGraphics::setmask()
{
  return END_OF_STEP;
}

/* ---------------------------------------------------------------------- */

void FixGraphics::init()
{
  int n = 0;
  for (auto &gi : items) {
    if (gi.style == SPHERE) {
      imgobjs[n] = DumpImage::SPHERE;
      imgparms[n][0] = gi.sphere.type;
      if (gi.sphere.xstr) {
        int ivar = input->variable->find(gi.sphere.xstr);
        if (ivar < 0)
          error->all(FLERR, Error::NOLASTLINE, "Variable name {} for fix graphics does not exist",
                     gi.sphere.xstr);
        if (input->variable->equalstyle(ivar) == 0)
          error->all(FLERR, Error::NOLASTLINE,
                     "Fix graphics variable {} is not equal-style variable", gi.sphere.xstr);
        gi.sphere.xvar = ivar;
      }
      if (gi.sphere.ystr) {
        int ivar = input->variable->find(gi.sphere.ystr);
        if (ivar < 0)
          error->all(FLERR, Error::NOLASTLINE, "Variable name {} for fix graphics does not exist",
                     gi.sphere.ystr);
        if (input->variable->equalstyle(ivar) == 0)
          error->all(FLERR, Error::NOLASTLINE,
                     "Fix graphics variable {} is not equal-style variable", gi.sphere.ystr);
        gi.sphere.yvar = ivar;
      }
      if (gi.sphere.zstr) {
        int ivar = input->variable->find(gi.sphere.zstr);
        if (ivar < 0)
          error->all(FLERR, Error::NOLASTLINE, "Variable name {} for fix graphics does not exist",
                     gi.sphere.zstr);
        if (input->variable->equalstyle(ivar) == 0)
          error->all(FLERR, Error::NOLASTLINE,
                     "Fix graphics variable {} is not equal-style variable", gi.sphere.zstr);
        gi.sphere.zvar = ivar;
      }
      if (gi.sphere.dstr) {
        int ivar = input->variable->find(gi.sphere.dstr);
        if (ivar < 0)
          error->all(FLERR, Error::NOLASTLINE, "Variable name {} for fix graphics does not exist",
                     gi.sphere.dstr);
        if (input->variable->equalstyle(ivar) == 0)
          error->all(FLERR, Error::NOLASTLINE,
                     "Fix graphics variable {} is not equal-style variable", gi.sphere.dstr);
        gi.sphere.dvar = ivar;
      }
      ++n;
    } else if (gi.style == CYLINDER) {
      imgobjs[n] = DumpImage::CYLINDER;
      imgparms[n][0] = gi.cylinder.type;
      if (gi.cylinder.x1str) {
        int ivar = input->variable->find(gi.cylinder.x1str);
        if (ivar < 0)
          error->all(FLERR, Error::NOLASTLINE, "Variable name {} for fix graphics does not exist",
                     gi.cylinder.x1str);
        if (input->variable->equalstyle(ivar) == 0)
          error->all(FLERR, Error::NOLASTLINE,
                     "Fix graphics variable {} is not equal-style variable", gi.cylinder.x1str);
        gi.cylinder.x1var = ivar;
      }
      if (gi.cylinder.y1str) {
        int ivar = input->variable->find(gi.cylinder.y1str);
        if (ivar < 0)
          error->all(FLERR, Error::NOLASTLINE, "Variable name {} for fix graphics does not exist",
                     gi.cylinder.y1str);
        if (input->variable->equalstyle(ivar) == 0)
          error->all(FLERR, Error::NOLASTLINE,
                     "Fix graphics variable {} is not equal-style variable", gi.cylinder.y1str);
        gi.cylinder.y1var = ivar;
      }
      if (gi.cylinder.z1str) {
        int ivar = input->variable->find(gi.cylinder.z1str);
        if (ivar < 0)
          error->all(FLERR, Error::NOLASTLINE, "Variable name {} for fix graphics does not exist",
                     gi.cylinder.z1str);
        if (input->variable->equalstyle(ivar) == 0)
          error->all(FLERR, Error::NOLASTLINE,
                     "Fix graphics variable {} is not equal-style variable", gi.cylinder.z1str);
        gi.cylinder.z1var = ivar;
      }
      if (gi.cylinder.x2str) {
        int ivar = input->variable->find(gi.cylinder.x2str);
        if (ivar < 0)
          error->all(FLERR, Error::NOLASTLINE, "Variable name {} for fix graphics does not exist",
                     gi.cylinder.x2str);
        if (input->variable->equalstyle(ivar) == 0)
          error->all(FLERR, Error::NOLASTLINE,
                     "Fix graphics variable {} is not equal-style variable", gi.cylinder.x2str);
        gi.cylinder.x2var = ivar;
      }
      if (gi.cylinder.y2str) {
        int ivar = input->variable->find(gi.cylinder.y2str);
        if (ivar < 0)
          error->all(FLERR, Error::NOLASTLINE, "Variable name {} for fix graphics does not exist",
                     gi.cylinder.y2str);
        if (input->variable->equalstyle(ivar) == 0)
          error->all(FLERR, Error::NOLASTLINE,
                     "Fix graphics variable {} is not equal-style variable", gi.cylinder.y2str);
        gi.cylinder.y2var = ivar;
      }
      if (gi.cylinder.z2str) {
        int ivar = input->variable->find(gi.cylinder.z2str);
        if (ivar < 0)
          error->all(FLERR, Error::NOLASTLINE, "Variable name {} for fix graphics does not exist",
                     gi.cylinder.z2str);
        if (input->variable->equalstyle(ivar) == 0)
          error->all(FLERR, Error::NOLASTLINE,
                     "Fix graphics variable {} is not equal-style variable", gi.cylinder.z2str);
        gi.cylinder.z2var = ivar;
      }
      if (gi.cylinder.dstr) {
        int ivar = input->variable->find(gi.cylinder.dstr);
        if (ivar < 0)
          error->all(FLERR, Error::NOLASTLINE, "Variable name {} for fix graphics does not exist",
                     gi.cylinder.dstr);
        if (input->variable->equalstyle(ivar) == 0)
          error->all(FLERR, Error::NOLASTLINE,
                     "Fix graphics variable {} is not equal-style variable", gi.cylinder.dstr);
        gi.cylinder.dvar = ivar;
      }
      ++n;
    } else if (gi.style == ARROW) {
      imgobjs[n] = DumpImage::CYLINDER;
      imgparms[n][0] = gi.arrow.type;
      if (gi.arrow.x1str) {
        int ivar = input->variable->find(gi.arrow.x1str);
        if (ivar < 0)
          error->all(FLERR, Error::NOLASTLINE, "Variable name {} for fix graphics does not exist",
                     gi.arrow.x1str);
        if (input->variable->equalstyle(ivar) == 0)
          error->all(FLERR, Error::NOLASTLINE,
                     "Fix graphics variable {} is not equal-style variable", gi.arrow.x1str);
        gi.arrow.x1var = ivar;
      }
      if (gi.arrow.y1str) {
        int ivar = input->variable->find(gi.arrow.y1str);
        if (ivar < 0)
          error->all(FLERR, Error::NOLASTLINE, "Variable name {} for fix graphics does not exist",
                     gi.arrow.y1str);
        if (input->variable->equalstyle(ivar) == 0)
          error->all(FLERR, Error::NOLASTLINE,
                     "Fix graphics variable {} is not equal-style variable", gi.arrow.y1str);
        gi.arrow.y1var = ivar;
      }
      if (gi.arrow.z1str) {
        int ivar = input->variable->find(gi.arrow.z1str);
        if (ivar < 0)
          error->all(FLERR, Error::NOLASTLINE, "Variable name {} for fix graphics does not exist",
                     gi.arrow.z1str);
        if (input->variable->equalstyle(ivar) == 0)
          error->all(FLERR, Error::NOLASTLINE,
                     "Fix graphics variable {} is not equal-style variable", gi.arrow.z1str);
        gi.arrow.z1var = ivar;
      }
      if (gi.arrow.x2str) {
        int ivar = input->variable->find(gi.arrow.x2str);
        if (ivar < 0)
          error->all(FLERR, Error::NOLASTLINE, "Variable name {} for fix graphics does not exist",
                     gi.arrow.x2str);
        if (input->variable->equalstyle(ivar) == 0)
          error->all(FLERR, Error::NOLASTLINE,
                     "Fix graphics variable {} is not equal-style variable", gi.arrow.x2str);
        gi.arrow.x2var = ivar;
      }
      if (gi.arrow.y2str) {
        int ivar = input->variable->find(gi.arrow.y2str);
        if (ivar < 0)
          error->all(FLERR, Error::NOLASTLINE, "Variable name {} for fix graphics does not exist",
                     gi.arrow.y2str);
        if (input->variable->equalstyle(ivar) == 0)
          error->all(FLERR, Error::NOLASTLINE,
                     "Fix graphics variable {} is not equal-style variable", gi.arrow.y2str);
        gi.arrow.y2var = ivar;
      }
      if (gi.arrow.z2str) {
        int ivar = input->variable->find(gi.arrow.z2str);
        if (ivar < 0)
          error->all(FLERR, Error::NOLASTLINE, "Variable name {} for fix graphics does not exist",
                     gi.arrow.z2str);
        if (input->variable->equalstyle(ivar) == 0)
          error->all(FLERR, Error::NOLASTLINE,
                     "Fix graphics variable {} is not equal-style variable", gi.arrow.z2str);
        gi.arrow.z2var = ivar;
      }
      if (gi.arrow.dstr) {
        int ivar = input->variable->find(gi.arrow.dstr);
        if (ivar < 0)
          error->all(FLERR, Error::NOLASTLINE, "Variable name {} for fix graphics does not exist",
                     gi.arrow.dstr);
        if (input->variable->equalstyle(ivar) == 0)
          error->all(FLERR, Error::NOLASTLINE,
                     "Fix graphics variable {} is not equal-style variable", gi.arrow.dstr);
        gi.arrow.dvar = ivar;
      }
      ++n;
      imgobjs[n] = DumpImage::CYLINDER;
      imgparms[n][0] = gi.arrow.type;
      ++n;
    } else if (gi.style == PROGBAR) {
      imgobjs[n] = DumpImage::CYLINDER;
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
      imgobjs[n] = DumpImage::CYLINDER;
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
        imgobjs[n] = DumpImage::CYLINDER;
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
      if (gi.progbar.pstr) {
        int ivar = input->variable->find(gi.progbar.pstr);
        if (ivar < 0)
          error->all(FLERR, Error::NOLASTLINE, "Variable name {} for fix graphics does not exist",
                     gi.progbar.pstr);
        if (input->variable->equalstyle(ivar) == 0)
          error->all(FLERR, Error::NOLASTLINE,
                     "Fix graphics variable {} is not equal-style variable", gi.progbar.pstr);
        gi.progbar.pvar = ivar;
      }
    }
  }
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixGraphics::end_of_step()
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
      if (gi.arrow.x1str) gi.arrow.tip[X] = input->variable->compute_equal(gi.arrow.x1var);
      if (gi.arrow.y1str) gi.arrow.tip[Y] = input->variable->compute_equal(gi.arrow.y1var);
      if (gi.arrow.z1str) gi.arrow.tip[Z] = input->variable->compute_equal(gi.arrow.z1var);
      if (gi.arrow.x2str) gi.arrow.bot[X] = input->variable->compute_equal(gi.arrow.x2var);
      if (gi.arrow.y2str) gi.arrow.bot[Y] = input->variable->compute_equal(gi.arrow.y2var);
      if (gi.arrow.z2str) gi.arrow.bot[Z] = input->variable->compute_equal(gi.arrow.z2var);
      if (gi.arrow.dstr) gi.arrow.diameter = 2.0 * input->variable->compute_equal(gi.arrow.dvar);

      double mid[3], vec[3];
      MathExtra::sub3(gi.arrow.bot, gi.arrow.tip, vec);
      MathExtra::scaleadd3(gi.arrow.ratio, vec, gi.arrow.tip, mid);
      imgparms[n][1] = gi.arrow.tip[X];
      imgparms[n][2] = gi.arrow.tip[Y];
      imgparms[n][3] = gi.arrow.tip[Z];
      imgparms[n][4] = mid[X];
      imgparms[n][5] = mid[Y];
      imgparms[n][6] = mid[Z];
      imgparms[n][7] = gi.arrow.diameter * (1.0 + 5.0 * gi.arrow.ratio);
      ++n;
      imgparms[n][1] = mid[X];
      imgparms[n][2] = mid[Y];
      imgparms[n][3] = mid[Z];
      imgparms[n][4] = gi.arrow.bot[X];
      imgparms[n][5] = gi.arrow.bot[Y];
      imgparms[n][6] = gi.arrow.bot[Z];
      imgparms[n][7] = gi.arrow.diameter;
      ++n;
    } else if (gi.style == PROGBAR) {
      ++n;
      if (gi.progbar.pstr) gi.progbar.progress = input->variable->compute_equal(gi.progbar.pvar);
      // bracket into [0.0;1.0] rather than throwing an error for just a viz item
      gi.progbar.progress = std::max(std::min(gi.progbar.progress, 1.0), 0.0);
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

int FixGraphics::image(int *&objs, double **&parms)
{
  objs = imgobjs;
  parms = imgparms;
  return numobjs;
}
