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

#include "dump_image.h"

#include "arg_info.h"
#include "atom.h"
#include "atom_vec.h"
#include "atom_vec_body.h"
#include "atom_vec_ellipsoid.h"
#include "atom_vec_line.h"
#include "atom_vec_tri.h"
#include "body.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "graphics.h"
#include "grid2d.h"
#include "grid3d.h"
#include "image.h"
#include "image_objects.h"
#include "input.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "molecule.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "output.h"
#include "platform.h"
#include "random_mars.h"
#include "region.h"
#include "region_block.h"
#include "region_cone.h"
#include "region_cylinder.h"
#include "region_ellipsoid.h"
#include "region_plane.h"
#include "region_prism.h"
#include "region_sphere.h"
#include "thermo.h"
#include "tokenizer.h"
#include "update.h"
#include "variable.h"

#include <cctype>
#include <cmath>
#include <cstring>

// clang-format on

using namespace LAMMPS_NS;
using MathConst::DEG2RAD;
using namespace ImageObjects;

namespace {
constexpr double BIG = 1.0e20;
enum { NUMERIC, ATOM, TYPE, ELEMENT, ATTRIBUTE, CONSTANT, INDEX };
enum { STATIC, DYNAMIC };
enum { NO = 0, YES = 1, AUTO = 2 };
enum { FILLED, FRAME, POINTS, TRANSPARENT };
enum { OFF = 0, CENTER, LOWERLEFT, LOWERRIGHT, UPPERLEFT, UPPERRIGHT };

//  convenience functions to change and restore lighting, assuming uncolored light

struct savedColors {
  double ambient;
  double key;
  double fill;
  double back;
};

savedColors reset_lighting(Image *image, double ambient, double key, double fill, double back)
{
  savedColors saved;
  saved.ambient = image->ambientColor[0];
  image->ambientColor[0] = image->ambientColor[1] = image->ambientColor[2] = ambient;
  saved.key = image->keyLightColor[0];
  image->keyLightColor[0] = image->keyLightColor[1] = image->keyLightColor[2] = key;
  saved.fill = image->fillLightColor[0];
  image->fillLightColor[0] = image->fillLightColor[1] = image->fillLightColor[2] = fill;
  saved.back = image->backLightColor[0];
  image->backLightColor[0] = image->backLightColor[1] = image->backLightColor[2] = back;
  return saved;
}

void restore_lighting(const savedColors &saved, Image *image)
{
  image->ambientColor[0] = image->ambientColor[1] = image->ambientColor[2] = saved.ambient;
  image->keyLightColor[0] = image->keyLightColor[1] = image->keyLightColor[2] = saved.key;
  image->fillLightColor[0] = image->fillLightColor[1] = image->fillLightColor[2] = saved.fill;
  image->backLightColor[0] = image->backLightColor[1] = image->backLightColor[2] = saved.back;
}

}    // namespace
// clang-format off

/* ---------------------------------------------------------------------- */

DumpImage::DumpImage(LAMMPS *lmp, int narg, char **arg) :
    DumpCustom(lmp, narg, arg), thetastr(nullptr), phistr(nullptr), cxstr(nullptr), cystr(nullptr),
    czstr(nullptr), upxstr(nullptr), upystr(nullptr), upzstr(nullptr), zoomstr(nullptr),
    diamtype(nullptr), diamelement(nullptr), bdiamtype(nullptr), colortype(nullptr),
    colorelement(nullptr), bcolortype(nullptr), aopacity(nullptr), bopacity(nullptr),
    grid2d(nullptr), grid3d(nullptr), id_grid_compute(nullptr), id_grid_fix(nullptr),
    grid_compute(nullptr), grid_fix(nullptr), gbuf(nullptr), avec_line(nullptr), avec_tri(nullptr),
    avec_ellipsoid(nullptr), avec_body(nullptr), image(nullptr), chooseghost(nullptr),
    bufcopy(nullptr)
{
  if (binary || multiproc) error->all(FLERR, 4, "Invalid dump image filename {}", filename);

  // force binary flag on to avoid corrupted output on Windows

  binary = 1;
  multifile_override = 0;

  // flag has_id as true to avoid bogus warnings about atom IDs for dump styles derived from DumpCustom

  has_id = true;

  // set filetype based on filename suffix

  if (utils::strmatch(filename, R"(\.jpg$)") || utils::strmatch(filename, R"(\.JPG$)") ||
      utils::strmatch(filename, R"(\.jpeg$)") || utils::strmatch(filename, R"(\.JPEG$)"))
    filetype = JPG;
  else if (utils::strmatch(filename, R"(\.png$)") || utils::strmatch(filename, R"(\.PNG$)"))
    filetype = PNG;
  else if (utils::strmatch(filename, R"(\.tga$)") || utils::strmatch(filename, R"(\.TGA$)"))
    filetype = TGA;
  else if (compressed && (utils::strmatch(filename, R"(\.jpg\.\w+$)") ||
                          utils::strmatch(filename, R"(\.JPG\.\w+$)") ||
                          utils::strmatch(filename, R"(\.jpeg\.\w+$)") ||
                          utils::strmatch(filename, R"(\.JPEG\.\w+$)")))
    error->all(FLERR, Error::NOLASTLINE, "Cannot use compression with JPEG images");
  else if (compressed && (utils::strmatch(filename, R"(\.png\.\w+$)") ||
                           utils::strmatch(filename, R"(\.PNG\.\w+$)")))
    error->all(FLERR, Error::NOLASTLINE, "Cannot use compression with PNG images");
  else if (compressed && (utils::strmatch(filename, R"(\.tga\.\w+$)") ||
                          utils::strmatch(filename, R"(\.TGA\.\w+$)")))
    error->all(FLERR, Error::NOLASTLINE, "Cannot use compression with TGA images");
  else filetype = PPM;

#ifndef LAMMPS_JPEG
  if (filetype == JPG)
    error->all(FLERR, Error::NOLASTLINE, "Support for writing images in JPEG format not included");
#endif
#ifndef LAMMPS_PNG
  if (filetype == PNG)
    error->all(FLERR, Error::NOLASTLINE, "Support for writing images in PNG format not included");
#endif

  // atom color,diameter settings

  if (nfield != 2)
    error->all(FLERR, Error::NOPOINTER,
               "Dump image command is missing attributes for color and size");

  acolor = ATTRIBUTE;
  if (strcmp(arg[5],"type") == 0) acolor = TYPE;
  else if (strcmp(arg[5],"element") == 0) acolor = ELEMENT;

  adiam = ATTRIBUTE;
  if (strcmp(arg[6],"type") == 0) adiam = TYPE;
  else if (strcmp(arg[6],"element") == 0) adiam = ELEMENT;

  // create Image class with two colormaps for atoms and grid cells
  // change defaults for 2d

  image = new Image(lmp,2);

  if (domain->dimension == 2) {
    image->theta = 0.0;
    image->phi = 0.0;
    image->up[0] = 0.0; image->up[1] = 1.0; image->up[2] = 0.0;
  }

  // set defaults for optional args

  atomflag = YES;
  gridflag = NO;
  lineflag = triflag = bodyflag = ellipsoidflag = NO;

  bcolor = ATOM;
  bdiam = NUMERIC;
  bdiamvalue = 0.5;
  if (atom->nbondtypes == 0) {
    bondflag = NO;
  } else {
    bondflag = YES;
  }

  cflag = STATIC;
  cx = cy = cz = 0.5;

  boxflag = YES;
  boxdiam = 0.02;
  axesflag = OFF;
  subboxflag = NO;
  boxopacity = 1.0;
  axesopacity = 1.0;
  subboxopacity = 1.0;

  // parse optional args

  int iarg = ioptional;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"atom") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR,"dump image atom", error);
      atomflag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;

    } else if (strcmp(arg[iarg],"adiam") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR,"dump image adiam", error);
      adiam = NUMERIC;
      adiamvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (adiamvalue <= 0.0)
        error->all(FLERR, iarg + 1, "Illegal dump image adiam value {}", adiamvalue);
      iarg += 2;

    } else if (strcmp(arg[iarg],"autobond") == 0) {
      if (iarg+3 > narg) utils::missing_cmd_args(FLERR,"dump image autobond", error);
      bondflag = AUTO;
      bcolor = ATOM;
      bondcutoff = utils::numeric(FLERR, arg[iarg+1],false, lmp);
      if (bondcutoff <= 0.0)
        error->all(FLERR, iarg + 1,"Illegal dump image autobond cutoff value {}", bondcutoff);
      bdiam = NUMERIC;
      bdiamvalue = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      if (bdiamvalue <= 0.0)
        error->all(FLERR, iarg + 2,"Illegal dump image autobond diameter value {}", bdiamvalue);
      iarg += 3;

    } else if (strcmp(arg[iarg],"bond") == 0) {
      if (iarg+3 > narg) utils::missing_cmd_args(FLERR,"dump image bond", error);
      if (atom->nbondtypes == 0)
        error->all(FLERR, iarg, "Dump image bond not allowed with no bond types defined");
      bondflag = YES;
      if (strcmp(arg[iarg+1],"none") == 0) bondflag = NO;
      else if (strcmp(arg[iarg+1],"atom") == 0) bcolor = ATOM;
      else if (strcmp(arg[iarg+1],"type") == 0) bcolor = TYPE;
      else error->all(FLERR, iarg + 1, "Unknown dump image bond color setting {}", arg[iarg + 1]);
      if (!islower(arg[iarg+2][0])) {
          bdiam = NUMERIC;
          bdiamvalue = utils::numeric(FLERR,arg[iarg+2],false,lmp);
          if (bdiamvalue <= 0.0)
            error->all(FLERR, iarg + 2,"Illegal dump image bond diameter value {}", bdiamvalue);
      } else if (strcmp(arg[iarg+2],"atom") == 0) bdiam = ATOM;
      else if (strcmp(arg[iarg+2],"type") == 0) bdiam = TYPE;
      else if (strcmp(arg[iarg+2],"none") == 0) bondflag = NO;
      else error->all(FLERR, iarg + 2, "Unknown dump image bond diameter setting {}", arg[iarg+2]);
      iarg += 3;

    } else if (strcmp(arg[iarg],"grid") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR,"dump image grid", error);
      gridflag = YES;

      char *id;
      int igrid,idata,index;
      int iflag = utils::check_grid_reference((char *) "Dump image", arg[iarg+1], nevery, id,
                                              igrid,idata,index,lmp);
      if (iflag < 0)
        error->all(FLERR, iarg+1, "Invalid grid reference {} in dump image command", arg[iarg+1]);

      if (iflag == ArgInfo::COMPUTE) {
        delete[] id_grid_compute;
        id_grid_compute = utils::strdup(id);
      } else if (iflag == ArgInfo::FIX) {
        delete[] id_grid_fix;
        id_grid_fix = utils::strdup(id);
      }
      delete[] id;
      grid_igrid = igrid;
      grid_idata = idata;
      grid_index = index;
      iarg += 2;

    } else if (strcmp(arg[iarg],"line") == 0) {
      if (iarg+3 > narg) utils::missing_cmd_args(FLERR,"dump image line", error);
      lineflag = YES;
      if (strcmp(arg[iarg+1],"type") == 0) lcolor = TYPE;
      else if (strcmp(arg[iarg+1],"atom") == 0) lcolor = ATOM;
      else if (strcmp(arg[iarg+1],"index") == 0) lcolor = INDEX;
      else
        error->all(FLERR, iarg+1, "Dump image line only supports color by type, atom, or index");
      ldiam = NUMERIC;
      ldiamvalue = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      iarg += 3;

    } else if (strcmp(arg[iarg],"tri") == 0) {
      if (iarg+4 > narg) utils::missing_cmd_args(FLERR,"dump image tri", error);
      triflag = YES;
      if (strcmp(arg[iarg+1],"type") == 0) tcolor = TYPE;
      else if (strcmp(arg[iarg+1],"atom") == 0) tcolor = ATOM;
      else if (strcmp(arg[iarg+1],"index") == 0) tcolor = INDEX;
      else
        error->all(FLERR, iarg+1, "Dump image tri only supports color by type, atom, or index");
      tstyle = utils::inumeric(FLERR,arg[iarg+2],false,lmp);
      tdiamvalue = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      iarg += 4;

    } else if (strcmp(arg[iarg],"ellipsoid") == 0) {
      if (iarg+5 > narg) utils::missing_cmd_args(FLERR,"dump image ellipsoid", error);
      ellipsoidflag = YES;
      if (strcmp(arg[iarg+1],"type") == 0) ecolor = TYPE;
      else if (strcmp(arg[iarg+1],"atom") == 0) ecolor = ATOM;
      else if (strcmp(arg[iarg+1],"index") == 0) ecolor = INDEX;
      else
        error->all(FLERR, iarg+1, "Dump image ellipsoid only supports color by type, atom, or index");
      estyle = utils::inumeric(FLERR,arg[iarg+2],false,lmp);
      // unless the wireframe setting is requested, we draw rounded triangles
      if (estyle !=2) estyle = 1;

      elevel = utils::inumeric(FLERR,arg[iarg+3],false,lmp);
      if (elevel == 0) elevel = 4; // default setting
      if ((elevel < 0) || (elevel > 6))
        error->all(FLERR, iarg+3, "Dump image ellipsoid invalid mesh refinement level {}", elevel);
      ediamvalue = utils::numeric(FLERR,arg[iarg+4],false,lmp);
      iarg += 5;

    } else if (strcmp(arg[iarg],"body") == 0) {
      if (iarg+4 > narg) utils::missing_cmd_args(FLERR,"dump image body", error);
      bodyflag = YES;
      if (strcmp(arg[iarg+1],"type") == 0) bodycolor = TYPE;
      else if (strcmp(arg[iarg+1],"atom") == 0) bodycolor = ATOM;
      else if (strcmp(arg[iarg+1],"index") == 0) bodycolor = INDEX;
      else
        error->all(FLERR, iarg+1, "Dump image body only supports color by type, atom, or index");
      bodyflag1 = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      bodyflag2 = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      iarg += 4;

    } else if (strcmp(arg[iarg],"compute") == 0) {
      if (iarg+5 > narg) utils::missing_cmd_args(FLERR,"dump image compute", error);
      std::string id_compute = arg[iarg+1];
      auto *computeptr = modify->get_compute_by_id(id_compute);
      if (!computeptr) error->all(FLERR, iarg+1, "Dump image compute ID {} does not exist", id_compute);
      int computecolor = TYPE;
      if (strcmp(arg[iarg+2],"type") == 0) computecolor = TYPE;
      else if (strcmp(arg[iarg+2],"element") == 0) computecolor = ELEMENT;
      else if (strcmp(arg[iarg+2],"const") == 0) computecolor = CONSTANT;
      else error->all(FLERR, iarg+2, "Unsupported color style for dump image compute {}", arg[iarg+2]);
      double computeflag1 = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      double computeflag2 = utils::numeric(FLERR,arg[iarg+4],false,lmp);
      objects.emplace_back(id_compute, computeptr, nullptr, computecolor, computeflag1, computeflag2,
                           image->color2rgb("white"));
      iarg += 5;

    } else if (strcmp(arg[iarg],"fix") == 0) {
      if (iarg+5 > narg) utils::missing_cmd_args(FLERR,"dump image fix", error);
      std::string id_fix = arg[iarg+1];
      auto *fixptr = modify->get_fix_by_id(id_fix);
      if (!fixptr) error->all(FLERR, iarg+1, "Dump image fix ID {} does not exist", id_fix);
      int fixcolor = TYPE;
      if (strcmp(arg[iarg+2],"type") == 0) fixcolor = TYPE;
      else if (strcmp(arg[iarg+2],"element") == 0) fixcolor = ELEMENT;
      else if (strcmp(arg[iarg+2],"const") == 0) fixcolor = CONSTANT;
      else error->all(FLERR, iarg+2, "Unsupported color style for dump image fix {}", arg[iarg+2]);
      double fixflag1 = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      double fixflag2 = utils::numeric(FLERR,arg[iarg+4],false,lmp);
      objects.emplace_back(id_fix, nullptr, fixptr, fixcolor, fixflag1, fixflag2,
                           image->color2rgb("white"));
      iarg += 5;

    } else if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+4 > narg) utils::missing_cmd_args(FLERR,"dump image region", error);
      auto *regptr = domain->get_region_by_id(arg[iarg+1]);
      if (!regptr)
        error->all(FLERR, iarg+1, "Unknown region {} for dump image", arg[iarg+1]);
      auto *regcolor = image->color2rgb(arg[iarg+2]);
      if (!regcolor)
        error->all(FLERR, iarg+2, "Unknown color {} for dump image", arg[iarg+2]);
      int drawstyle = FILLED;
      if (strcmp(arg[iarg+3],"filled") == 0) drawstyle = FILLED;
      else if (strcmp(arg[iarg+3],"frame") == 0) drawstyle = FRAME;
      else if (strcmp(arg[iarg+3],"points") == 0) drawstyle = POINTS;
      else if (strcmp(arg[iarg+3],"transparent") == 0) drawstyle = TRANSPARENT;
      else error->all(FLERR, iarg+3, "Unknown region draw style {}", arg[iarg+3]);
      double framediam = 0.5;
      double opacity = 1.0;
      int npoints = 0;
      if (drawstyle == FRAME) {
        if (iarg+5 > narg) utils::missing_cmd_args(FLERR,"dump image region frame", error);
        framediam = utils::numeric(FLERR, arg[iarg+4], false, lmp);
        if (framediam <= 0.0)
          error->all(FLERR, iarg+4, "Dump image region frame diameter must be > 0.0");
        ++iarg;
      } else if (drawstyle == POINTS) {
        if (iarg+6 > narg) utils::missing_cmd_args(FLERR,"dump image region points", error);
        npoints = utils::inumeric(FLERR, arg[iarg+4], false, lmp);
        if (npoints < 1)
          error->all(FLERR, iarg+4, "Dump image region number of points must be > 0");
        framediam = utils::numeric(FLERR, arg[iarg+5], false, lmp);
        if (framediam <= 0.0)
          error->all(FLERR, iarg+5, "Dump image region point diameter must be > 0.0");
        iarg += 2;
      } else if (drawstyle == TRANSPARENT) {
        if (iarg+5 > narg) utils::missing_cmd_args(FLERR,"dump image region transparent", error);
        opacity = utils::numeric(FLERR, arg[iarg+4], false, lmp);
        if ((opacity < 0.0) || (opacity > 1.0))
          error->all(FLERR, iarg+4, "Dump image region opacity must be in the range 0.0 to 1.0");

        ++iarg;
      }
      iarg += 4;
      regions.emplace_back(regptr->id, regptr, regcolor, drawstyle, framediam, opacity, npoints);

    } else if (strcmp(arg[iarg],"size") == 0) {
      if (iarg+3 > narg) utils::missing_cmd_args(FLERR,"dump image size", error);
      int width = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      int height = utils::inumeric(FLERR,arg[iarg+2],false,lmp);
      if (width <= 0 || height <= 0)
        error->all(FLERR, Error::NOPOINTER, "Illegal dump image dimensions");
      if (image->fsaa) {
        image->width = width*2;
        image->height = height*2;
      } else {
        image->width = width;
        image->height = height;
      }
      iarg += 3;

    } else if (strcmp(arg[iarg],"view") == 0) {
      if (iarg+3 > narg) utils::missing_cmd_args(FLERR,"dump image view", error);
      if (utils::strmatch(arg[iarg+1],"^v_")) {
        delete[] thetastr;
        thetastr = utils::strdup(arg[iarg+1]+2);
      } else {
        const double theta = utils::numeric(FLERR,arg[iarg+1],false,lmp);
        if (theta < 0.0 || theta > 180.0)
          error->all(FLERR, iarg+1, "Invalid dump image theta value");
        image->theta = DEG2RAD * theta;
      }
      if (utils::strmatch(arg[iarg+2],"^v_")) {
        delete[] phistr;
        phistr = utils::strdup(arg[iarg+2]+2);
      } else {
        image->phi = DEG2RAD * utils::numeric(FLERR,arg[iarg+2],false,lmp);
      }
      iarg += 3;

    } else if (strcmp(arg[iarg],"center") == 0) {
      if (iarg+5 > narg) utils::missing_cmd_args(FLERR,"dump image center", error);
      if (strcmp(arg[iarg+1],"s") == 0) cflag = STATIC;
      else if (strcmp(arg[iarg+1],"d") == 0) cflag = DYNAMIC;
      else error->all(FLERR, iarg+1, "Unknown dump image center style {}", arg[iarg+1]);
      if (utils::strmatch(arg[iarg+2],"^v_")) {
        delete[] cxstr;
        cxstr = utils::strdup(arg[iarg+2]+2);
        cflag = DYNAMIC;
      } else cx = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      if (utils::strmatch(arg[iarg+3],"^v_")) {
        delete[] cystr;
        cystr = utils::strdup(arg[iarg+3]+2);
        cflag = DYNAMIC;
      } else cy = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      if (utils::strmatch(arg[iarg+4],"^v_")) {
        delete[] czstr;
        czstr = utils::strdup(arg[iarg+4]+2);
        cflag = DYNAMIC;
      } else cz = utils::numeric(FLERR,arg[iarg+4],false,lmp);
      iarg += 5;

    } else if (strcmp(arg[iarg],"up") == 0) {
      if (iarg+4 > narg) utils::missing_cmd_args(FLERR,"dump image up", error);
      if (utils::strmatch(arg[iarg+1],"^v_")) {
        delete[] upxstr;
        upxstr = utils::strdup(arg[iarg+1]+2);
      } else image->up[0] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (utils::strmatch(arg[iarg+2],"^v_")) {
        delete[] upystr;
        upystr = utils::strdup(arg[iarg+2]+2);
      } else image->up[1] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      if (utils::strmatch(arg[iarg+3],"^v_")) {
        delete[] upzstr;
        upzstr = utils::strdup(arg[iarg+3]+2);
      } else image->up[2] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      iarg += 4;

    } else if (strcmp(arg[iarg],"zoom") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR,"dump image zoom", error);
      if (utils::strmatch(arg[iarg+1],"^v_")) {
        delete[] zoomstr;
        zoomstr = utils::strdup(arg[iarg+1]+2);
      } else {
        double zoom = utils::numeric(FLERR,arg[iarg+1],false,lmp);
        if (zoom <= 0.0) error->all(FLERR, iarg+1, "Invail dump image zoom value {}", zoom);
        image->zoom = zoom;
      }
      iarg += 2;

    } else if (strcmp(arg[iarg],"box") == 0) {
      if (iarg+3 > narg) utils::missing_cmd_args(FLERR,"dump image box", error);
      boxflag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      boxdiam = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      if (boxdiam < 0.0) error->all(FLERR, iarg+2, "Invalid dump image box diameter {}", boxdiam);
      iarg += 3;

    } else if (strcmp(arg[iarg],"axes") == 0) {
      if (iarg+3 > narg) utils::missing_cmd_args(FLERR,"dump image axes", error);
      if (strcmp(arg[iarg+1],"no") == 0) {
        axesflag = OFF;
      } else if (strcmp(arg[iarg+1],"center") == 0) {
        axesflag = CENTER;
      } else if (strcmp(arg[iarg+1],"lowerleft") == 0) {
        axesflag = LOWERLEFT;
      } else if (strcmp(arg[iarg+1],"yes") == 0) {
        axesflag = LOWERLEFT;
      } else if (strcmp(arg[iarg+1],"lowerright") == 0) {
        axesflag = LOWERRIGHT;
      } else if (strcmp(arg[iarg+1],"upperleft") == 0) {
        axesflag = UPPERLEFT;
      } else if (strcmp(arg[iarg+1],"upperright") == 0) {
        axesflag = UPPERRIGHT;
      } else {
        error->all(FLERR, iarg+1, "Unknown axes location {}", arg[iarg+1]);
      }

      if (axesflag) {
        axeslen = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        axesdiam = utils::numeric(FLERR,arg[iarg+3],false,lmp);
        if (axeslen <= 0.0)
          error->all(FLERR, iarg+2, "Invalid dump image axes length {}", axeslen);
        if (axesdiam <= 0.0)
          error->all(FLERR,"Invalid dump image axes diameter {}", axesdiam);
      }
      iarg += 4;

    } else if (strcmp(arg[iarg],"subbox") == 0) {
      if (iarg+3 > narg) utils::missing_cmd_args(FLERR,"dump image subbox", error);
      subboxflag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      subboxdiam = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      if (subboxdiam < 0.0)
        error->all(FLERR,"Invalid dump image subbox diameter {}", subboxdiam);
      iarg += 3;

    } else if (strcmp(arg[iarg],"shiny") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR,"dump image shiny", error);
      double shiny = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (shiny < 0.0 || shiny > 1.0)
        error->all(FLERR, iarg+1, "Invalid dump image shiny value {}", shiny);
      image->shiny = shiny;
      iarg += 2;

    } else if (strcmp(arg[iarg],"fsaa") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR,"dump image fsaa", error);
      int aa = utils::logical(FLERR, arg[iarg+1], false, lmp);
      if (aa) {
        if (!image->fsaa) {
          image->width = image->width*2;
          image->height = image->height*2;
        }
      } else {
        if (image->fsaa) {
          image->width = image->width/2;
          image->height = image->height/2;
        }
      }
      image->fsaa = aa;
      iarg += 2;

    } else if (strcmp(arg[iarg],"ssao") == 0) {
      if (iarg+4 > narg) utils::missing_cmd_args(FLERR,"dump image ssao", error);
      image->ssao = utils::logical(FLERR,arg[iarg+1],false,lmp);
      int seed = utils::inumeric(FLERR,arg[iarg+2],false,lmp);
      if (seed <= 0) error->all(FLERR, iarg + 2, "Invalid dump image ssao seed {}", seed);
      image->seed = seed;
      double ssaoint = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      if (ssaoint < 0.0 || ssaoint > 1.0)
        error->all(FLERR, iarg+3, "Invalid dump image ssao intensity value {}", ssaoint);
      image->ssaoint = ssaoint;
      iarg += 4;

    } else error->all(FLERR, iarg, "Unknown dump image keyword {}", arg[iarg]);
  }

  // error checks and setup for lineflag, triflag, bodyflag

  if (lineflag) {
    avec_line = dynamic_cast<AtomVecLine *>(atom->style_match("line"));
    if (!avec_line)
      error->all(FLERR, Error::NOLASTLINE, "Dump image line requires atom style line");
  }
  if (triflag) {
    avec_tri = dynamic_cast<AtomVecTri *>(atom->style_match("tri"));
    if (!avec_tri)
      error->all(FLERR, Error::NOLASTLINE, "Dump image tri requires atom style tri");
  }

  if (ellipsoidflag) {
    avec_ellipsoid = dynamic_cast<AtomVecEllipsoid *>(atom->style_match("ellipsoid"));
    if (!avec_ellipsoid)
      error->all(FLERR, Error::NOLASTLINE, "Dump image ellipsoid requires atom style ellipsoid");
  }
  if (bodyflag) {
    avec_body = dynamic_cast<AtomVecBody *>(atom->style_match("body"));
    if (!avec_body)
      error->all(FLERR, Error::NOLASTLINE, "Dump image body yes requires atom style body");
  }

  extraflag = 0;
  if (lineflag || triflag || ellipsoidflag || bodyflag) extraflag = 1;

  // allocate image buffer now that image size is known

  image->buffers();

  // communication needed for bonds colored by atoms

  if (bondflag) {
    if (bcolor == ATOM || bdiam == ATOM) comm_forward = 3;
    else comm_forward = 1;
  }

  // additional defaults for dump_modify options

  diamtype = new double[ntypes+1];
  diamelement = new double[ntypes+1];
  colortype = new double*[ntypes+1];
  colorelement = new double*[ntypes+1];
  aopacity = new double[ntypes+1];

  for (int i = 1; i <= ntypes; i++) {
    diamtype[i] = 1.0;
    aopacity[i] = 1.0;
    if (i % 6 == 1) colortype[i] = image->color2rgb("red");
    else if (i % 6 == 2) colortype[i] = image->color2rgb("green");
    else if (i % 6 == 3) colortype[i] = image->color2rgb("blue");
    else if (i % 6 == 4) colortype[i] = image->color2rgb("yellow");
    else if (i % 6 == 5) colortype[i] = image->color2rgb("cyan");
    else if (i % 6 == 0) colortype[i] = image->color2rgb("magenta");
  }

  if (bondflag == YES) {
    bdiamtype = new double[atom->nbondtypes+1];
    bcolortype = new double*[atom->nbondtypes+1];
    bopacity = new double[atom->nbondtypes+1];
    for (int i = 1; i <= atom->nbondtypes; i++) {
      bdiamtype[i] = 0.5;
      bopacity[i] = 1.0;
      if (i % 6 == 1) bcolortype[i] = image->color2rgb("red");
      else if (i % 6 == 2) bcolortype[i] = image->color2rgb("green");
      else if (i % 6 == 3) bcolortype[i] = image->color2rgb("blue");
      else if (i % 6 == 4) bcolortype[i] = image->color2rgb("yellow");
      else if (i % 6 == 5) bcolortype[i] = image->color2rgb("cyan");
      else if (i % 6 == 0) bcolortype[i] = image->color2rgb("magenta");
    }
  }

  // viewflag = DYNAMIC if any view parameter is dynamic

  viewflag = STATIC;
  if (thetastr || phistr || cflag == DYNAMIC ||
      upxstr || upystr || upzstr || zoomstr) viewflag = DYNAMIC;

  box_bounds();
  if (cflag == STATIC) box_center();
  if (viewflag == STATIC) view_params();

  // local data

  maxbufcopy = 0;
  maxgrid = 0;
}

/* ---------------------------------------------------------------------- */

DumpImage::~DumpImage()
{
  delete image;
  output->thermo->set_image_fname("");

  delete[] diamtype;
  delete[] diamelement;
  delete[] colortype;
  delete[] colorelement;
  delete[] aopacity;
  delete[] bdiamtype;
  delete[] bcolortype;
  delete[] bopacity;
  memory->destroy(chooseghost);
  memory->destroy(bufcopy);
  memory->destroy(gbuf);

  delete[] upxstr;
  delete[] upystr;
  delete[] upzstr;
  delete[] zoomstr;
  delete[] thetastr;
  delete[] phistr;
  delete[] cxstr;
  delete[] cystr;
  delete[] czstr;

  delete[] id_grid_compute;
  delete[] id_grid_fix;
}

/* ---------------------------------------------------------------------- */

void DumpImage::init_style()
{
  if (multifile == 0 && !multifile_override)
    error->all(FLERR, Error::NOLASTLINE,
               "Dump image requires file name with '*' requesting one snapshot per file");
  if (sort_flag) error->all(FLERR, Error::NOLASTLINE, "Dump image cannot perform sorting");

  DumpCustom::init_style();

  // for grid output, find current ptr for compute or fix
  // check that fix frequency is acceptable

  if (gridflag) {
    if (id_grid_compute) {
      grid_compute = modify->get_compute_by_id(id_grid_compute);
      if (!grid_compute)
        error->all(FLERR, Error::NOLASTLINE,
                   "Could not find dump image grid compute ID {}", id_grid_compute);
    } else if (id_grid_fix) {
      grid_fix = modify->get_fix_by_id(id_grid_fix);
      if (!grid_fix)
        error->all(FLERR, Error::NOLASTLINE, "Could not find dump image fix ID {}",id_grid_fix);
      if (nevery % grid_fix->peratom_freq)
        error->all(FLERR, Error::NOLASTLINE, "Dump image and grid fix {} not computed at "
                   "compatible times{}", id_grid_fix, utils::errorurl(7));
    }
  }

  // check image variables

  if (thetastr) {
    thetavar = input->variable->find(thetastr);
    if (thetavar < 0)
      error->all(FLERR, Error::NOLASTLINE, "Variable name for dump image theta does not exist");
    if (!input->variable->equalstyle(thetavar))
      error->all(FLERR, Error::NOLASTLINE, "Variable for dump image theta is invalid style");
  }
  if (phistr) {
    phivar = input->variable->find(phistr);
    if (phivar < 0)
      error->all(FLERR, Error::NOLASTLINE, "Variable name for dump image phi does not exist");
    if (!input->variable->equalstyle(phivar))
      error->all(FLERR, Error::NOLASTLINE, "Variable for dump image phi is invalid style");
  }
  if (cxstr) {
    cxvar = input->variable->find(cxstr);
    if (cxvar < 0)
      error->all(FLERR, Error::NOLASTLINE, "Variable name for dump image center does not exist");
    if (!input->variable->equalstyle(cxvar))
      error->all(FLERR, Error::NOLASTLINE, "Variable for dump image center is invalid style");
  }
  if (cystr) {
    cyvar = input->variable->find(cystr);
    if (cyvar < 0)
      error->all(FLERR, Error::NOLASTLINE, "Variable name for dump image center does not exist");
    if (!input->variable->equalstyle(cyvar))
      error->all(FLERR, Error::NOLASTLINE, "Variable for dump image center is invalid style");
  }
  if (czstr) {
    czvar = input->variable->find(czstr);
    if (czvar < 0)
      error->all(FLERR, Error::NOLASTLINE, "Variable name for dump image center does not exist");
    if (!input->variable->equalstyle(czvar))
      error->all(FLERR, Error::NOLASTLINE, "Variable for dump image center is invalid style");
  }
  if (upxstr) {
    upxvar = input->variable->find(upxstr);
    if (upxvar < 0)
      error->all(FLERR, Error::NOLASTLINE, "Variable name for dump image center does not exist");
    if (!input->variable->equalstyle(upxvar))
      error->all(FLERR, Error::NOLASTLINE, "Variable for dump image center is invalid style");
  }
  if (upystr) {
    upyvar = input->variable->find(upystr);
    if (upyvar < 0)
      error->all(FLERR, Error::NOLASTLINE, "Variable name for dump image center does not exist");
    if (!input->variable->equalstyle(upyvar))
      error->all(FLERR, Error::NOLASTLINE, "Variable for dump image center is invalid style");
  }
  if (upzstr) {
    upzvar = input->variable->find(upzstr);
    if (upzvar < 0)
      error->all(FLERR, Error::NOLASTLINE, "Variable name for dump image center does not exist");
    if (!input->variable->equalstyle(upzvar))
      error->all(FLERR, Error::NOLASTLINE, "Variable for dump image center is invalid style");
  }
  if (zoomstr) {
    zoomvar = input->variable->find(zoomstr);
    if (zoomvar < 0)
      error->all(FLERR, Error::NOLASTLINE, "Variable name for dump image zoom does not exist");
    if (!input->variable->equalstyle(zoomvar))
      error->all(FLERR, Error::NOLASTLINE, "Variable for dump image zoom is invalid style");
  }

  // set up type -> element mapping

  for (int i = 1; i <= ntypes; i++) {
    colorelement[i] = image->element2color(typenames[i]);
    if (colorelement[i] == nullptr)
      error->all(FLERR, Error::NOLASTLINE, "Invalid dump image element name");
  }

  for (int i = 1; i <= ntypes; i++) {
    diamelement[i] = image->element2diam(typenames[i]);
    if (diamelement[i] == 0.0)
      error->all(FLERR, Error::NOLASTLINE, "Invalid dump image element name");
  }

  if (bondflag == AUTO) {
    if (neighbor->style == Neighbor::MULTI)
      error->all(FLERR, Error::NOLASTLINE,
                 "Dump image autobond requires neighbor style 'bin' or 'nsq'");
    if (force->pair == nullptr)
      error->all(FLERR, "Dump image autobond requires a pair style to be defined");
    if (bondcutoff > neighbor->cutneighmax)
      error->all(FLERR, "Dump image autobond cutoff > neighbor cutoff");
    if ((bondcutoff > neighbor->cutneighmin) && (comm->me == 0))
      error->warning(FLERR, "Dump image autobond cutoff > minimum neighbor cutoff");
    if ((domain->xperiodic && (bondcutoff > domain->xprd)) ||
        (domain->yperiodic && (bondcutoff > domain->yprd)) ||
        ((domain->dimension == 3) && domain->zperiodic && (bondcutoff > domain->zprd)))
      error->all(FLERR, "Dump image autobond cutoff is larger than periodic domain");
  }

  // check if computes and fixes with visualization info still exist
  for (auto &iobj : objects) {
    if (iobj.cptr) {
      auto *computeptr = modify->get_compute_by_id(iobj.id);
      if (!computeptr)
        error->all(FLERR, Error::NOLASTLINE, "Compute ID {} for dump image does not exist", iobj.id);
      iobj.cptr = computeptr;

    } else if (iobj.fptr) {
      auto *fixptr = modify->get_fix_by_id(iobj.id);
      if (!fixptr)
        error->all(FLERR, Error::NOLASTLINE, "Fix ID {} for dump image does not exist", iobj.id);
      iobj.fptr = fixptr;

      // check if fix data for dump image is available at the required steps
      int nfreq = fixptr->global_freq;
      if ((update->ntimestep != 0) && ((nfreq == 0) || (nevery % nfreq)))
        error->all(FLERR, Error::NOLASTLINE,
                   "Dump {} and fix {} are not executed at compatible timesteps {}",
                   style, fixptr->style, utils::errorurl(7));
    }
  }
}

/* ---------------------------------------------------------------------- */

void DumpImage::write()
{
  // open new file

  openfile();

  // reset box center and view parameters if dynamic

  box_bounds();
  if (cflag == DYNAMIC) box_center();
  if (viewflag == DYNAMIC) view_params();

  // nme = # of atoms this proc will contribute to dump

  nme = count();

  if (nme > maxbuf) {
    maxbuf = nme;
    memory->destroy(buf);
    memory->create(buf,maxbuf*size_one,"dump:buf");
  }

  // pack atom buf with color & diameter

  pack(nullptr);

  // set minmax color range if using dynamic atom color map

  if (acolor == ATTRIBUTE && image->map_dynamic(0)) {
    double two[2],twoall[2];
    double lo = BIG;
    double hi = -BIG;
    int m = 0;
    for (int i = 0; i < nchoose; i++) {
      lo = MIN(lo,buf[m]);
      hi = MAX(hi,buf[m]);
      m += size_one;
    }
    two[0] = -lo;
    two[1] = hi;
    MPI_Allreduce(two,twoall,2,MPI_DOUBLE,MPI_MAX,world);
    int flag = image->map_minmax(0,-twoall[0],twoall[1]);
    if (flag) error->all(FLERR,"Invalid atom color map min/max values");
  }

  // pack grid gbuf with grid cell values
  // ngrid = # of grid cells this proc owns

  if (gridflag) {
    if (domain->dimension == 2) {
      if (grid_compute)
        grid2d = (Grid2d *) grid_compute->get_grid_by_index(grid_igrid);
      else if (grid_fix)
        grid2d = (Grid2d *) grid_fix->get_grid_by_index(grid_igrid);
      grid2d->get_size(nxgrid,nygrid);
      grid2d->get_bounds_owned(nxlo_in,nxhi_in,nylo_in,nyhi_in);
      ngrid = (nxhi_in-nxlo_in+1) * (nyhi_in-nylo_in+1);
    } else {
      if (grid_compute)
        grid3d = (Grid3d *) grid_compute->get_grid_by_index(grid_igrid);
      else if (grid_fix)
        grid3d = (Grid3d *) grid_fix->get_grid_by_index(grid_igrid);
      grid3d->get_size(nxgrid,nygrid,nzgrid);
      grid3d->get_bounds_owned(nxlo_in,nxhi_in,nylo_in,nyhi_in,nzlo_in,nzhi_in);
      ngrid = (nxhi_in-nxlo_in+1) * (nyhi_in-nylo_in+1) * (nzhi_in-nzlo_in+1);
    }

    // ensure gbuf is large enough

    if (ngrid > maxgrid) {
      memory->destroy(gbuf);
      maxgrid = ngrid;
      memory->create(gbuf,maxgrid,"dump/image:gbuf");
    }

    // invoke Compute for per-grid quantities
    // cannot invoke before first run, otherwise invoke if necessary

    if (grid_compute) {
      if (!grid_compute->is_initialized())
        error->all(FLERR,"Grid compute ID {} used in dump image cannot be invoked "
                   "before initialization by a run", grid_compute->id);
      if (!(grid_compute->invoked_flag & Compute::INVOKED_PERGRID)) {
        grid_compute->compute_pergrid();
        grid_compute->invoked_flag |= Compute::INVOKED_PERGRID;
      }
    }

    // access grid data and load gbuf

    if (domain->dimension == 2) {
      if (grid_index == 0) {
        double **vec2d;
        if (grid_compute)
          vec2d = (double **)
            grid_compute->get_griddata_by_index(grid_idata);
        else if (grid_fix)
          vec2d = (double **)
          grid_fix->get_griddata_by_index(grid_idata);
        int n = 0;
        for (int iy = nylo_in; iy <= nyhi_in; iy++)
          for (int ix = nxlo_in; ix <= nxhi_in; ix++)
            gbuf[n++] = vec2d[iy][ix];
      } else {
        double ***array2d;
        if (grid_compute)
          array2d = (double ***)
            grid_compute->get_griddata_by_index(grid_idata);
        else if (grid_fix)
          array2d = (double ***)
            grid_fix->get_griddata_by_index(grid_idata);
        int index = grid_index - 1;
        int n = 0;
        for (int iy = nylo_in; iy <= nyhi_in; iy++)
          for (int ix = nxlo_in; ix <= nxhi_in; ix++)
            gbuf[n++] = array2d[iy][ix][index];
      }

    } else if (domain->dimension == 3) {
      if (grid_index == 0) {
        double ***vec3d;
        if (grid_compute)
          vec3d = (double ***)
            grid_compute->get_griddata_by_index(grid_idata);
        else if (grid_fix)
          vec3d = (double ***)
            grid_fix->get_griddata_by_index(grid_idata);
        int n = 0;
        for (int iz = nzlo_in; iz <= nzhi_in; iz++)
          for (int iy = nylo_in; iy <= nyhi_in; iy++)
            for (int ix = nxlo_in; ix <= nxhi_in; ix++)
              gbuf[n++] = vec3d[iz][iy][ix];

        }
    } else {
      double ****array3d;
      if (grid_compute)
        array3d = (double ****)
          grid_compute->get_griddata_by_index(grid_idata);
      else if (grid_fix)
        array3d = (double ****)
          grid_fix->get_griddata_by_index(grid_idata);
      int index = grid_index - 1;
      int n = 0;
      for (int iz = nzlo_in; iz <= nzhi_in; iz++)
        for (int iy = nylo_in; iy <= nyhi_in; iy++)
          for (int ix = nxlo_in; ix <= nxhi_in; ix++)
            gbuf[n++] = array3d[iz][iy][ix][index];
      }
  }

  // set minmax color range if using dynamic grid color map

  if (gridflag && image->map_dynamic(1)) {
    double two[2],twoall[2];
    double lo = BIG;
    double hi = -BIG;
    for (int i = 0; i < ngrid; i++) {
      lo = MIN(lo,gbuf[i]);
      hi = MAX(hi,gbuf[i]);
    }
    two[0] = -lo;
    two[1] = hi;
    MPI_Allreduce(two,twoall,2,MPI_DOUBLE,MPI_MAX,world);
    int flag = image->map_minmax(1,-twoall[0],twoall[1]);
    if (flag) error->all(FLERR,"Invalid grid color map min/max values");
  }

  // create image on each proc, then merge them

  image->clear();
  create_image();
  image->merge();

  // write image file

  if (me == 0) {
    if (filetype == JPG) image->write_JPG(fp);
    else if (filetype == PNG) image->write_PNG(fp);
    else if (filetype == TGA) image->write_TGA(fp);
    else image->write_PPM(fp);
    if (multifile) {
      fp = nullptr;             // implicitly close file

      // cache last dump image filename for access through library interface.
      // update only *after* the file has been written so there will be no invalid read.
      // have to recreate the substitution done within openfile().

      output->thermo->set_image_fname(utils::star_subst(filename, update->ntimestep, padflag));
    }
  }
}

/* ----------------------------------------------------------------------
   simulation box bounds
------------------------------------------------------------------------- */

void DumpImage::box_bounds()
{
  if (domain->triclinic == 0) {
    boxxlo = domain->boxlo[0];
    boxxhi = domain->boxhi[0];
    boxylo = domain->boxlo[1];
    boxyhi = domain->boxhi[1];
    boxzlo = domain->boxlo[2];
    boxzhi = domain->boxhi[2];
  } else {
    boxxlo = domain->boxlo_bound[0];
    boxxhi = domain->boxhi_bound[0];
    boxylo = domain->boxlo_bound[1];
    boxyhi = domain->boxhi_bound[1];
    boxzlo = domain->boxlo_bound[2];
    boxzhi = domain->boxhi_bound[2];
    boxxy = domain->xy;
    boxxz = domain->xz;
    boxyz = domain->yz;
  }
}

/* ----------------------------------------------------------------------
   reset view parameters
   called once from constructor if view is STATIC
   called every snapshot from write() if view is DYNAMIC
------------------------------------------------------------------------- */

void DumpImage::box_center()
{
  if (cxstr) cx = input->variable->compute_equal(cxvar);
  if (cystr) cy = input->variable->compute_equal(cyvar);
  if (czstr) cz = input->variable->compute_equal(czvar);

  image->xctr = boxxlo + cx*(boxxhi-boxxlo);
  image->yctr = boxylo + cy*(boxyhi-boxylo);
  image->zctr = boxzlo + cz*(boxzhi-boxzlo);
}

/* ----------------------------------------------------------------------
   reset view parameters in Image class
   called once from constructor if view is STATIC
   called every snapshot from write() if view is DYNAMIC
------------------------------------------------------------------------- */

void DumpImage::view_params()
{
  // view direction theta and phi

  if (thetastr) {
    const double theta = input->variable->compute_equal(thetavar);
    if (theta < 0.0 || theta > 180.0)
      error->all(FLERR,"Invalid dump image theta value");
    image->theta = DEG2RAD * theta;
  }

  if (phistr) image->phi = DEG2RAD * input->variable->compute_equal(phivar);

  // up vector

  if (upxstr) image->up[0] = input->variable->compute_equal(upxvar);
  if (upystr) image->up[1] = input->variable->compute_equal(upyvar);
  if (upzstr) image->up[2] = input->variable->compute_equal(upzvar);

  // zoom

  if (zoomstr) image->zoom = input->variable->compute_equal(zoomvar);
  if (image->zoom <= 0.0) error->all(FLERR,"Invalid dump image zoom value");

  // remainder of view setup is internal to Image class

  image->view_params(boxxlo,boxxhi,boxylo,boxyhi,boxzlo,boxzhi);
}

/* ----------------------------------------------------------------------
   create image for all data this proc owns
   all procs draw simulation box edges if requested
   every drawn pixel has depth so merge can decide which to keep
------------------------------------------------------------------------- */

void DumpImage::create_image()
{
  int i,j,k,m,n,itype,atom1,atom2,imol,iatom,btype,ibonus,drawflag;
  tagint tagprev;
  double diameter,opacity,delx,dely,delz;
  int *bodyvec;
  double **bodyarray;
  double *color,*color1,*color2;
  double *p1,*p2,*p3;
  double pt1[3],pt2[3],pt3[3];
  double mat[3][3];

  // set defaults to avoid accidental uninitialized accesses
  diameter = 1.0;
  opacity = 1.0;

  // render my atoms

  if (atomflag) {
    double **x = atom->x;
    int *type = atom->type;
    int *line = atom->line;
    int *tri = atom->tri;
    int *ellipsoid = atom->ellipsoid;
    int *body = atom->body;

    m = 0;
    for (i = 0; i < nchoose; i++) {
      j = clist[i];

      itype = type[j];
      opacity = aopacity[itype];
      if (acolor == TYPE) {
        color = colortype[itype];
      } else if (acolor == ELEMENT) {
        color = colorelement[itype];
      } else if (acolor == ATTRIBUTE) {
        color = image->map_value2color(0,buf[m]);
      } else color = image->color2rgb("white");

      if (adiam == NUMERIC) {
        diameter = adiamvalue;
      } else if (adiam == TYPE) {
        diameter = diamtype[itype];
      } else if (adiam == ELEMENT) {
        diameter = diamelement[itype];
      } else if (adiam == ATTRIBUTE) {
        diameter = buf[m+1];
      }

      // do not draw if line,tri,body keywords enabled and atom is one of those

      drawflag = 1;
      if (extraflag) {
        if (lineflag && line[j] >= 0) drawflag = 0;
        if (triflag && tri[j] >= 0) drawflag = 0;
        if (ellipsoidflag && ellipsoid[j] >= 0) drawflag = 0;
        if (bodyflag && body[j] >= 0) drawflag = 0;
      }

      if (drawflag) image->draw_sphere(x[j],color,diameter,opacity);

      m += size_one;
    }
  }

  // render my grid cells
  // 2 triangles for 2d rectangle, 12 triangles for 3d cube surface
  // grid_cell_corners_2d/3d calculates orthogonal vs triclinic corner pts
  // for 3d, outward normals on all 6 faces

  if (gridflag) {

    // reset lighting for flat surfaces to make them brighter

    auto saved = reset_lighting(image, 0.9, 0.3, 0.3, 0.3);

    int n = 0;
    if (domain->dimension == 2) {
      for (int iy = nylo_in; iy <= nyhi_in; iy++)
        for (int ix = nxlo_in; ix <= nxhi_in; ix++) {
          grid_cell_corners_2d(ix,iy);
          color = image->map_value2color(1,gbuf[n++]);
          image->draw_triangle(gcorners[0],gcorners[1],gcorners[3],color);
          image->draw_triangle(gcorners[0],gcorners[3],gcorners[2],color);
        }
    } else {
      for (int iz = nzlo_in; iz <= nzhi_in; iz++)
        for (int iy = nylo_in; iy <= nyhi_in; iy++)
          for (int ix = nxlo_in; ix <= nxhi_in; ix++) {
            grid_cell_corners_3d(ix,iy,iz);
            color = image->map_value2color(1,gbuf[n++]);
            // lower x face
            image->draw_triangle(gcorners[0],gcorners[4],gcorners[6],color);
            image->draw_triangle(gcorners[0],gcorners[6],gcorners[2],color);
            // upper x face
            image->draw_triangle(gcorners[1],gcorners[5],gcorners[7],color);
            image->draw_triangle(gcorners[1],gcorners[7],gcorners[3],color);
            // lower y face
            image->draw_triangle(gcorners[0],gcorners[1],gcorners[5],color);
            image->draw_triangle(gcorners[0],gcorners[5],gcorners[4],color);
            // upper y face
            image->draw_triangle(gcorners[2],gcorners[6],gcorners[7],color);
            image->draw_triangle(gcorners[2],gcorners[7],gcorners[3],color);
            // lower z face
            image->draw_triangle(gcorners[0],gcorners[2],gcorners[3],color);
            image->draw_triangle(gcorners[0],gcorners[3],gcorners[1],color);
            // upper z face
            image->draw_triangle(gcorners[4],gcorners[5],gcorners[7],color);
            image->draw_triangle(gcorners[4],gcorners[7],gcorners[6],color);
          }
    }

    // restore lighting for curved objects

    restore_lighting(saved, image);
  }

  // render atoms that are lines

  if (lineflag) {
    double length,theta,dx,dy;
    double **x = atom->x;
    int *line = atom->line;
    int *type = atom->type;

    m = 0;
    for (i = 0; i < nchoose; i++) {
      j = clist[i];
      if (line[j] < 0) continue;

      itype = type[j];
      opacity = aopacity[itype];
      if (lcolor == TYPE) {
        color = colortype[itype];
      } else if (lcolor == INDEX) {
        itype = (line[j] % atom->ntypes) + 1;
        color = colortype[itype];
      } else if (lcolor == ATOM) {
        if (acolor == TYPE) {
          color = colortype[itype];
        } else if (acolor == ELEMENT) {
          color = colorelement[itype];
        } else if (acolor == ATTRIBUTE) {
          color = image->map_value2color(0,buf[m]);
        } else {
          color = image->color2rgb("white");
        }
      } else {
        color = image->color2rgb("white");
      }

      if (ldiam == NUMERIC) {
        diameter = ldiamvalue;
      }

      length = avec_line->bonus[line[j]].length;
      theta = avec_line->bonus[line[j]].theta;
      dx = 0.5*length*cos(theta);
      dy = 0.5*length*sin(theta);

      pt1[0] = x[j][0] + dx;
      pt1[1] = x[j][1] + dy;
      pt1[2] = 0.0;
      pt2[0] = x[j][0] - dx;
      pt2[1] = x[j][1] - dy;
      pt2[2] = 0.0;

      image->draw_cylinder(pt1,pt2,color,ldiamvalue,3,opacity);

      m += size_one;
    }
  }

  // render atoms that are triangles
  // tstyle = 1 for tri only, 2 for edges only, 3 for both

  if (triflag) {
    int tridraw = 1;
    if (tstyle == 2) tridraw = 0;
    int edgedraw = 1;
    if (tstyle == 1) edgedraw = 0;

    double **x = atom->x;
    int *tri = atom->tri;
    int *type = atom->type;

    m = 0;
    for (i = 0; i < nchoose; i++) {
      j = clist[i];
      if (tri[j] < 0) continue;

      itype = type[j];
      opacity = aopacity[itype];
      if (tcolor == TYPE) {
        color = colortype[itype];
      } else if (tcolor == INDEX) {
        itype = (tri[j] % atom->ntypes) + 1;
        color = colortype[itype];
      } else if (tcolor == ATOM) {
        if (acolor == TYPE) {
          color = colortype[itype];
        } else if (acolor == ELEMENT) {
          color = colorelement[itype];
        } else if (acolor == ATTRIBUTE) {
          color = image->map_value2color(0,buf[m]);
        } else {
          color = image->color2rgb("white");
        }
      } else {
        color = image->color2rgb("white");
      }

      MathExtra::quat_to_mat(avec_tri->bonus[tri[j]].quat,mat);
      MathExtra::matvec(mat,avec_tri->bonus[tri[j]].c1,pt1);
      MathExtra::matvec(mat,avec_tri->bonus[tri[j]].c2,pt2);
      MathExtra::matvec(mat,avec_tri->bonus[tri[j]].c3,pt3);
      MathExtra::add3(pt1,x[j],pt1);
      MathExtra::add3(pt2,x[j],pt2);
      MathExtra::add3(pt3,x[j],pt3);

      if (tridraw) {
        // brighten flat surfaces a bit
        auto saved = reset_lighting(image, 0.3, 0.8, 0.45, 0.8);

        image->draw_triangle(pt1,pt2,pt3,color,opacity);

        // restore previous settings
        restore_lighting(saved, image);
      }
      if (edgedraw) {
        image->draw_cylinder(pt1,pt2,color,tdiamvalue,3,opacity);
        image->draw_cylinder(pt2,pt3,color,tdiamvalue,3,opacity);
        image->draw_cylinder(pt3,pt1,color,tdiamvalue,3,opacity);
      }
      m += size_one;
    }
  }

  // render atoms that are ellipsoids
  // estyle = 1 for tri only, estyle 2 for wireframe only

  if (ellipsoidflag) {
    double **x = atom->x;
    int *ellipsoid = atom->ellipsoid;
    int *type = atom->type;
    m = 0;

    for (i = 0; i < nchoose; i++) {
      j = clist[i];
      if (ellipsoid[j] < 0) continue;

      itype = type[j];
      opacity = aopacity[itype];
      if (ecolor == TYPE) {
        color = colortype[itype];
      } else if (ecolor == INDEX) {
        itype = (ellipsoid[j] % atom->ntypes) + 1;
        color = colortype[itype];
      } else if (ecolor == ATOM) {
        if (acolor == TYPE) {
          color = colortype[itype];
        } else if (acolor == ELEMENT) {
          color = colorelement[itype];
        } else if (acolor == ATTRIBUTE) {
          color = image->map_value2color(0,buf[m]);
        } else {
          color = image->color2rgb("white");
        }
      } else {
        color = image->color2rgb("white");
      }
      EllipsoidObj e(elevel);
      if (avec_ellipsoid->bonus_super) {
        auto *bonus = avec_ellipsoid->bonus_super;
        e.draw(image, estyle, color, x[j], bonus[ellipsoid[j]].shape, bonus[ellipsoid[j]].quat,
               ediamvalue, opacity, bonus[ellipsoid[j]].block);
      } else {
        auto *bonus = avec_ellipsoid->bonus;
        e.draw(image, estyle, color, x[j], bonus[ellipsoid[j]].shape,bonus[ellipsoid[j]].quat,
               ediamvalue, opacity, nullptr);
      }
      m += size_one;
    }
  }

  // render atoms that are bodies

  if (bodyflag) {
    Body *bptr = avec_body->bptr;
    int *body = atom->body;
    int *type = atom->type;

    m = 0;
    for (i = 0; i < nchoose; i++) {
      j = clist[i];
      if (body[j] < 0) continue;

      itype = type[j];
      opacity = aopacity[itype];
      if (bodycolor == TYPE) {
        color = colortype[itype];
      } else if (bodycolor == INDEX) {
        itype = (body[j] % atom->ntypes) + 1;
        color = colortype[itype];
      } else if (bodycolor == ATOM) {
        if (acolor == TYPE) {
          color = colortype[itype];
        } else if (acolor == ELEMENT) {
          color = colorelement[itype];
        } else if (acolor == ATTRIBUTE) {
          color = image->map_value2color(0,buf[m]);
        } else {
          color = image->color2rgb("white");
        }
      } else {
        color = image->color2rgb("white");
      }

      ibonus = body[j];
      n = bptr->image(ibonus,bodyflag1,bodyflag2,bodyvec,bodyarray);
      for (k = 0; k < n; k++) {
        if (bodyvec[k] == Graphics::SPHERE)
          image->draw_sphere(bodyarray[k],color,bodyarray[k][3],opacity);
        else if (bodyvec[k] == Graphics::LINE)
          image->draw_cylinder(&bodyarray[k][0],&bodyarray[k][3],color,bodyarray[k][6],3,opacity);
        else if (bodyvec[k] == Graphics::TRI) {
          // brighten flat surfaces somewhat
          auto saved = reset_lighting(image, 0.3, 0.8, 0.45, 0.8);
          image->draw_triangle(&bodyarray[k][0],&bodyarray[k][3],&bodyarray[k][6],color,opacity);
          // restore previous settings
          restore_lighting(saved, image);
        } else if (bodyvec[k] == Graphics::TRINORM) {
          // brighten surfaces a little bit
          auto saved = reset_lighting(image, 0.6, 0.3, 0.5, 0.7);
          image->draw_trinorm(&bodyarray[k][0],&bodyarray[k][3],&bodyarray[k][6],
                              &bodyarray[k][9],&bodyarray[k][12],&bodyarray[k][15],
                              color,color,color,opacity);
          // restore previous settings
          restore_lighting(saved, image);
        }
      }

      m += size_one;
    }
  }

  // render bonds for my atoms
  // both atoms in bond must be selected for bond to be rendered
  // if newton_bond is off, only render bond once
  // render bond in 2 pieces if crosses periodic boundary
  // if bond is deleted (type = 0), do not render
  // if bond is turned off (type < 0), still render

  if (bondflag == YES) {
    double **x = atom->x;
    tagint *tag = atom->tag;
    tagint **bond_atom = atom->bond_atom;
    int **bond_type = atom->bond_type;
    int *num_bond = atom->num_bond;
    int *molindex = atom->molindex;
    int *molatom = atom->molatom;
    int *type = atom->type;
    int nlocal = atom->nlocal;
    int newton_bond = force->newton_bond;
    int molecular = atom->molecular;
    Molecule **onemols = atom->avec->onemols;

    // communicate choose flag for ghost atoms to know if they are selected
    // if bcolor/bdiam = ATOM, setup bufcopy to comm atom color/diam attributes

    if (atom->nmax > maxbufcopy) {
      maxbufcopy = atom->nmax;
      memory->destroy(chooseghost);
      memory->create(chooseghost,maxbufcopy,"dump:chooseghost");
      if (comm_forward == 3) {
        memory->destroy(bufcopy);
        memory->create(bufcopy,maxbufcopy,2,"dump:bufcopy");
      }
    }

    for (i = 0; i < nlocal; i++) chooseghost[i] = choose[i];

    if (comm_forward == 3) {
      for (i = 0; i < nlocal; i++) bufcopy[i][0] = bufcopy[i][1] = 0.0;
      m = 0;
      for (i = 0; i < nchoose; i++) {
        j = clist[i];
        bufcopy[j][0] = buf[m];
        bufcopy[j][1] = buf[m+1];
        m += size_one;
      }
    }

    comm->forward_comm(this);

    for (i = 0; i < nchoose; i++) {
      atom1 = clist[i];
      if (molecular == Atom::MOLECULAR) n = num_bond[atom1];
      else {
        if (molindex[atom1] < 0) continue;
        imol = molindex[atom1];
        iatom = molatom[atom1];
        n = onemols[imol]->num_bond[iatom];
      }

      for (m = 0; m < n; m++) {
        if (molecular == Atom::MOLECULAR) {
          btype = bond_type[atom1][m];
          atom2 = atom->map(bond_atom[atom1][m]);
        } else {
          tagprev = tag[i] - iatom - 1;
          btype = atom->map(onemols[imol]->bond_type[iatom][m]);
          atom2 = atom->map(onemols[imol]->bond_atom[iatom][m]+tagprev);
        }

        if (atom2 < 0 || !chooseghost[atom2]) continue;
        if (newton_bond == 0 && tag[atom1] > tag[atom2]) continue;
        if (btype == 0) continue;
        if (btype < 0) btype = -btype;

        if (bcolor == ATOM) {
          if (acolor == TYPE) {
            color1 = colortype[type[atom1]];
            color2 = colortype[type[atom2]];
          } else if (acolor == ELEMENT) {
            color1 = colorelement[type[atom1]];
            color2 = colorelement[type[atom2]];
          } else if (acolor == ATTRIBUTE) {
            color1 = image->map_value2color(0,bufcopy[atom1][0]);
            color2 = image->map_value2color(0,bufcopy[atom2][0]);
          } else {
            color1 = image->color2rgb("white");
            color2 = image->color2rgb("white");
          }
        } else if (bcolor == TYPE) {
          color = bcolortype[btype];
        }

        if (bdiam == NUMERIC) {
          diameter = bdiamvalue;
        } else if (bdiam == ATOM) {
          if (adiam == NUMERIC) {
            diameter = adiamvalue;
          } else if (adiam == TYPE) {
            diameter = MIN(diamtype[type[atom1]],diamtype[type[atom2]]);
          } else if (adiam == ELEMENT) {
            diameter = MIN(diamelement[type[atom1]],diamelement[type[atom2]]);
          } else if (adiam == ATTRIBUTE) {
            diameter = MIN(bufcopy[atom1][1],bufcopy[atom2][1]);
          }
        } else if (bdiam == TYPE) {
          diameter = bdiamtype[btype];
        }

        // draw cylinder in 2 pieces if bcolor = ATOM
        // or bond crosses periodic boundary

        delx = x[atom2][0] - x[atom1][0];
        dely = x[atom2][1] - x[atom1][1];
        delz = x[atom2][2] - x[atom1][2];

        if (bcolor == ATOM || domain->minimum_image_check(delx,dely,delz)) {
          domain->minimum_image(FLERR, delx,dely,delz);
          double xmid[3];
          xmid[0] = x[atom1][0] + 0.5*delx;
          xmid[1] = x[atom1][1] + 0.5*dely;
          xmid[2] = x[atom1][2] + 0.5*delz;
          if (bcolor == ATOM)
            image->draw_cylinder(x[atom1],xmid,color1,diameter,3,aopacity[type[atom1]]);
          else image->draw_cylinder(x[atom1],xmid,color,diameter,3,bopacity[btype]);
          xmid[0] = x[atom2][0] - 0.5*delx;
          xmid[1] = x[atom2][1] - 0.5*dely;
          xmid[2] = x[atom2][2] - 0.5*delz;
          if (bcolor == ATOM)
            image->draw_cylinder(xmid,x[atom2],color2,diameter,3,aopacity[type[atom1]]);
          else image->draw_cylinder(xmid,x[atom2],color,diameter,3,bopacity[btype]);

        } else image->draw_cylinder(x[atom1],x[atom2],color,diameter,3,bopacity[btype]);
      }
    }
  }

  // render dynamic bonds for my atoms
  // both atoms in pair must be selected and closer than bond cutoff for bond to be rendered
  // if newton_bond is off, only render bond once
  // render bond in 2 pieces if crosses periodic boundary

  if (bondflag == AUTO) {
    // grab some suitable neighbor list, if available
    auto *list = neighbor->get_best_pair_list();
    if (!list) {
      if (comm->me == 0)
        error->warning(FLERR, "No suitable existing neighbor list for dump image autobond found");
    } else {
      int nlocal = atom->nlocal;

      // communicate choose flag for ghost atoms to know if they are selected
      // if bcolor/bdiam = ATOM, setup bufcopy to comm atom color/diam attributes

      if (atom->nmax > maxbufcopy) {
        maxbufcopy = atom->nmax;
        memory->destroy(chooseghost);
        memory->create(chooseghost,maxbufcopy,"dump:chooseghost");
        if (comm_forward == 3) {
          memory->destroy(bufcopy);
          memory->create(bufcopy,maxbufcopy,2,"dump:bufcopy");
        }
      }
      for (i = 0; i < nlocal; i++) chooseghost[i] = choose[i];

      if (comm_forward == 3) {
        for (i = 0; i < nlocal; i++) bufcopy[i][0] = bufcopy[i][1] = 0.0;
        m = 0;
        for (i = 0; i < nchoose; i++) {
          j = clist[i];
          bufcopy[j][0] = buf[m];
          bufcopy[j][1] = buf[m+1];
          m += size_one;
        }
      }

      comm->forward_comm(this);

      double **x = atom->x;
      double *mass = atom->mass;
      double *rmass = atom->rmass;
      tagint *tag = atom->tag;
      int *type = atom->type;

      int *numneigh = list->numneigh;
      int **firstneigh = list->firstneigh;
      const double cutsq = bondcutoff * bondcutoff;
      int newton_pair = force->newton_pair;
      bool checkmass = false;

      // check for element by mass only for units "real" or "metal"
      if ((strcmp(update->unit_style, "real") == 0) || (strcmp(update->unit_style, "metal") == 0))
        checkmass = true;

      for (int ii = 0; ii < nchoose; ii++) {
        atom1 = clist[ii];
        const double xtmp = x[atom1][0];
        const double ytmp = x[atom1][1];
        const double ztmp = x[atom1][2];

        // loop over neighbors
        auto *jlist = firstneigh[atom1];
        const int jnum = numneigh[atom1];
        for (int jj = 0; jj < jnum; ++jj) {
          atom2 = jlist[jj] & NEIGHMASK;
          if (!chooseghost[atom2]) continue;
          if ((newton_pair == 0) && (tag[atom1] > tag[atom2])) continue;
          // skip hydrogen-hydrogen bonds for units real or metal based on their mass
          // this is primarily for water and alkyl groups
          if (checkmass) {
            if (rmass) {
              if ((rmass[atom1] < 3.0) && (rmass[atom2] < 3.0)) continue;
            } else {
              if ((mass[type[atom1]] < 3.0) && (mass[type[atom2]] < 3.0)) continue;
            }
          }
          double dx = x[atom2][0] - xtmp;
          double dy = x[atom2][1] - ytmp;
          double dz = x[atom2][2] - ztmp;
          const double rsq = dx*dx + dy*dy + dz*dz;

          if (rsq < cutsq) {
            if (acolor == TYPE) {
              color1 = colortype[type[atom1]];
              color2 = colortype[type[atom2]];
            } else if (acolor == ELEMENT) {
              color1 = colorelement[type[atom1]];
              color2 = colorelement[type[atom2]];
            } else if (acolor == ATTRIBUTE) {
              color1 = image->map_value2color(0,bufcopy[atom1][0]);
              color2 = image->map_value2color(0,bufcopy[atom2][0]);
            } else {
              color1 = image->color2rgb("white");
              color2 = image->color2rgb("white");
            }
            diameter = bdiamvalue;

            // draw cylinder in 2 pieces if bcolor = ATOM
            // or bond crosses periodic boundary
            double xmid[3];

            domain->minimum_image(FLERR,dx,dy,dz);
            xmid[0] = x[atom1][0] + 0.5*dx;
            xmid[1] = x[atom1][1] + 0.5*dy;
            xmid[2] = x[atom1][2] + 0.5*dz;
            image->draw_cylinder(x[atom1],xmid,color1,diameter,3,aopacity[type[atom1]]);

            xmid[0] = x[atom2][0] - 0.5*dx;
            xmid[1] = x[atom2][1] - 0.5*dy;
            xmid[2] = x[atom2][2] - 0.5*dz;
            image->draw_cylinder(xmid,x[atom2],color2,diameter,3,aopacity[type[atom2]]);
          }
        }
      }
    }
  }

  // clang-format on
  // render objects provided by computes and fixes

  for (const auto &iobj : objects) {
    int *objvec = nullptr;
    double **objarray = nullptr;
    const int ntypes = atom->ntypes;
    if (iobj.cptr) {
      n = iobj.cptr->compute_image(objvec, objarray);
    } else if (iobj.fptr) {
      n = iobj.fptr->image(objvec, objarray);
    } else {
      continue;    // ignore objects without pointer
    }

    for (i = 0; i < n; i++) {
      if (!objvec || !objarray) continue;

      // set color and transparency
      double *white = image->color2rgb("white");
      opacity = iobj.opacity;
      if (iobj.colorstyle == TYPE) {
        itype = static_cast<int>(objarray[i][0] - 1.0) % ntypes + 1;
        color = (itype > 0) ? colortype[itype] : white;
      } else if (iobj.colorstyle == ELEMENT) {
        itype = static_cast<int>(objarray[i][0] - 1.0) % ntypes + 1;
        color = (itype > 0) ? colorelement[itype] : white;
      } else if (iobj.colorstyle == CONSTANT) {
        color = iobj.rgb;
      } else {
        color = white;
        opacity = 1.0;
      }

      if (objvec[i] == Graphics::SPHERE) {
        diameter = objarray[i][4];
        if (objarray[i][4] < 0) {
          if (adiam == NUMERIC) {
            diameter = adiamvalue;
          } else if (adiam == TYPE) {
            diameter = diamtype[itype];
          } else if (adiam == ELEMENT) {
            diameter = diamelement[itype];
          }
        }
        image->draw_sphere(&objarray[i][1], color, diameter + iobj.flag2, opacity);
      } else if (objvec[i] == Graphics::LINE) {
        // @sjplimp for consistency this should be:
        // image->draw_cylinder(&objarray[i][1],&objarray[i][4],color,iobj.flag2,iobj.flag1);
        image->draw_cylinder(&objarray[i][1], &objarray[i][4], color, iobj.flag1, 3, opacity);

      } else if (objvec[i] == Graphics::TRI) {    // don't render surface meshes in 2d
        if (domain->dimension == 3) {
          p1 = &objarray[i][1];
          p2 = &objarray[i][4];
          p3 = &objarray[i][7];
          if (static_cast<int>(iobj.flag1) % 2) {
            image->draw_triangle(p1, p2, p3, color, opacity);
          } else {
            image->draw_cylinder(p1, p2, color, iobj.flag2, 3, opacity);
            image->draw_cylinder(p2, p3, color, iobj.flag2, 3, opacity);
            image->draw_cylinder(p3, p1, color, iobj.flag2, 3, opacity);
          }
        }

      } else if (objvec[i] == Graphics::TRINORM) {    // don't render surface meshes in 2d
        if (domain->dimension == 3) {
          double *color2 = nullptr;
          double *color3 = nullptr;
          if (iobj.colorstyle == TYPE) {
            itype = static_cast<int>(objarray[i][1] - 1.0) % ntypes + 1;
            color2 = (itype > 0) ? colortype[itype] : white;
            itype = static_cast<int>(objarray[i][2] - 1.0) % ntypes + 1;
            color3 = (itype > 0) ? colortype[itype] : white;
          } else if (iobj.colorstyle == ELEMENT) {
            itype = static_cast<int>(objarray[i][1] - 1.0) % ntypes + 1;
            color2 = (itype > 0) ? colorelement[itype] : white;
            itype = static_cast<int>(objarray[i][2] - 1.0) % ntypes + 1;
            color3 = (itype > 0) ? colorelement[itype] : white;
          } else if (iobj.colorstyle == CONSTANT) {
            color2 = iobj.rgb;
            color3 = iobj.rgb;
          } else {
            color2 = color3 = image->color2rgb("white");
            opacity = 1.0;
          }

          p1 = &objarray[i][3];
          p2 = &objarray[i][6];
          p3 = &objarray[i][9];
          if (static_cast<int>(iobj.flag1) % 2) {
            image->draw_trinorm(p1, p2, p3, &objarray[i][12], &objarray[i][15], &objarray[i][18],
                                color, color2, color3, opacity);
          } else {
            image->draw_cylinder(p1, p2, color, iobj.flag2, 3, opacity);
            image->draw_cylinder(p2, p3, color, iobj.flag2, 3, opacity);
            image->draw_cylinder(p3, p1, color, iobj.flag2, 3, opacity);
          }
        }

      } else if (objvec[i] == Graphics::CYLINDER) {
        image->draw_cylinder(&objarray[i][1], &objarray[i][4], color, objarray[i][7] + iobj.flag2,
                             (int) iobj.flag1, opacity);

      } else if (objvec[i] == Graphics::TRIANGLE) {
        image->draw_triangle(&objarray[i][1], &objarray[i][4], &objarray[i][7], color, opacity);

      } else if (objvec[i] == Graphics::ARROW) {
        ArrowObj a(objarray[i][9]);
        a.draw(image, color, &objarray[i][1], objarray[i][7] + iobj.flag1, &objarray[i][4],
               objarray[i][8] + iobj.flag2, opacity);

      } else if (objvec[i] == Graphics::CONE) {
        ConeObj c(1.0, objarray[i][7] + iobj.flag2, objarray[i][8] + iobj.flag2,
                  (int) objarray[i][9]);
        c.draw(image, vec3{objarray[i][1], objarray[i][2], objarray[i][3]},
               vec3{objarray[i][4], objarray[i][5], objarray[i][6]}, color, opacity);

      } else if (objvec[i] == Graphics::PIXMAP) {
        // get pointer to pixmap buffer and get background transparency color
        const auto *pixmap = (const unsigned char *) ubuf(objarray[i][6]).i;    // NOLINT
        double transcolor[3] = {objarray[i][7], objarray[i][8], objarray[i][9]};
        if (iobj.flag1 == 0.0)    // coordinates are in box coordinates
          image->draw_pixmap(&objarray[i][1], (int) objarray[i][4], (int) objarray[i][5], pixmap,
                             transcolor, objarray[i][10], opacity);
        else    // coordinates are in image coordinates, ignore z
          image->draw_pixmap((int) objarray[i][1], (int) objarray[i][2], (int) objarray[i][4],
                             (int) objarray[i][5], pixmap, transcolor, objarray[i][10], opacity);

      } else if (objvec[i] == Graphics::BOND) {
        int type1 = static_cast<int>(objarray[i][0] - 1.0) % ntypes + 1;
        int type2 = static_cast<int>(objarray[i][1] - 1.0) % ntypes + 1;
        double *color1, *color2;
        opacity = iobj.opacity;
        if (iobj.colorstyle == TYPE) {
          color1 = colortype[type1];
          color2 = colortype[type2];
        } else if (iobj.colorstyle == ELEMENT) {
          color1 = colorelement[type1];
          color2 = colorelement[type2];
        } else if (iobj.colorstyle == CONSTANT) {
          color1 = iobj.rgb;
          color2 = iobj.rgb;
        } else {
          color1 = image->color2rgb("white");
          color2 = image->color2rgb("white");
        }

        diameter = 0.5;
        if (bdiam == ATOM) {
          if (adiam == NUMERIC) {
            diameter = adiamvalue;
          } else if (adiam == TYPE) {
            diameter = MIN(diamtype[type1], diamtype[type2]);
          } else if (adiam == ELEMENT) {
            diameter = MIN(diamelement[type1], diamelement[type2]);
          } else if (adiam == ATTRIBUTE) {
            diameter = MIN(bufcopy[atom1][1], bufcopy[atom2][1]);
          }
        } else {
          diameter = bdiamvalue;
        }
        // bond diameter adjustment from dump image command line
        diameter += iobj.flag2;

        // draw bond cylinder in 2 pieces

        int capflag = (iobj.flag1 != 0.0) ? 3 : 0;
        delx = objarray[i][5] - objarray[i][2];
        dely = objarray[i][6] - objarray[i][3];
        delz = objarray[i][7] - objarray[i][4];

        domain->minimum_image(FLERR, delx, dely, delz);
        double xmid[3];
        xmid[0] = objarray[i][2] + 0.5 * delx;
        xmid[1] = objarray[i][3] + 0.5 * dely;
        xmid[2] = objarray[i][4] + 0.5 * delz;
        image->draw_cylinder(&objarray[i][2], xmid, color1, diameter, capflag, opacity);
        xmid[0] = objarray[i][5] - 0.5 * delx;
        xmid[1] = objarray[i][6] - 0.5 * dely;
        xmid[2] = objarray[i][7] - 0.5 * delz;
        image->draw_cylinder(xmid, &objarray[i][5], color2, diameter, capflag, opacity);
      }
    }
  }

  // render regions, if supported

  for (auto &reg : regions) {

    // check if region has changed or went away
    auto *ptr = domain->get_region_by_id(reg.id);
    if (!ptr) error->all(FLERR, "Dump image region {} does not exist", reg.id);
    reg.ptr = ptr;

    // update internal variables

    reg.ptr->prematch();

    // for POINTS style we have the same code for all region styles

    if (reg.style == POINTS) {
      int seed = (int) (platform::walltime() * 1000000) % 1000000;
      RanMars rand(lmp, seed);
      double pos[3];

      double xoff = domain->boxlo[0];
      double yoff = domain->boxlo[1];
      double zoff = domain->boxlo[2];
      double xlen = domain->xprd;
      double ylen = domain->yprd;
      double zlen = domain->zprd;
      if (domain->triclinic) {
        xoff = domain->boxlo_bound[0];
        yoff = domain->boxlo_bound[1];
        zoff = domain->boxlo_bound[2];
        xlen = domain->boxhi_bound[0] - domain->boxlo_bound[0];
        ylen = domain->boxhi_bound[1] - domain->boxlo_bound[1];
        zlen = domain->boxhi_bound[2] - domain->boxlo_bound[2];
      }

      // temporarily turn off the region's open faces flag since
      // otherwise Region::match() will match for *any* position
      int saved_flag = reg.ptr->openflag;
      reg.ptr->openflag = 0;

      for (int i = 0; i < reg.npoints; ++i) {
        pos[0] = rand.uniform() * xlen + xoff;
        pos[1] = rand.uniform() * ylen + yoff;
        pos[2] = rand.uniform() * zlen + zoff;
        if (reg.ptr->match(pos[0], pos[1], pos[2]))
          image->draw_sphere(pos, reg.color, reg.diameter);
      }

      // restore saved open face flag
      reg.ptr->openflag = saved_flag;

    } else {
      std::string regstyle = reg.ptr->style;

      // prism differs from block only in the positions of the corners; handle them together
      if ((regstyle == "block") || (regstyle == "prism")) {

        // set the positions of the corners
        using cornerdata = std::array<vec3, 8>;
        cornerdata corners;

        if (regstyle == "block") {
          auto *myreg = dynamic_cast<RegBlock *>(reg.ptr);
          // inconsistent style. should not happen.
          if (!myreg) continue;

          // clamp region boundaries to box boundaries
          double xlo = MAX(myreg->xlo, domain->boxlo[0]);
          double ylo = MAX(myreg->ylo, domain->boxlo[1]);
          double zlo = MAX(myreg->zlo, domain->boxlo[2]);
          double xhi = MIN(myreg->xhi, domain->boxhi[0]);
          double yhi = MIN(myreg->yhi, domain->boxhi[1]);
          double zhi = MIN(myreg->zhi, domain->boxhi[2]);

          corners = cornerdata{vec3{xlo, ylo, zlo}, vec3{xlo, ylo, zhi}, vec3{xlo, yhi, zhi},
                               vec3{xlo, yhi, zlo}, vec3{xhi, ylo, zlo}, vec3{xhi, ylo, zhi},
                               vec3{xhi, yhi, zhi}, vec3{xhi, yhi, zlo}};
        }

        if (regstyle == "prism") {
          auto *myreg = dynamic_cast<RegPrism *>(reg.ptr);
          // inconsistent style. should not happen.
          if (!myreg) continue;

          // clamp region boundaries to box boundaries
          double xlo = MAX(myreg->xlo, domain->boxlo[0]);
          double ylo = MAX(myreg->ylo, domain->boxlo[1]);
          double zlo = MAX(myreg->zlo, domain->boxlo[2]);
          double xhi = MIN(myreg->xhi, domain->boxhi[0]);
          double yhi = MIN(myreg->yhi, domain->boxhi[1]);
          double zhi = MIN(myreg->zhi, domain->boxhi[2]);

          corners = cornerdata{vec3{xlo, ylo, zlo},
                               vec3{xlo + myreg->xz, ylo + myreg->yz, zhi},
                               vec3{xlo + myreg->xy + myreg->xz, yhi + myreg->yz, zhi},
                               vec3{xlo + myreg->xy, yhi, zlo},
                               vec3{xhi, ylo, zlo},
                               vec3{xhi + myreg->xz, ylo + myreg->yz, zhi},
                               vec3{xhi + myreg->xy + myreg->xz, yhi + myreg->yz, zhi},
                               vec3{xhi + myreg->xy, yhi, zlo}};
        }

        for (int i = 0; i < 8; ++i)
          reg.ptr->forward_transform(corners[i][0], corners[i][1], corners[i][2]);

#define DRAW_CYLINDER(j, k) \
  image->draw_cylinder(corners[j].data(), corners[k].data(), reg.color, reg.diameter, 3, 1.0)
#define DRAW_TRIANGLE(j, k, l) \
  image->draw_triangle(corners[j].data(), corners[k].data(), corners[l].data(), reg.color, opacity)

        if (reg.style == FRAME) {
          DRAW_CYLINDER(0, 1);
          DRAW_CYLINDER(1, 2);
          DRAW_CYLINDER(0, 3);
          DRAW_CYLINDER(2, 3);
          DRAW_CYLINDER(0, 4);
          DRAW_CYLINDER(1, 5);
          DRAW_CYLINDER(2, 6);
          DRAW_CYLINDER(3, 7);
          DRAW_CYLINDER(4, 5);
          DRAW_CYLINDER(5, 6);
          DRAW_CYLINDER(4, 7);
          DRAW_CYLINDER(6, 7);
        } else if ((reg.style == FILLED) || (reg.style == TRANSPARENT)) {
          opacity = (reg.style == TRANSPARENT) ? reg.opacity : 1.0;
          if (!reg.ptr->open_faces[0]) {
            DRAW_TRIANGLE(0, 1, 2);
            DRAW_TRIANGLE(2, 3, 0);
          }
          if (!reg.ptr->open_faces[1]) {
            DRAW_TRIANGLE(4, 5, 6);
            DRAW_TRIANGLE(6, 7, 4);
          }
          if (!reg.ptr->open_faces[2]) {
            DRAW_TRIANGLE(0, 4, 7);
            DRAW_TRIANGLE(7, 3, 0);
          }
          if (!reg.ptr->open_faces[3]) {
            DRAW_TRIANGLE(1, 2, 6);
            DRAW_TRIANGLE(6, 5, 1);
          }
          if (!reg.ptr->open_faces[4]) {
            DRAW_TRIANGLE(0, 1, 5);
            DRAW_TRIANGLE(5, 4, 0);
          }
          if (!reg.ptr->open_faces[5]) {
            DRAW_TRIANGLE(3, 2, 6);
            DRAW_TRIANGLE(6, 7, 3);
          }
        }
#undef DRAW_CYLINDER
#undef DRAW_TRIANGLE

        // cylinder is just a special case of cone so we can handle them together
      } else if ((regstyle == "cone") || (regstyle == "cylinder")) {

        // construct triangle mesh and store direction and hi/lo center coordinates
        vec3 lo, hi;
        double length = 0.0, radiuslo = 0.0, radiushi = 0.0;
        double xdir = 0.0, ydir = 0.0, zdir = 0.0;

        if (regstyle == "cone") {
          auto *myreg = dynamic_cast<RegCone *>(reg.ptr);
          // inconsistent style. should not happen.
          if (!myreg) continue;

          radiuslo = myreg->radiuslo;
          radiushi = myreg->radiushi;
          if (myreg->axis == 'x') {
            xdir = 1.0;
            double xlo = MAX(myreg->lo, domain->boxlo[0]);
            double xhi = MIN(myreg->hi, domain->boxhi[0]);
            length = xhi - xlo;
            lo = {xlo, myreg->c1, myreg->c2};
            hi = {xhi, myreg->c1, myreg->c2};
          } else if (myreg->axis == 'y') {
            ydir = 1.0;
            double ylo = MAX(myreg->lo, domain->boxlo[1]);
            double yhi = MIN(myreg->hi, domain->boxhi[1]);
            length = yhi - ylo;
            lo = {myreg->c1, ylo, myreg->c2};
            hi = {myreg->c1, yhi, myreg->c2};
          } else {    // myreg->axis == 'z'
            zdir = 1.0;
            double zlo = MAX(myreg->lo, domain->boxlo[2]);
            double zhi = MIN(myreg->hi, domain->boxhi[2]);
            length = zhi - zlo;
            lo = {myreg->c1, myreg->c2, zlo};
            hi = {myreg->c1, myreg->c2, zhi};
          }
        }

        double radius = 0.0;
        if (regstyle == "cylinder") {
          auto *myreg = dynamic_cast<RegCylinder *>(reg.ptr);
          // inconsistent style. should not happen.
          if (!myreg) continue;

          length = myreg->hi - myreg->lo;
          radiuslo = myreg->radius;
          radiushi = myreg->radius;
          if (myreg->axis == 'x') {
            xdir = 1.0;
            double xlo = MAX(myreg->lo, domain->boxlo[0]);
            double xhi = MIN(myreg->hi, domain->boxhi[0]);
            length = xhi - xlo;
            lo = {xlo, myreg->c1, myreg->c2};
            hi = {xhi, myreg->c1, myreg->c2};
          } else if (myreg->axis == 'y') {
            ydir = 1.0;
            double ylo = MAX(myreg->lo, domain->boxlo[1]);
            double yhi = MIN(myreg->hi, domain->boxhi[1]);
            length = yhi - ylo;
            lo = {myreg->c1, ylo, myreg->c2};
            hi = {myreg->c1, yhi, myreg->c2};
          } else {    // myreg->axis == 'z'
            zdir = 1.0;
            double zlo = MAX(myreg->lo, domain->boxlo[2]);
            double zhi = MIN(myreg->hi, domain->boxhi[2]);
            length = zhi - zlo;
            lo = {myreg->c1, myreg->c2, zlo};
            hi = {myreg->c1, myreg->c2, zhi};
          }
        }

        // compute direction vector and positions of lo/hi centers

        vec3 dir{xdir, ydir, zdir};
        vec3 mid = 0.5 * (hi + lo);

        // determine which faces to draw
        int faceflag = 0;
        if (reg.style == FRAME) {
          faceflag = 4;    // no caps
        } else {
          if (!reg.ptr->open_faces[0]) faceflag |= 1;
          if (!reg.ptr->open_faces[1]) faceflag |= 2;
          if (!reg.ptr->open_faces[2]) faceflag |= 4;
        }

        // determine draw style flags
        int drawflag = (reg.style == FRAME) ? 2 : 1;
        opacity = (reg.style == TRANSPARENT) ? reg.opacity : 1.0;

        ConeObj c(length, radiushi, radiuslo, faceflag,
                  (reg.style == FRAME) ? RESOLUTION : 2 * RESOLUTION);
        c.draw(image, drawflag, dir, mid, reg.color, reg.ptr, reg.diameter, opacity);

        // draw an additional uncapped cylinder for a smoother image for filled cylinders only
        if ((regstyle == "cylinder") && ((reg.style == FILLED) || (reg.style == TRANSPARENT))) {

          reg.ptr->forward_transform(lo[0], lo[1], lo[2]);
          reg.ptr->forward_transform(hi[0], hi[1], hi[2]);
          if (!reg.ptr->open_faces[2])
            image->draw_cylinder(lo.data(), hi.data(), reg.color, 2.0 * radius, 0, opacity);
        }

      } else if (regstyle == "ellipsoid") {
        auto *myreg = dynamic_cast<RegEllipsoid *>(reg.ptr);
        // inconsistent style. should not happen.
        if (!myreg) continue;

        double center[3] = {myreg->xc, myreg->yc, myreg->zc};
        double radius[3] = {myreg->a, myreg->b, myreg->c};
        if (reg.style == FRAME) {
          // we use a moderate mesh level for wireframe so we can better see what is inside
          EllipsoidObj(4).draw(image, 2, reg.color, center, radius, reg.ptr, reg.diameter, 1.0);
        } else if ((reg.style == FILLED) || (reg.style == TRANSPARENT)) {
          // we get a smoother surface with a one step higher mesh level (and 4x more triangles)
          EllipsoidObj(5).draw(image, 1, reg.color, center, radius, reg.ptr, reg.diameter,
                               (reg.style == TRANSPARENT) ? reg.opacity : 1.0);
        }

      } else if (regstyle == "plane") {
        auto *myreg = dynamic_cast<RegPlane *>(reg.ptr);
        // inconsistent style. should not happen.
        if (!myreg) continue;

        // determine draw style flags
        int flag = 1;
        opacity = 1.0;
        if (reg.style == FRAME) {
          flag = 2;
        } else if (reg.style == TRANSPARENT) {
          opacity = reg.opacity;
        }
        double center[3] = {myreg->xp, myreg->yp, myreg->zp};
        double scale = MathExtra::len3(domain->prd);
        if (domain->triclinic)
          PlaneObj(6).draw(image, flag, reg.color, center, myreg->normal, domain->boxlo_bound,
                           domain->boxhi_bound, scale, myreg, reg.diameter, opacity);
        else
          PlaneObj(6).draw(image, flag, reg.color, center, myreg->normal, domain->boxlo,
                           domain->boxhi, scale, myreg, reg.diameter, opacity);

      } else if (regstyle == "sphere") {
        auto *myreg = dynamic_cast<RegSphere *>(reg.ptr);
        // inconsistent style. should not happen.
        if (!myreg) continue;

        double center[3] = {myreg->xc, myreg->yc, myreg->zc};
        if (reg.style == FRAME) {
          // use wireframe mode of ellipsoid with three identical radii at a medium mesh level
          double radius[3] = {myreg->radius, myreg->radius, myreg->radius};
          EllipsoidObj(4).draw(image, 2, reg.color, center, radius, reg.ptr, reg.diameter, 1.0);
        } else if ((reg.style == FILLED) || (reg.style == TRANSPARENT)) {
          opacity = (reg.style == TRANSPARENT) ? reg.opacity : 1.0;
          myreg->forward_transform(center[0], center[1], center[2]);
          image->draw_sphere(center, reg.color, 2.0 * myreg->radius, opacity);
        }

      } else {
        if (comm->me == 0)
          error->warning(FLERR, "Region style {} is not yet supported by dump image", regstyle);
      }
    }
  }
  // clang-format off

  // render outline of my sub-box, orthogonal or triclinic

  if (subboxflag) {
    diameter = MIN(boxxhi-boxxlo,boxyhi-boxylo);
    if (domain->dimension == 3) diameter = MIN(diameter,boxzhi-boxzlo);
    diameter *= subboxdiam;

    double *sublo = domain->sublo;
    double *subhi = domain->subhi;

    double (*boxcorners)[3];
    double box[8][3];
    if (domain->triclinic == 0) {
      box[0][0] = sublo[0]; box[0][1] = sublo[1]; box[0][2] = sublo[2];
      box[1][0] = subhi[0]; box[1][1] = sublo[1]; box[1][2] = sublo[2];
      box[2][0] = sublo[0]; box[2][1] = subhi[1]; box[2][2] = sublo[2];
      box[3][0] = subhi[0]; box[3][1] = subhi[1]; box[3][2] = sublo[2];
      box[4][0] = sublo[0]; box[4][1] = sublo[1]; box[4][2] = subhi[2];
      box[5][0] = subhi[0]; box[5][1] = sublo[1]; box[5][2] = subhi[2];
      box[6][0] = sublo[0]; box[6][1] = subhi[1]; box[6][2] = subhi[2];
      box[7][0] = subhi[0]; box[7][1] = subhi[1]; box[7][2] = subhi[2];
      boxcorners = box;
    } else {
      domain->subbox_corners();
      boxcorners = domain->corners;
    }

    image->draw_box(boxcorners,diameter,subboxopacity);
  }

  // render outline of simulation box, orthogonal or triclinic

  if (boxflag) {
    diameter = MIN(boxxhi-boxxlo,boxyhi-boxylo);
    if (domain->dimension == 3) diameter = MIN(diameter,boxzhi-boxzlo);
    diameter *= boxdiam;

    double (*boxcorners)[3];
    double box[8][3];
    if (domain->triclinic == 0) {
      box[0][0] = boxxlo; box[0][1] = boxylo; box[0][2] = boxzlo;
      box[1][0] = boxxhi; box[1][1] = boxylo; box[1][2] = boxzlo;
      box[2][0] = boxxlo; box[2][1] = boxyhi; box[2][2] = boxzlo;
      box[3][0] = boxxhi; box[3][1] = boxyhi; box[3][2] = boxzlo;
      box[4][0] = boxxlo; box[4][1] = boxylo; box[4][2] = boxzhi;
      box[5][0] = boxxhi; box[5][1] = boxylo; box[5][2] = boxzhi;
      box[6][0] = boxxlo; box[6][1] = boxyhi; box[6][2] = boxzhi;
      box[7][0] = boxxhi; box[7][1] = boxyhi; box[7][2] = boxzhi;
      boxcorners = box;
    } else {
      domain->box_corners();
      boxcorners = domain->corners;
    }

    image->draw_box(boxcorners,diameter,boxopacity);
  }

  // render XYZ axes in red/green/blue
  // offset by 10% of box size and scale by axeslen

  if (axesflag) {
    diameter = MIN(boxxhi-boxxlo,boxyhi-boxylo);
    if (domain->dimension == 3) diameter = MIN(diameter,boxzhi-boxzlo);
    diameter *= axesdiam;

    double (*boxcorners)[3];
    double axes[4][3];
    if (domain->triclinic == 0) {
      axes[0][0] = boxxlo; axes[0][1] = boxylo; axes[0][2] = boxzlo;
      axes[1][0] = boxxhi; axes[1][1] = boxylo; axes[1][2] = boxzlo;
      axes[2][0] = boxxlo; axes[2][1] = boxyhi; axes[2][2] = boxzlo;
      axes[3][0] = boxxlo; axes[3][1] = boxylo; axes[3][2] = boxzhi;
    } else {
      domain->box_corners();
      boxcorners = domain->corners;
      axes[0][0] = boxcorners[0][0];
      axes[0][1] = boxcorners[0][1];
      axes[0][2] = boxcorners[0][2];
      axes[1][0] = boxcorners[1][0];
      axes[1][1] = boxcorners[1][1];
      axes[1][2] = boxcorners[1][2];
      axes[2][0] = boxcorners[2][0];
      axes[2][1] = boxcorners[2][1];
      axes[2][2] = boxcorners[2][2];
      axes[3][0] = boxcorners[4][0];
      axes[3][1] = boxcorners[4][1];
      axes[3][2] = boxcorners[4][2];
    }

    double offset = MAX(boxxhi-boxxlo,boxyhi-boxylo);
    if (domain->dimension == 3) offset = MAX(offset,boxzhi-boxzlo);

    if (axesflag == CENTER)
      offset *= 0.5;
    else
      offset *= 0.1;
    if (axesflag == CENTER) {
      axes[0][0] += offset; axes[0][1] += offset; axes[0][2] += offset;
      axes[1][0] += offset; axes[1][1] += offset; axes[1][2] += offset;
      axes[2][0] += offset; axes[2][1] += offset; axes[2][2] += offset;
      axes[3][0] += offset; axes[3][1] += offset; axes[3][2] += offset;
    } else if (axesflag == LOWERLEFT) {
      axes[0][0] -= offset; axes[0][1] -= offset; axes[0][2] -= offset;
      axes[1][0] -= offset; axes[1][1] -= offset; axes[1][2] -= offset;
      axes[2][0] -= offset; axes[2][1] -= offset; axes[2][2] -= offset;
      axes[3][0] -= offset; axes[3][1] -= offset; axes[3][2] -= offset;
    } else if ((domain->dimension == 3) && (axesflag == LOWERRIGHT)) {
      axes[0][0] -= offset; axes[0][1] += 8.0*offset; axes[0][2] -= offset;
      axes[1][0] -= offset; axes[1][1] += 8.0*offset; axes[1][2] -= offset;
      axes[2][0] -= offset; axes[2][1] += 8.0*offset; axes[2][2] -= offset;
      axes[3][0] -= offset; axes[3][1] += 8.0*offset; axes[3][2] -= offset;
    } else if  ((domain->dimension == 3) && (axesflag == UPPERLEFT)) {
      axes[0][0] -= offset; axes[0][1] -= offset; axes[0][2] += 8.0*offset;
      axes[1][0] -= offset; axes[1][1] -= offset; axes[1][2] += 8.0*offset;
      axes[2][0] -= offset; axes[2][1] -= offset; axes[2][2] += 8.0*offset;
      axes[3][0] -= offset; axes[3][1] -= offset; axes[3][2] += 8.0*offset;
    } else if  ((domain->dimension == 3) && (axesflag == UPPERRIGHT)) {
      axes[0][0] -= offset; axes[0][1] += 8.0*offset; axes[0][2] += 8.0*offset;
      axes[1][0] -= offset; axes[1][1] += 8.0*offset; axes[1][2] += 8.0*offset;
      axes[2][0] -= offset; axes[2][1] += 8.0*offset; axes[2][2] += 8.0*offset;
      axes[3][0] -= offset; axes[3][1] += 8.0*offset; axes[3][2] += 8.0*offset;
    } else if ((domain->dimension == 2) && (axesflag == LOWERRIGHT)) {
      axes[0][0] += 8.0*offset; axes[0][1] -= offset; axes[0][2] -= offset;
      axes[1][0] += 8.0*offset; axes[1][1] -= offset; axes[1][2] -= offset;
      axes[2][0] += 8.0*offset; axes[2][1] -= offset; axes[2][2] -= offset;
      axes[3][0] += 8.0*offset; axes[3][1] -= offset; axes[3][2] -= offset;
    } else if  ((domain->dimension == 2) && (axesflag == UPPERLEFT)) {
      axes[0][0] -= offset; axes[0][1] += 8.0*offset; axes[0][2] -= offset;
      axes[1][0] -= offset; axes[1][1] += 8.0*offset; axes[1][2] -= offset;
      axes[2][0] -= offset; axes[2][1] += 8.0*offset; axes[2][2] -= offset;
      axes[3][0] -= offset; axes[3][1] += 8.0*offset; axes[3][2] -= offset;
    } else if  ((domain->dimension == 2) && (axesflag == UPPERRIGHT)) {
      axes[0][0] += 8.0*offset; axes[0][1] += 8.0*offset; axes[0][2] -= offset;
      axes[1][0] += 8.0*offset; axes[1][1] += 8.0*offset; axes[1][2] -= offset;
      axes[2][0] += 8.0*offset; axes[2][1] += 8.0*offset; axes[2][2] -= offset;
      axes[3][0] += 8.0*offset; axes[3][1] += 8.0*offset; axes[3][2] -= offset;
    }

    axes[1][0] = axes[0][0] + axeslen*(axes[1][0]-axes[0][0]);
    axes[1][1] = axes[0][1] + axeslen*(axes[1][1]-axes[0][1]);
    axes[1][2] = axes[0][2] + axeslen*(axes[1][2]-axes[0][2]);
    axes[2][0] = axes[0][0] + axeslen*(axes[2][0]-axes[0][0]);
    axes[2][1] = axes[0][1] + axeslen*(axes[2][1]-axes[0][1]);
    axes[2][2] = axes[0][2] + axeslen*(axes[2][2]-axes[0][2]);
    axes[3][0] = axes[0][0] + axeslen*(axes[3][0]-axes[0][0]);
    axes[3][1] = axes[0][1] + axeslen*(axes[3][1]-axes[0][1]);
    axes[3][2] = axes[0][2] + axeslen*(axes[3][2]-axes[0][2]);

    image->draw_axes(axes,diameter,axesopacity);
  }
}

/* ---------------------------------------------------------------------- */

void DumpImage::grid_cell_corners_2d(int ix, int iy)
{
  double *boxlo = domain->boxlo;
  double *prd = domain->prd;

  if (!domain->triclinic) {
    double xdelta = prd[0] / nxgrid;
    double ydelta = prd[1] / nygrid;

    int n = 0;
    for (int y = 0; y < 2; y++)
      for (int x = 0; x < 2; x++) {
        gcorners[n][0] = boxlo[0] + (ix+x) * xdelta;
        gcorners[n][1] = boxlo[1] + (iy+y) * ydelta;
        gcorners[n][2] = 0.0;
        n++;
      }

  } else {
    double lamda[3];
    lamda[2] = 0.0;

    double dx = 1.0 / nxgrid;
    double dy = 1.0 / nygrid;

    int n = 0;
    for (int y = 0; y < 2; y++)
      for (int x = 0; x < 2; x++) {
        lamda[0] = (ix+x) * dx;
        lamda[1] = (iy+y) * dy;
        domain->lamda2x(lamda,gcorners[n]);
        n++;
      }
  }
}

/* ---------------------------------------------------------------------- */

void DumpImage::grid_cell_corners_3d(int ix, int iy, int iz)
{
  double *boxlo = domain->boxlo;
  double *prd = domain->prd;

  if (!domain->triclinic) {
    double xdelta = prd[0] / nxgrid;
    double ydelta = prd[1] / nygrid;
    double zdelta = prd[2] / nzgrid;

    int n = 0;
    for (int z = 0; z < 2; z++)
      for (int y = 0; y < 2; y++)
        for (int x = 0; x < 2; x++) {
          gcorners[n][0] = boxlo[0] + (ix+x) * xdelta;
          gcorners[n][1] = boxlo[1] + (iy+y) * ydelta;
          gcorners[n][2] = boxlo[2] + (iz+z) * zdelta;
          n++;
        }

  } else {
    double lamda[3];

    double dx = 1.0 / nxgrid;
    double dy = 1.0 / nygrid;
    double dz = 1.0 / nzgrid;

    int n = 0;
    for (int z = 0; z < 2; z++)
      for (int y = 0; y < 2; y++)
        for (int x = 0; x < 2; x++) {
          lamda[0] = (ix+x) * dx;
          lamda[1] = (iy+y) * dy;
          lamda[2] = (iz+z) * dz;
          domain->lamda2x(lamda,gcorners[n]);
          n++;
        }
  }
}

/* ---------------------------------------------------------------------- */

int DumpImage::pack_forward_comm(int n, int *list, double *buf,
                                 int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,m;

  m = 0;

  if (comm_forward == 1) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = chooseghost[j];
    }
  } else {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = chooseghost[j];
      buf[m++] = bufcopy[j][0];
      buf[m++] = bufcopy[j][1];
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void DumpImage::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;

  if (comm_forward == 1)
    for (i = first; i < last; i++) chooseghost[i] = static_cast<int>(buf[m++]);
  else {
    for (i = first; i < last; i++) {
      chooseghost[i] = static_cast<int>(buf[m++]);
      bufcopy[i][0] = buf[m++];
      bufcopy[i][1] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

void *DumpImage::extract(const char *str, int &dim)
{
  dim=0;
  if (strcmp(str,"image") == 0) {
    return image;
  }
  return nullptr;
}

/* ---------------------------------------------------------------------- */

int DumpImage::modify_param(int narg, char **arg)
{
  int n = DumpCustom::modify_param(narg,arg);
  if (n) return n;

  // determine offset in list of arguments for error pointer. also handle the no match case.
  int argoff = 0;
  while (input && input->arg[argoff] && input->arg[argoff+1]) {
    if ((strcmp(input->arg[argoff], arg[0]) == 0) && (strcmp(input->arg[argoff+1], arg[1]) == 0))
      break;
    ++argoff;
  }

  if (strcmp(arg[0],"acolor") == 0) {
    if (narg < 3) utils::missing_cmd_args(FLERR, "dump_modify acolor", error);
    int nlo,nhi;
    utils::bounds_typelabel(FLERR,arg[1],1,atom->ntypes,nlo,nhi,lmp,Atom::ATOM);

    // get list of colors
    // assign colors in round-robin fashion to types

    auto colors = Tokenizer(arg[2],"/").as_vector();
    const int ncolors = colors.size();

    int m = 0;
    for (int i = nlo; i <= nhi; i++) {
      colortype[i] = image->color2rgb(colors[m%ncolors].c_str());
      if (colortype[i] == nullptr)
        error->all(FLERR,argoff+2,"Invalid color in dump_modify acolor command {}", arg[2]);
      m++;
    }
    return 3;
  }

  if (strcmp(arg[0],"adiam") == 0) {
    if (narg < 3) utils::missing_cmd_args(FLERR, "dump_modify adiam", error);
    int nlo,nhi;
    utils::bounds_typelabel(FLERR,arg[1],1,atom->ntypes,nlo,nhi,lmp,Atom::ATOM);
    double diam = utils::numeric(FLERR,arg[2],false,lmp);
    if (diam <= 0.0)
      error->all(FLERR,argoff+2,"Illegal diameter {} in dump_modify adiam command", diam);
    for (int i = nlo; i <= nhi; i++) diamtype[i] = diam;
    return 3;
  }

  if (strcmp(arg[0],"atrans") == 0) {
    if (narg < 3) utils::missing_cmd_args(FLERR, "dump_modify atrans", error);
    int nlo,nhi;
    utils::bounds_typelabel(FLERR,arg[1],1,atom->ntypes,nlo,nhi,lmp,Atom::ATOM);
    double opacity = utils::numeric(FLERR,arg[2],false,lmp);
    if ((opacity < 0.0) || (opacity > 1.0))
      error->all(FLERR,argoff+2,"Illegal transparency {} in dump_modify atrans command", opacity);
    for (int i = nlo; i <= nhi; i++) aopacity[i] = opacity;
    return 3;
  }

  if ((strcmp(arg[0],"amap") == 0) || (strcmp(arg[0],"gmap") == 0)) {
    if (narg < 6) utils::missing_cmd_args(FLERR, "dump_modify amap/gmap", error);
    if (strlen(arg[3]) != 2)
      error->all(FLERR,argoff+3, "Incorrect dump_modify amap/gmap colormap style {}", arg[3]);
    int factor = 0;
    if (arg[3][0] == 's') factor = 1;
    else if (arg[3][0] == 'c') factor = 2;
    else if (arg[3][0] == 'd') factor = 3;
    else error->all(FLERR,argoff+3, "Unknown dump_modify amap/gmap color map type {}", arg[3][0]);
    int nentry = utils::inumeric(FLERR,arg[5],false,lmp);
    if (nentry < 1)
      error->all(FLERR, argoff+5, "Must have at least 1 color map entry for dump_modify amap/gmap");
    n = 6 + factor*nentry;
    if (narg < n)  utils::missing_cmd_args(FLERR, "dump_modify amap/gmap", error);
    int flag = 0;
    if (strcmp(arg[0], "amap") == 0) flag = image->map_reset(0, n-1, &arg[1]);
    if (strcmp(arg[0], "gmap") == 0) flag = image->map_reset(1, n-1, &arg[1]);
    if (flag)
      error->all(FLERR, argoff+flag, "Invalid map settings in dump_modify amap/gmap command");

    return n;
  }

  if (strcmp(arg[0],"bcolor") == 0) {
    if (narg < 3) utils::missing_cmd_args(FLERR, "dump_modify bcolor", error);
    if (atom->nbondtypes == 0)
      error->all(FLERR,"Dump modify bcolor not allowed with no bond types");
    // ignore if bonds are not displayed
    if (bondflag == NO) return 3;
    int nlo,nhi;
    utils::bounds_typelabel(FLERR,arg[1],1,atom->nbondtypes,nlo,nhi,lmp,Atom::BOND);

    // process list of ncount colornames separated by '/'
    // assign colors in round-robin fashion to bond types

    auto colors = Tokenizer(arg[2],"/").as_vector();
    const int ncolors = colors.size();

    int m = 0;
    for (int i = nlo; i <= nhi; i++) {
      bcolortype[i] = image->color2rgb(colors[m%ncolors].c_str());
      if (bcolortype[i] == nullptr)
        error->all(FLERR, argoff + 2, "Invalid color in dump_modify bcolor command");
      m++;
    }
    return 3;
  }

  if (strcmp(arg[0],"bdiam") == 0) {
    if (narg < 3) utils::missing_cmd_args(FLERR, "dump_modify bdiam", error);
    if (atom->nbondtypes == 0)
      error->all(FLERR, argoff, "Dump modify bdiam not allowed with no bond types");
    // ignore if bonds are not displayed
    if (bondflag == NO) return 3;
    int nlo,nhi;
    utils::bounds_typelabel(FLERR,arg[1],1,atom->nbondtypes,nlo,nhi,lmp,Atom::BOND);
    double diam = utils::numeric(FLERR,arg[2],false,lmp);
    if (diam <= 0.0)
      error->all(FLERR, argoff + 2, "Invalid bdiam diameter in dump_modify command");
    for (int i = nlo; i <= nhi; i++) bdiamtype[i] = diam;
    return 3;
  }

  if (strcmp(arg[0],"btrans") == 0) {
    if (narg < 3) utils::missing_cmd_args(FLERR, "dump_modify btrans", error);
    if (atom->nbondtypes == 0)
      error->all(FLERR,"Dump modify btrans not allowed with no bond types");
    // ignore if bonds are not displayed
    if (bondflag == NO) return 3;
    int nlo,nhi;
    utils::bounds_typelabel(FLERR,arg[1],1,atom->nbondtypes,nlo,nhi,lmp,Atom::BOND);
    double opacity = utils::numeric(FLERR,arg[2],false,lmp);
    if ((opacity < 0.0) || (opacity > 1.0))
      error->all(FLERR, argoff + 2, "Invalid btrans opacity in dump_modify command");
    for (int i = nlo; i <= nhi; i++) bopacity[i] = opacity;
    return 3;
  }

  if (strcmp(arg[0],"backcolor") == 0) {
    if (narg < 2) utils::missing_cmd_args(FLERR, "dump_modify backcolor", error);
    double *color = image->color2rgb(arg[1]);
    if (color == nullptr) error->all(FLERR,"Invalid color in dump_modify command");
    image->background[0] = static_cast<int>(color[0]*255.0);
    image->background[1] = static_cast<int>(color[1]*255.0);
    image->background[2] = static_cast<int>(color[2]*255.0);
    return 2;
  }

  if (strcmp(arg[0],"backcolor2") == 0) {
    if (narg < 2) utils::missing_cmd_args(FLERR, "dump_modify backcolor2", error);
    if (strcmp(arg[1],"none") == 0) {
      image->background2[0] = -1;
      image->background2[1] = -1;
      image->background2[2] = -1;
    } else {
      double *color = image->color2rgb(arg[1]);
      if (color == nullptr)
        error->all(FLERR, argoff + 1, "Invalid color in dump_modify backcolor2 command");
      image->background2[0] = static_cast<int>(color[0]*255.0);
      image->background2[1] = static_cast<int>(color[1]*255.0);
      image->background2[2] = static_cast<int>(color[2]*255.0);
    }
    return 2;
  }

  if (strcmp(arg[0],"boxcolor") == 0) {
    if (narg < 2) utils::missing_cmd_args(FLERR, "dump_modify boxcolor", error);
    image->boxcolor = image->color2rgb(arg[1]);
    if (image->boxcolor == nullptr)
      error->all(FLERR, argoff + 1, "Invalid color in dump_modify boxcolor command");
    return 2;
  }

  if (strcmp(arg[0],"boxtrans") == 0) {
    if (narg < 2) utils::missing_cmd_args(FLERR, "dump_modify boxtrans", error);
    boxopacity = utils::numeric(FLERR,arg[1],false,lmp);
    if ((boxopacity < 0.0) || (boxopacity > 1.0))
      error->all(FLERR, argoff + 1, "Invalid boxtrans opacity in dump_modify command");
    return 2;
  }

  if (strcmp(arg[0],"subboxtrans") == 0) {
    if (narg < 2) utils::missing_cmd_args(FLERR, "dump_modify subboxtrans", error);
    subboxopacity = utils::numeric(FLERR,arg[1],false,lmp);
    if ((subboxopacity < 0.0) || (subboxopacity > 1.0))
      error->all(FLERR, argoff + 1, "Invalid subboxtrans opacity in dump_modify command");
    return 2;
  }

  if (strcmp(arg[0],"axestrans") == 0) {
    if (narg < 2) utils::missing_cmd_args(FLERR, "dump_modify axestrans", error);
    axesopacity = utils::numeric(FLERR,arg[1],false,lmp);
    if ((axesopacity < 0.0) || (axesopacity > 1.0))
      error->all(FLERR, argoff + 1, "Invalid axestrans opacity in dump_modify command");
    return 2;
  }

  if (strcmp(arg[0],"color") == 0) {
    if (narg < 3) utils::missing_cmd_args(FLERR, "dump_modify color", error);
    if (utils::strmatch(arg[2], "^0x[0-9a-fA-F]+$")) {
      char *ptr = nullptr;
      auto val = strtol(arg[2], &ptr, 16);
      double rval =  ((val & 0xff0000) >> 16) / 255.0;
      double gval =  ((val & 0x00ff00) >> 8) / 255.0;
      double bval =  (val & 0x0000ff) / 255.0;
      int flag = image->addcolor(arg[1], rval, gval, bval);
      if (flag) error->all(FLERR, argoff + 1 + flag, "Incorrect dump_modify color command");
      return 3;
    } else {
      if (narg < 5) utils::missing_cmd_args(FLERR, "dump_modify color", error);
      int flag = image->addcolor(arg[1],utils::numeric(FLERR,arg[2],false,lmp),
                                 utils::numeric(FLERR,arg[3],false,lmp),
                                 utils::numeric(FLERR,arg[4],false,lmp));
      if (flag) error->all(FLERR, argoff + 1 + flag, "Incorrect dump_modify color command");
      return 5;
    }
  }

  if (strcmp(arg[0],"ccolor") == 0) {
    if (narg < 3) utils::missing_cmd_args(FLERR, "dump_modify ccolor", error);
    auto *color =  image->color2rgb(arg[2]);
    if (!color) error->all(FLERR, argoff + 2, "Unknown color for dump_modify ccolor: {}", arg[2]);
    bool match = false;
    for (auto &iobj : objects) {
      if (iobj.cptr && (iobj.id == arg[1])) {
        iobj.rgb = color;
        match = true;
      }
    }
    if (!match) error->all(FLERR, argoff + 1, "Compute ID {} is not included in dump {}", arg[1], id);
    return 3;
  }

  if (strcmp(arg[0],"ctrans") == 0) {
    if (narg < 3) utils::missing_cmd_args(FLERR, "dump_modify ctrans", error);
    double opacity = utils::numeric(FLERR,arg[2],false,lmp);
    if ((opacity < 0.0) || (opacity > 1.0))
      error->all(FLERR, argoff + 2, "Illegal transparency {} for dump_modify ctrans", arg[2]);
    bool match = false;
    for (auto &iobj : objects) {
      if (iobj.cptr && (iobj.id == arg[1])) {
        iobj.opacity = opacity;
        match = true;
      }
    }
    if (!match) error->all(FLERR, argoff + 1, "Compute ID {} is not included in dump {}", arg[1], id);
    return 3;
  }

  if (strcmp(arg[0],"fcolor") == 0) {
    if (narg < 3) utils::missing_cmd_args(FLERR, "dump_modify fcolor", error);
    auto *color =  image->color2rgb(arg[2]);
    if (!color) error->all(FLERR, argoff + 2, "Unknown color for dump_modify fcolor: {}", arg[2]);
    bool match = false;
    for (auto &iobj : objects) {
      if (iobj.fptr && (iobj.id == arg[1])) {
        iobj.rgb = color;
        match = true;
      }
    }
    if (!match) error->all(FLERR, argoff + 1, "Fix ID {} is not included in dump {}", arg[1], id);
    return 3;
  }

  if (strcmp(arg[0],"ftrans") == 0) {
    if (narg < 3) utils::missing_cmd_args(FLERR, "dump_modify ftrans", error);
    double opacity = utils::numeric(FLERR,arg[2],false,lmp);
    if ((opacity < 0.0) || (opacity > 1.0))
      error->all(FLERR, argoff + 2, "Illegal transparency {} for dump_modify ftrans", arg[2]);
    bool match = false;
    for (auto &iobj : objects) {
      if (iobj.fptr && (iobj.id == arg[1])) {
        iobj.opacity = opacity;
        match = true;
      }
    }
    if (!match) error->all(FLERR, argoff + 1, "Fix ID {} is not included in dump {}", arg[1], id);
    return 3;
  }

  return 0;
}
