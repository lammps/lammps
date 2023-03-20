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
#include "atom_vec_line.h"
#include "atom_vec_tri.h"
#include "body.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "grid2d.h"
#include "grid3d.h"
#include "image.h"
#include "input.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "molecule.h"
#include "tokenizer.h"
#include "update.h"
#include "variable.h"

#include <cmath>
#include <cctype>
#include <cstring>

using namespace LAMMPS_NS;
using MathConst::DEG2RAD;

#define BIG 1.0e20

enum{NUMERIC,ATOM,TYPE,ELEMENT,ATTRIBUTE};
enum{SPHERE,LINE,TRI};           // also in some Body and Fix child classes
enum{STATIC,DYNAMIC};
enum{NO=0,YES=1};

/* ---------------------------------------------------------------------- */

DumpImage::DumpImage(LAMMPS *lmp, int narg, char **arg) :
  DumpCustom(lmp, narg, arg), thetastr(nullptr), phistr(nullptr), cxstr(nullptr),
  cystr(nullptr), czstr(nullptr), upxstr(nullptr), upystr(nullptr), upzstr(nullptr),
  zoomstr(nullptr), diamtype(nullptr), diamelement(nullptr),
  bdiamtype(nullptr), colortype(nullptr), colorelement(nullptr), bcolortype(nullptr),
  avec_line(nullptr), avec_tri(nullptr), avec_body(nullptr), fixptr(nullptr), image(nullptr),
  chooseghost(nullptr), bufcopy(nullptr)
{
  if (binary || multiproc) error->all(FLERR,"Invalid dump image filename");

  // force binary flag on to avoid corrupted output on Windows

  binary = 1;
  multifile_override = 0;

  // flag has_id as true to avoid bogus warnings about atom IDs for dump styles derived from DumpCustom

  has_id = true;

  // set filetype based on filename suffix

  int n = strlen(filename);
  if (strlen(filename) > 4 && strcmp(&filename[n-4],".jpg") == 0)
    filetype = JPG;
  else if (strlen(filename) > 4 && strcmp(&filename[n-4],".JPG") == 0)
    filetype = JPG;
  else if (strlen(filename) > 5 && strcmp(&filename[n-5],".jpeg") == 0)
    filetype = JPG;
  else if (strlen(filename) > 5 && strcmp(&filename[n-5],".JPEG") == 0)
    filetype = JPG;
  else if (strlen(filename) > 4 && strcmp(&filename[n-4],".png") == 0)
    filetype = PNG;
  else if (strlen(filename) > 4 && strcmp(&filename[n-4],".PNG") == 0)
    filetype = PNG;
  else filetype = PPM;

#ifndef LAMMPS_JPEG
  if (filetype == JPG)
    error->all(FLERR,"Support for writing images in JPEG format not included");
#endif
#ifndef LAMMPS_PNG
  if (filetype == PNG)
    error->all(FLERR,"Support for writing images in PNG format not included");
#endif

  // atom color,diameter settings

  if (nfield != 2) error->all(FLERR,"Illegal dump image command");

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
  lineflag = triflag = bodyflag = fixflag = NO;
  id_grid_compute = id_grid_fix = nullptr;

  if (atom->nbondtypes == 0) bondflag = NO;
  else {
    bondflag = YES;
    bcolor = ATOM;
    bdiam = NUMERIC;
    bdiamvalue = 0.5;
  }
  char *fixID = nullptr;

  thetastr = phistr = nullptr;
  cflag = STATIC;
  cx = cy = cz = 0.5;
  cxstr = cystr = czstr = nullptr;

  upxstr = upystr = upzstr = nullptr;
  zoomstr = nullptr;
  boxflag = YES;
  boxdiam = 0.02;
  axesflag = NO;
  subboxflag = NO;

  // parse optional args

  int iarg = ioptional;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"atom") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump image command");
      atomflag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;

    } else if (strcmp(arg[iarg],"adiam") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump image command");
      adiam = NUMERIC;
      adiamvalue = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (adiamvalue <= 0.0) error->all(FLERR,"Illegal dump image command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"bond") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal dump image command");
      if (atom->nbondtypes == 0)
        error->all(FLERR,"Dump image bond not allowed with no bond types");
      bondflag = YES;
      if (strcmp(arg[iarg+1],"none") == 0) bondflag = NO;
      else if (strcmp(arg[iarg+1],"atom") == 0) bcolor = ATOM;
      else if (strcmp(arg[iarg+1],"type") == 0) bcolor = TYPE;
      else error->all(FLERR,"Illegal dump image command");
      if (!islower(arg[iarg+2][0])) {
          bdiam = NUMERIC;
          bdiamvalue = utils::numeric(FLERR,arg[iarg+2],false,lmp);
          if (bdiamvalue <= 0.0) error->all(FLERR,"Illegal dump image command");
      } else if (strcmp(arg[iarg+2],"atom") == 0) bdiam = ATOM;
      else if (strcmp(arg[iarg+2],"type") == 0) bdiam = TYPE;
      else if (strcmp(arg[iarg+2],"none") == 0) bondflag = NO;
      else error->all(FLERR,"Illegal dump image command");
      iarg += 3;

    } else if (strcmp(arg[iarg],"grid") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump image command");
      gridflag = YES;

      char *id;
      int igrid,idata,index;
      int iflag =
        utils::check_grid_reference((char *) "Dump image",
                                    arg[iarg+1],nevery,id,
                                    igrid,idata,index,lmp);
      if (iflag < 0) error->all(FLERR,"Invalid grid reference in dump image command");

      if (iflag == ArgInfo::COMPUTE) id_grid_compute = utils::strdup(id);
      else if (iflag == ArgInfo::FIX) id_grid_fix = utils::strdup(id);
      delete [] id;
      grid_igrid = igrid;
      grid_idata = idata;
      grid_index = index;
      iarg += 2;

    } else if (strcmp(arg[iarg],"line") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal dump image command");
      lineflag = YES;
      if (strcmp(arg[iarg+1],"type") == 0) lcolor = TYPE;
      else error->all(FLERR,"Illegal dump image command");
      ldiam = NUMERIC;
      ldiamvalue = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      iarg += 3;

    } else if (strcmp(arg[iarg],"tri") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal dump image command");
      triflag = YES;
      if (strcmp(arg[iarg+1],"type") == 0) tcolor = TYPE;
      else error->all(FLERR,"Illegal dump image command");
      tstyle = utils::inumeric(FLERR,arg[iarg+2],false,lmp);
      tdiamvalue = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      iarg += 4;

    } else if (strcmp(arg[iarg],"body") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal dump image command");
      bodyflag = YES;
      if (strcmp(arg[iarg+1],"type") == 0) bodycolor = TYPE;
      else error->all(FLERR,"Illegal dump image command");
      bodyflag1 = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      bodyflag2 = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      iarg += 4;

    } else if (strcmp(arg[iarg],"fix") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal dump image command");
      fixflag = YES;
      fixID = arg[iarg+1];
      if (strcmp(arg[iarg+2],"type") == 0) fixcolor = TYPE;
      else error->all(FLERR,"Illegal dump image command");
      fixflag1 = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      fixflag2 = utils::numeric(FLERR,arg[iarg+4],false,lmp);
      iarg += 5;

    } else if (strcmp(arg[iarg],"size") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal dump image command");
      int width = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      int height = utils::inumeric(FLERR,arg[iarg+2],false,lmp);
      if (width <= 0 || height <= 0)
        error->all(FLERR,"Illegal dump image command");
      image->width = width;
      image->height = height;
      iarg += 3;

    } else if (strcmp(arg[iarg],"view") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal dump image command");
      if (utils::strmatch(arg[iarg+1],"^v_")) {
        thetastr = utils::strdup(arg[iarg+1]+2);
      } else {
        const double theta = utils::numeric(FLERR,arg[iarg+1],false,lmp);
        if (theta < 0.0 || theta > 180.0)
          error->all(FLERR,"Invalid dump image theta value");
        image->theta = DEG2RAD * theta;
      }
      if (utils::strmatch(arg[iarg+2],"^v_")) {
        phistr = utils::strdup(arg[iarg+2]+2);
      } else {
        image->phi = DEG2RAD * utils::numeric(FLERR,arg[iarg+2],false,lmp);
      }
      iarg += 3;

    } else if (strcmp(arg[iarg],"center") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal dump image command");
      if (strcmp(arg[iarg+1],"s") == 0) cflag = STATIC;
      else if (strcmp(arg[iarg+1],"d") == 0) cflag = DYNAMIC;
      else error->all(FLERR,"Illegal dump image command");
      if (utils::strmatch(arg[iarg+2],"^v_")) {
        cxstr = utils::strdup(arg[iarg+2]+2);
        cflag = DYNAMIC;
      } else cx = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      if (utils::strmatch(arg[iarg+3],"^v_")) {
        cystr = utils::strdup(arg[iarg+3]+2);
        cflag = DYNAMIC;
      } else cy = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      if (utils::strmatch(arg[iarg+4],"^v_")) {
        czstr = utils::strdup(arg[iarg+4]+2);
        cflag = DYNAMIC;
      } else cz = utils::numeric(FLERR,arg[iarg+4],false,lmp);
      iarg += 5;

    } else if (strcmp(arg[iarg],"up") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal dump image command");
      if (utils::strmatch(arg[iarg+1],"^v_")) {
        upxstr = utils::strdup(arg[iarg+1]+2);
      } else image->up[0] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (utils::strmatch(arg[iarg+2],"^v_")) {
        upystr = utils::strdup(arg[iarg+2]+2);
      } else image->up[1] = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      if (utils::strmatch(arg[iarg+3],"^v_")) {
        upzstr = utils::strdup(arg[iarg+3]+2);
      } else image->up[2] = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      iarg += 4;

    } else if (strcmp(arg[iarg],"zoom") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump image command");
      if (utils::strmatch(arg[iarg+1],"^v_")) {
        zoomstr = utils::strdup(arg[iarg+1]+2);
      } else {
        double zoom = utils::numeric(FLERR,arg[iarg+1],false,lmp);
        if (zoom <= 0.0) error->all(FLERR,"Illegal dump image command");
        image->zoom = zoom;
      }
      iarg += 2;

    } else if (strcmp(arg[iarg],"box") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal dump image command");
      boxflag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      boxdiam = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      if (boxdiam < 0.0) error->all(FLERR,"Illegal dump image command");
      iarg += 3;

    } else if (strcmp(arg[iarg],"axes") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal dump image command");
      axesflag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      axeslen = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      axesdiam = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      if (axeslen < 0.0 || axesdiam < 0.0)
        error->all(FLERR,"Illegal dump image command");
      iarg += 4;

    } else if (strcmp(arg[iarg],"subbox") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal dump image command");
      subboxflag = utils::logical(FLERR,arg[iarg+1],false,lmp);
      subboxdiam = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      if (subboxdiam < 0.0) error->all(FLERR,"Illegal dump image command");
      iarg += 3;

    } else if (strcmp(arg[iarg],"shiny") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump image command");
      double shiny = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (shiny < 0.0 || shiny > 1.0)
        error->all(FLERR,"Illegal dump image command");
      image->shiny = shiny;
      iarg += 2;

    } else if (strcmp(arg[iarg],"ssao") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal dump image command");
      image->ssao = utils::logical(FLERR,arg[iarg+1],false,lmp);
      int seed = utils::inumeric(FLERR,arg[iarg+2],false,lmp);
      if (seed <= 0) error->all(FLERR,"Illegal dump image command");
      image->seed = seed;
      double ssaoint = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      if (ssaoint < 0.0 || ssaoint > 1.0)
        error->all(FLERR,"Illegal dump image command");
      image->ssaoint = ssaoint;
      iarg += 4;

    } else error->all(FLERR,"Illegal dump image command");
  }

  // error checks and setup for lineflag, triflag, bodyflag, fixflag

  if (lineflag) {
    avec_line = dynamic_cast<AtomVecLine *>(atom->style_match("line"));
    if (!avec_line)
      error->all(FLERR,"Dump image line requires atom style line");
  }
  if (triflag) {
    avec_tri = dynamic_cast<AtomVecTri *>(atom->style_match("tri"));
    if (!avec_tri)
      error->all(FLERR,"Dump image tri requires atom style tri");
  }
  if (bodyflag) {
    avec_body = dynamic_cast<AtomVecBody *>(atom->style_match("body"));
    if (!avec_body)
      error->all(FLERR,"Dump image body yes requires atom style body");
  }

  extraflag = 0;
  if (lineflag || triflag || bodyflag) extraflag = 1;

  if (fixflag) {
    fixptr = modify->get_fix_by_id(fixID);
    if (!fixptr) error->all(FLERR,"Fix ID {} for dump image does not exist", fixID);

  }

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

  for (int i = 1; i <= ntypes; i++) {
    diamtype[i] = 1.0;
    if (i % 6 == 1) colortype[i] = image->color2rgb("red");
    else if (i % 6 == 2) colortype[i] = image->color2rgb("green");
    else if (i % 6 == 3) colortype[i] = image->color2rgb("blue");
    else if (i % 6 == 4) colortype[i] = image->color2rgb("yellow");
    else if (i % 6 == 5) colortype[i] = image->color2rgb("aqua");
    else if (i % 6 == 0) colortype[i] = image->color2rgb("cyan");
  }

  if (bondflag) {
    bdiamtype = new double[atom->nbondtypes+1];
    bcolortype = new double*[atom->nbondtypes+1];
    for (int i = 1; i <= atom->nbondtypes; i++) {
      bdiamtype[i] = 0.5;
      if (i % 6 == 1) bcolortype[i] = image->color2rgb("red");
      else if (i % 6 == 2) bcolortype[i] = image->color2rgb("green");
      else if (i % 6 == 3) bcolortype[i] = image->color2rgb("blue");
      else if (i % 6 == 4) bcolortype[i] = image->color2rgb("yellow");
      else if (i % 6 == 5) bcolortype[i] = image->color2rgb("aqua");
      else if (i % 6 == 0) bcolortype[i] = image->color2rgb("cyan");
    }
  } else {
    bdiamtype = nullptr;
    bcolortype = nullptr;
  }

  // viewflag = DYNAMIC if any view parameter is dynamic

  viewflag = STATIC;
  if (thetastr || phistr || cflag == DYNAMIC ||
      upxstr || upystr || upzstr || zoomstr) viewflag = DYNAMIC;

  box_bounds();
  if (cflag == STATIC) box_center();
  if (viewflag == STATIC) view_params();

  // local data

  grid_compute = nullptr;
  grid_fix = nullptr;

  maxbufcopy = 0;
  chooseghost = nullptr;
  bufcopy = nullptr;

  maxgrid = 0;
  gbuf = nullptr;
}

/* ---------------------------------------------------------------------- */

DumpImage::~DumpImage()
{
  delete image;

  delete [] diamtype;
  delete [] diamelement;
  delete [] colortype;
  delete [] colorelement;
  delete [] bdiamtype;
  delete [] bcolortype;
  memory->destroy(chooseghost);
  memory->destroy(bufcopy);
  memory->destroy(gbuf);

  delete [] id_grid_compute;
  delete [] id_grid_fix;
}

/* ---------------------------------------------------------------------- */

void DumpImage::init_style()
{
  if (multifile == 0 && !multifile_override)
    error->all(FLERR,"Dump image requires one snapshot per file");
  if (sort_flag) error->all(FLERR,"Dump image cannot perform sorting");

  DumpCustom::init_style();

  // for grid output, find current ptr for compute or fix
  // check that fix frequency is acceptable

  if (gridflag) {
    if (id_grid_compute) {
      grid_compute = modify->get_compute_by_id(id_grid_compute);
      if (!grid_compute)
        error->all(FLERR,"Could not find dump image grid compute ID {}",id_grid_compute);
    } else if (id_grid_fix) {
      grid_fix = modify->get_fix_by_id(id_grid_fix);
      if (!grid_fix) error->all(FLERR,"Could not find dump image fix ID {}",id_grid_fix);
      if (nevery % grid_fix->peratom_freq)
        error->all(FLERR,"Dump image and grid fix not computed at compatible times");
    }
  }

  // check image variables

  if (thetastr) {
    thetavar = input->variable->find(thetastr);
    if (thetavar < 0)
      error->all(FLERR,"Variable name for dump image theta does not exist");
    if (!input->variable->equalstyle(thetavar))
      error->all(FLERR,"Variable for dump image theta is invalid style");
  }
  if (phistr) {
    phivar = input->variable->find(phistr);
    if (phivar < 0)
      error->all(FLERR,"Variable name for dump image phi does not exist");
    if (!input->variable->equalstyle(phivar))
      error->all(FLERR,"Variable for dump image phi is invalid style");
  }
  if (cxstr) {
    cxvar = input->variable->find(cxstr);
    if (cxvar < 0)
      error->all(FLERR,"Variable name for dump image center does not exist");
    if (!input->variable->equalstyle(cxvar))
      error->all(FLERR,"Variable for dump image center is invalid style");
  }
  if (cystr) {
    cyvar = input->variable->find(cystr);
    if (cyvar < 0)
      error->all(FLERR,"Variable name for dump image center does not exist");
    if (!input->variable->equalstyle(cyvar))
      error->all(FLERR,"Variable for dump image center is invalid style");
  }
  if (czstr) {
    czvar = input->variable->find(czstr);
    if (czvar < 0)
      error->all(FLERR,"Variable name for dump image center does not exist");
    if (!input->variable->equalstyle(czvar))
      error->all(FLERR,"Variable for dump image center is invalid style");
  }
  if (upxstr) {
    upxvar = input->variable->find(upxstr);
    if (upxvar < 0)
      error->all(FLERR,"Variable name for dump image center does not exist");
    if (!input->variable->equalstyle(upxvar))
      error->all(FLERR,"Variable for dump image center is invalid style");
  }
  if (upystr) {
    upyvar = input->variable->find(upystr);
    if (upyvar < 0)
      error->all(FLERR,"Variable name for dump image center does not exist");
    if (!input->variable->equalstyle(upyvar))
      error->all(FLERR,"Variable for dump image center is invalid style");
  }
  if (upzstr) {
    upzvar = input->variable->find(upzstr);
    if (upzvar < 0)
      error->all(FLERR,"Variable name for dump image center does not exist");
    if (!input->variable->equalstyle(upzvar))
      error->all(FLERR,"Variable for dump image center is invalid style");
  }
  if (zoomstr) {
    zoomvar = input->variable->find(zoomstr);
    if (zoomvar < 0)
      error->all(FLERR,"Variable name for dump image zoom does not exist");
    if (!input->variable->equalstyle(zoomvar))
      error->all(FLERR,"Variable for dump image zoom is invalid style");
  }

  // set up type -> element mapping

  if (atomflag && acolor == ELEMENT) {
    for (int i = 1; i <= ntypes; i++) {
      colorelement[i] = image->element2color(typenames[i]);
      if (colorelement[i] == nullptr)
        error->all(FLERR,"Invalid dump image element name");
    }
  }

  if (atomflag && adiam == ELEMENT) {
    for (int i = 1; i <= ntypes; i++) {
      diamelement[i] = image->element2diam(typenames[i]);
      if (diamelement[i] == 0.0)
        error->all(FLERR,"Invalid dump image element name");
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
    // only if within a run or minimize
    // else require the compute is current
    // this prevents the compute from being invoked by the WriteDump class

    if (grid_compute) {
      if (update->whichflag == 0) {
        if (grid_compute->invoked_pergrid != update->ntimestep)
          error->all(FLERR,"Grid compute {} used in dump image between runs is not current",
                     grid_compute->id);
      } else {
        if (!(grid_compute->invoked_flag & Compute::INVOKED_PERGRID)) {
          grid_compute->compute_pergrid();
          grid_compute->invoked_flag |= Compute::INVOKED_PERGRID;
        }
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
    else image->write_PPM(fp);
    if (multifile) {
      fclose(fp);
      fp = nullptr;
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
  double diameter,delx,dely,delz;
  int *bodyvec,*fixvec;
  double **bodyarray,**fixarray;
  double *color,*color1,*color2;
  double *p1,*p2,*p3;
  double xmid[3],pt1[3],pt2[3],pt3[3];
  double mat[3][3];

  // render my atoms

  if (atomflag) {
    double **x = atom->x;
    int *line = atom->line;
    int *tri = atom->tri;
    int *body = atom->body;

    m = 0;
    for (i = 0; i < nchoose; i++) {
      j = clist[i];

      if (acolor == TYPE) {
        itype = static_cast<int> (buf[m]);
        color = colortype[itype];
      } else if (acolor == ELEMENT) {
        itype = static_cast<int> (buf[m]);
        color = colorelement[itype];
      } else if (acolor == ATTRIBUTE) {
        color = image->map_value2color(0,buf[m]);
      } else color = image->color2rgb("white");

      if (adiam == NUMERIC) {
        diameter = adiamvalue;
      } else if (adiam == TYPE) {
        itype = static_cast<int> (buf[m+1]);
        diameter = diamtype[itype];
      } else if (adiam == ELEMENT) {
        itype = static_cast<int> (buf[m+1]);
        diameter = diamelement[itype];
      } else if (adiam == ATTRIBUTE) {
        diameter = buf[m+1];
      }

      // do not draw if line,tri,body keywords enabled and atom is one of those

      drawflag = 1;
      if (extraflag) {
        if (lineflag && line[j] >= 0) drawflag = 0;
        if (triflag && tri[j] >= 0) drawflag = 0;
        if (bodyflag && body[j] >= 0) drawflag = 0;
      }

      if (drawflag) image->draw_sphere(x[j],color,diameter);

      m += size_one;
    }
  }

  // render my grid cells
  // 2 triangles for 2d rectangle, 12 triangles for 3d cube surface
  // grid_cell_corners_2d/3d calculates orthogonal vs triclinic corner pts
  // for 3d, outward normals on all 6 faces

  if (gridflag) {

    // reset lighting for flat surfaces to make them brighter

    image->ambientColor[0] = image->ambientColor[1] = image->ambientColor[2] = 0.9;
    image->keyLightColor[0] = image->keyLightColor[1] = image->keyLightColor[2] = 0.3;
    image->fillLightColor[0] = image->fillLightColor[1] = image->fillLightColor[2] = 0.3;
    image->backLightColor[0] = image->backLightColor[1] = image->backLightColor[2] = 0.3;

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

    image->ambientColor[0] = image->ambientColor[1] = image->ambientColor[2] = 0.0;
    image->keyLightColor[0] = image->keyLightColor[1] = image->keyLightColor[2] = 0.9;
    image->fillLightColor[0] = image->fillLightColor[1] = image->fillLightColor[2] = 0.45;
    image->backLightColor[0] = image->backLightColor[1] = image->backLightColor[2] = 0.9;
  }

  // render atoms that are lines

  if (lineflag) {
    double length,theta,dx,dy;
    double **x = atom->x;
    int *line = atom->line;
    int *type = atom->type;

    for (i = 0; i < nchoose; i++) {
      j = clist[i];
      if (line[j] < 0) continue;

      if (lcolor == TYPE) {
        color = colortype[type[j]];
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

      image->draw_cylinder(pt1,pt2,color,ldiamvalue,3);
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

    for (i = 0; i < nchoose; i++) {
      j = clist[i];
      if (tri[j] < 0) continue;

      if (tcolor == TYPE) {
        color = colortype[type[j]];
      }

      MathExtra::quat_to_mat(avec_tri->bonus[tri[i]].quat,mat);
      MathExtra::matvec(mat,avec_tri->bonus[tri[i]].c1,pt1);
      MathExtra::matvec(mat,avec_tri->bonus[tri[i]].c2,pt2);
      MathExtra::matvec(mat,avec_tri->bonus[tri[i]].c3,pt3);
      MathExtra::add3(pt1,x[i],pt1);
      MathExtra::add3(pt2,x[i],pt2);
      MathExtra::add3(pt3,x[i],pt3);

      if (tridraw) image->draw_triangle(pt1,pt2,pt3,color);
      if (edgedraw) {
        image->draw_cylinder(pt1,pt2,color,tdiamvalue,3);
        image->draw_cylinder(pt2,pt3,color,tdiamvalue,3);
        image->draw_cylinder(pt3,pt1,color,tdiamvalue,3);
      }
    }
  }

  // render atoms that are bodies

  if (bodyflag) {
    Body *bptr = avec_body->bptr;
    int *body = atom->body;

    m = 0;
    for (i = 0; i < nchoose; i++) {
      j = clist[i];
      if (body[j] < 0) continue;

      if (bodycolor == TYPE) {
        itype = static_cast<int> (buf[m]);
        color = colortype[itype];
      }

      ibonus = body[i];
      n = bptr->image(ibonus,bodyflag1,bodyflag2,bodyvec,bodyarray);
      for (k = 0; k < n; k++) {
        if (bodyvec[k] == SPHERE)
          image->draw_sphere(bodyarray[k],color,bodyarray[k][3]);
        else if (bodyvec[k] == LINE)
          image->draw_cylinder(&bodyarray[k][0],&bodyarray[k][3],
                               color,bodyarray[k][6],3);
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

  if (bondflag) {
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
          itype = btype;
          if (itype < 0) itype = -itype;
          color = bcolortype[itype];
        }

        if (bdiam == NUMERIC) {
          diameter = bdiamvalue;
        } else if (bdiam == ATOM) {
          if (adiam == NUMERIC) {
            diameter = adiamvalue;
          } else if (adiam == TYPE) {
            diameter = MIN(diamtype[type[atom1]],diamtype[type[atom1]]);
          } else if (adiam == ELEMENT) {
            diameter = MIN(diamelement[type[atom1]],diamelement[type[atom1]]);
          } else if (adiam == ATTRIBUTE) {
            diameter = MIN(bufcopy[atom1][1],bufcopy[atom2][1]);
          }
        } else if (bdiam == TYPE) {
          itype = btype;
          if (itype < 0) itype = -itype;
          diameter = bdiamtype[itype];
        }

        // draw cylinder in 2 pieces if bcolor = ATOM
        // or bond crosses periodic boundary

        delx = x[atom2][0] - x[atom1][0];
        dely = x[atom2][1] - x[atom1][1];
        delz = x[atom2][2] - x[atom1][2];

        if (bcolor == ATOM || domain->minimum_image_check(delx,dely,delz)) {
          domain->minimum_image(delx,dely,delz);
          xmid[0] = x[atom1][0] + 0.5*delx;
          xmid[1] = x[atom1][1] + 0.5*dely;
          xmid[2] = x[atom1][2] + 0.5*delz;
          if (bcolor == ATOM)
            image->draw_cylinder(x[atom1],xmid,color1,diameter,3);
          else image->draw_cylinder(x[atom1],xmid,color,diameter,3);
          xmid[0] = x[atom2][0] - 0.5*delx;
          xmid[1] = x[atom2][1] - 0.5*dely;
          xmid[2] = x[atom2][2] - 0.5*delz;
          if (bcolor == ATOM)
            image->draw_cylinder(xmid,x[atom2],color2,diameter,3);
          else image->draw_cylinder(xmid,x[atom2],color,diameter,3);

        } else image->draw_cylinder(x[atom1],x[atom2],color,diameter,3);
      }
    }
  }

  // render objects provided by a fix

  if (fixflag) {
    int tridraw=0,edgedraw=0;
    if (domain->dimension == 3) {
      tridraw = 1;
      edgedraw = 1;
      if ((int) fixflag1 == 2) tridraw = 0;
      if ((int) fixflag1 == 1) edgedraw = 0;
    }

    n = fixptr->image(fixvec,fixarray);

    for (i = 0; i < n; i++) {
      if (fixvec[i] == SPHERE) {
        // no fix draws spheres yet
      } else if (fixvec[i] == LINE) {
        if (fixcolor == TYPE) {
          itype = static_cast<int> (fixarray[i][0]);
          color = colortype[itype];
        }
        image->draw_cylinder(&fixarray[i][1],&fixarray[i][4],
                             color,fixflag1,3);
      } else if (fixvec[i] == TRI) {
        if (fixcolor == TYPE) {
          itype = static_cast<int> (fixarray[i][0]);
          color = colortype[itype];
        }
        p1 = &fixarray[i][1];
        p2 = &fixarray[i][4];
        p3 = &fixarray[i][7];
        if (tridraw) image->draw_triangle(p1,p2,p3,color);
        if (edgedraw) {
          image->draw_cylinder(p1,p2,color,fixflag2,3);
          image->draw_cylinder(p2,p3,color,fixflag2,3);
          image->draw_cylinder(p3,p1,color,fixflag2,3);
        }
      }
    }
  }

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

    image->draw_box(boxcorners,diameter);
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

    image->draw_box(boxcorners,diameter);
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
    offset *= 0.1;
    axes[0][0] -= offset; axes[0][1] -= offset; axes[0][2] -= offset;
    axes[1][0] -= offset; axes[1][1] -= offset; axes[1][2] -= offset;
    axes[2][0] -= offset; axes[2][1] -= offset; axes[2][2] -= offset;
    axes[3][0] -= offset; axes[3][1] -= offset; axes[3][2] -= offset;

    axes[1][0] = axes[0][0] + axeslen*(axes[1][0]-axes[0][0]);
    axes[1][1] = axes[0][1] + axeslen*(axes[1][1]-axes[0][1]);
    axes[1][2] = axes[0][2] + axeslen*(axes[1][2]-axes[0][2]);
    axes[2][0] = axes[0][0] + axeslen*(axes[2][0]-axes[0][0]);
    axes[2][1] = axes[0][1] + axeslen*(axes[2][1]-axes[0][1]);
    axes[2][2] = axes[0][2] + axeslen*(axes[2][2]-axes[0][2]);
    axes[3][0] = axes[0][0] + axeslen*(axes[3][0]-axes[0][0]);
    axes[3][1] = axes[0][1] + axeslen*(axes[3][1]-axes[0][1]);
    axes[3][2] = axes[0][2] + axeslen*(axes[3][2]-axes[0][2]);

    image->draw_axes(axes,diameter);
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
    for (i = first; i < last; i++) chooseghost[i] = static_cast<int> (buf[m++]);
  else {
    for (i = first; i < last; i++) {
      chooseghost[i] = static_cast<int> (buf[m++]);
      bufcopy[i][0] = buf[m++];
      bufcopy[i][1] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

int DumpImage::modify_param(int narg, char **arg)
{
  int n = DumpCustom::modify_param(narg,arg);
  if (n) return n;

  if (strcmp(arg[0],"acolor") == 0) {
    if (narg < 3) error->all(FLERR,"Illegal dump_modify command");
    int nlo,nhi;
    utils::bounds(FLERR,arg[1],1,atom->ntypes,nlo,nhi,error);

    // get list of colors
    // assign colors in round-robin fashion to types

    auto colors = Tokenizer(arg[2],"/").as_vector();
    const int ncolors = colors.size();

    int m = 0;
    for (int i = nlo; i <= nhi; i++) {
      colortype[i] = image->color2rgb(colors[m%ncolors].c_str());
      if (colortype[i] == nullptr)
        error->all(FLERR,"Invalid color in dump_modify command");
      m++;
    }
    return 3;
  }

  if (strcmp(arg[0],"adiam") == 0) {
    if (narg < 3) error->all(FLERR,"Illegal dump_modify command");
    int nlo,nhi;
    utils::bounds(FLERR,arg[1],1,atom->ntypes,nlo,nhi,error);
    double diam = utils::numeric(FLERR,arg[2],false,lmp);
    if (diam <= 0.0) error->all(FLERR,"Illegal dump_modify command");
    for (int i = nlo; i <= nhi; i++) diamtype[i] = diam;
    return 3;
  }

  if ((strcmp(arg[0],"amap") == 0) || (strcmp(arg[0],"gmap") == 0)) {
    if (narg < 6) error->all(FLERR,"Illegal dump_modify command");
    if (strlen(arg[3]) != 2) error->all(FLERR,"Illegal dump_modify command");
    int factor = 0;
    if (arg[3][0] == 's') factor = 1;
    else if (arg[3][0] == 'c') factor = 2;
    else if (arg[3][0] == 'd') factor = 3;
    else error->all(FLERR,"Illegal dump_modify command");
    int nentry = utils::inumeric(FLERR,arg[5],false,lmp);
    if (nentry < 1) error->all(FLERR,"Illegal dump_modify command");
    n = 6 + factor*nentry;
    if (narg < n) error->all(FLERR,"Illegal dump_modify command");
    int flag;
    if (strcmp(arg[0],"amap") == 0) flag = image->map_reset(0,n-1,&arg[1]);
    if (strcmp(arg[0],"gmap") == 0) flag = image->map_reset(1,n-1,&arg[1]);
    if (flag) error->all(FLERR,"Illegal dump_modify command");
    return n;
  }

  if (strcmp(arg[0],"bcolor") == 0) {
    if (narg < 3) error->all(FLERR,"Illegal dump_modify command");
    if (atom->nbondtypes == 0)
      error->all(FLERR,"Dump modify bcolor not allowed with no bond types");
    int nlo,nhi;
    utils::bounds(FLERR,arg[1],1,atom->nbondtypes,nlo,nhi,error);

    // process list of ncount colornames separated by '/'
    // assign colors in round-robin fashion to bond types

    auto colors = Tokenizer(arg[2],"/").as_vector();
    const int ncolors = colors.size();

    int m = 0;
    for (int i = nlo; i <= nhi; i++) {
      bcolortype[i] = image->color2rgb(colors[m%ncolors].c_str());
      if (bcolortype[i] == nullptr)
        error->all(FLERR,"Invalid color in dump_modify command");
      m++;
    }
    return 3;
  }

  if (strcmp(arg[0],"bdiam") == 0) {
    if (narg < 3) error->all(FLERR,"Illegal dump_modify command");
    if (atom->nbondtypes == 0)
      error->all(FLERR,"Dump modify bdiam not allowed with no bond types");
    int nlo,nhi;
    utils::bounds(FLERR,arg[1],1,atom->nbondtypes,nlo,nhi,error);
    double diam = utils::numeric(FLERR,arg[2],false,lmp);
    if (diam <= 0.0) error->all(FLERR,"Illegal dump_modify command");
    for (int i = nlo; i <= nhi; i++) bdiamtype[i] = diam;
    return 3;
  }

  if (strcmp(arg[0],"backcolor") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    double *color = image->color2rgb(arg[1]);
    if (color == nullptr) error->all(FLERR,"Invalid color in dump_modify command");
    image->background[0] = static_cast<int> (color[0]*255.0);
    image->background[1] = static_cast<int> (color[1]*255.0);
    image->background[2] = static_cast<int> (color[2]*255.0);
    return 2;
  }

  if (strcmp(arg[0],"boxcolor") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    image->boxcolor = image->color2rgb(arg[1]);
    if (image->boxcolor == nullptr)
      error->all(FLERR,"Invalid color in dump_modify command");
    return 2;
  }

  if (strcmp(arg[0],"color") == 0) {
    if (narg < 5) error->all(FLERR,"Illegal dump_modify command");
    int flag = image->addcolor(arg[1],utils::numeric(FLERR,arg[2],false,lmp),
                               utils::numeric(FLERR,arg[3],false,lmp),
                               utils::numeric(FLERR,arg[4],false,lmp));
    if (flag) error->all(FLERR,"Illegal dump_modify command");
    return 5;
  }

  return 0;
}
