/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "ctype.h"
#include "stdlib.h"
#include "string.h"
#include "dump_image.h"
#include "image.h"
#include "atom.h"
#include "atom_vec.h"
#include "molecule.h"
#include "domain.h"
#include "group.h"
#include "force.h"
#include "comm.h"
#include "input.h"
#include "variable.h"
#include "math_const.h"
#include "error.h"
#include "memory.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define BIG 1.0e20

enum{NUMERIC,ATOM,TYPE,ELEMENT,ATTRIBUTE};
enum{STATIC,DYNAMIC};
enum{NO,YES};

/* ---------------------------------------------------------------------- */

DumpImage::DumpImage(LAMMPS *lmp, int narg, char **arg) :
  DumpCustom(lmp, narg, arg)
{
  if (binary || multiproc) error->all(FLERR,"Invalid dump image filename");

  // force binary flag on to avoid corrupted output on Windows

  binary = 1;
  multifile_override = 0;

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

  // create Image class with single colormap for atoms
  // change defaults for 2d

  image = new Image(lmp,1);

  if (domain->dimension == 2) {
    image->theta = 0.0;
    image->phi = 0.0;
    image->up[0] = 0.0; image->up[1] = 1.0; image->up[2] = 0.0;
  }

  // set defaults for optional args

  atomflag = YES;
  if (atom->nbondtypes == 0) bondflag = NO;
  else {
    bondflag = YES;
    bcolor = ATOM;
    bdiam = NUMERIC;
    bdiamvalue = 0.5;
  }

  thetastr = phistr = NULL;
  cflag = STATIC;
  cx = cy = cz = 0.5;
  cxstr = cystr = czstr = NULL;

  upxstr = upystr = upzstr = NULL;
  zoomstr = NULL;
  perspstr = NULL;
  boxflag = YES;
  boxdiam = 0.02;
  axesflag = NO;

  // parse optional args

  int iarg = ioptional;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"adiam") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump image command");
      adiam = NUMERIC;
      adiamvalue = force->numeric(FLERR,arg[iarg+1]);
      if (adiamvalue <= 0.0) error->all(FLERR,"Illegal dump image command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"atom") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump image command");
      if (strcmp(arg[iarg+1],"yes") == 0) atomflag = YES;
      else if (strcmp(arg[iarg+1],"no") == 0) atomflag = NO;
      else error->all(FLERR,"Illegal dump image command");
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
          bdiamvalue = force->numeric(FLERR,arg[iarg+2]);
          if (bdiamvalue <= 0.0) error->all(FLERR,"Illegal dump image command");
      } else if (strcmp(arg[iarg+2],"atom") == 0) bdiam = ATOM;
      else if (strcmp(arg[iarg+2],"type") == 0) bdiam = TYPE;
      else if (strcmp(arg[iarg+2],"none") == 0) bondflag = NO;
      else error->all(FLERR,"Illegal dump image command");
      iarg += 3;

    } else if (strcmp(arg[iarg],"size") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal dump image command");
      int width = force->inumeric(FLERR,arg[iarg+1]);
      int height = force->inumeric(FLERR,arg[iarg+2]);
      if (width <= 0 || height <= 0)
        error->all(FLERR,"Illegal dump image command");
      image->width = width;
      image->height = height;
      iarg += 3;

    } else if (strcmp(arg[iarg],"view") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal dump image command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
        int n = strlen(&arg[iarg+1][2]) + 1;
        thetastr = new char[n];
        strcpy(thetastr,&arg[iarg+1][2]);
      } else {
        double theta = force->numeric(FLERR,arg[iarg+1]);
        if (theta < 0.0 || theta > 180.0)
          error->all(FLERR,"Invalid dump image theta value");
        theta *= MY_PI/180.0;
        image->theta = theta;
      }
      if (strstr(arg[iarg+2],"v_") == arg[iarg+2]) {
        int n = strlen(&arg[iarg+2][2]) + 1;
        phistr = new char[n];
        strcpy(phistr,&arg[iarg+2][2]);
      } else {
        double phi = force->numeric(FLERR,arg[iarg+2]);
        phi *= MY_PI/180.0;
        image->phi = phi;
      }
      iarg += 3;

    } else if (strcmp(arg[iarg],"center") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal dump image command");
      if (strcmp(arg[iarg+1],"s") == 0) cflag = STATIC;
      else if (strcmp(arg[iarg+1],"d") == 0) cflag = DYNAMIC;
      else error->all(FLERR,"Illegal dump image command");
      if (strstr(arg[iarg+2],"v_") == arg[iarg+2]) {
        int n = strlen(&arg[iarg+2][2]) + 1;
        cxstr = new char[n];
        strcpy(cxstr,&arg[iarg+2][2]);
        cflag = DYNAMIC;
      } else cx = force->numeric(FLERR,arg[iarg+2]);
      if (strstr(arg[iarg+3],"v_") == arg[iarg+3]) {
        int n = strlen(&arg[iarg+3][2]) + 1;
        cystr = new char[n];
        strcpy(cystr,&arg[iarg+3][2]);
        cflag = DYNAMIC;
      } else cy = force->numeric(FLERR,arg[iarg+3]);
      if (strstr(arg[iarg+4],"v_") == arg[iarg+4]) {
        int n = strlen(&arg[iarg+4][2]) + 1;
        czstr = new char[n];
        strcpy(czstr,&arg[iarg+4][2]);
        cflag = DYNAMIC;
      } else cz = force->numeric(FLERR,arg[iarg+4]);
      iarg += 5;

    } else if (strcmp(arg[iarg],"up") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal dump image command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
        int n = strlen(&arg[iarg+1][2]) + 1;
        upxstr = new char[n];
        strcpy(upxstr,&arg[iarg+1][2]);
      } else image->up[0] = force->numeric(FLERR,arg[iarg+1]);
      if (strstr(arg[iarg+2],"v_") == arg[iarg+2]) {
        int n = strlen(&arg[iarg+2][2]) + 1;
        upystr = new char[n];
        strcpy(upystr,&arg[iarg+2][2]);
      } else image->up[1] = force->numeric(FLERR,arg[iarg+2]);
      if (strstr(arg[iarg+3],"v_") == arg[iarg+3]) {
        int n = strlen(&arg[iarg+3][2]) + 1;
        upzstr = new char[n];
        strcpy(upzstr,&arg[iarg+3][2]);
      } else image->up[2] = force->numeric(FLERR,arg[iarg+3]);
      iarg += 4;

    } else if (strcmp(arg[iarg],"zoom") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump image command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
        int n = strlen(&arg[iarg+1][2]) + 1;
        zoomstr = new char[n];
        strcpy(zoomstr,&arg[iarg+1][2]);
      } else {
        double zoom = force->numeric(FLERR,arg[iarg+1]);
        if (zoom <= 0.0) error->all(FLERR,"Illegal dump image command");
        image->zoom = zoom;
      }
      iarg += 2;

    } else if (strcmp(arg[iarg],"persp") == 0) {
      error->all(FLERR,"Dump image persp option is not yet supported");
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump image command");
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) {
        int n = strlen(&arg[iarg+1][2]) + 1;
        perspstr = new char[n];
        strcpy(perspstr,&arg[iarg+1][2]);
      } else {
        double persp = force->numeric(FLERR,arg[iarg+1]);
        if (persp < 0.0) error->all(FLERR,"Illegal dump image command");
        image->persp = persp;
      }
      iarg += 2;

    } else if (strcmp(arg[iarg],"box") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal dump image command");
      if (strcmp(arg[iarg+1],"yes") == 0) boxflag = YES;
      else if (strcmp(arg[iarg+1],"no") == 0) boxflag = NO;
      else error->all(FLERR,"Illegal dump image command");
      boxdiam = force->numeric(FLERR,arg[iarg+2]);
      if (boxdiam < 0.0) error->all(FLERR,"Illegal dump image command");
      iarg += 3;

    } else if (strcmp(arg[iarg],"axes") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal dump image command");
      if (strcmp(arg[iarg+1],"yes") == 0) axesflag = YES;
      else if (strcmp(arg[iarg+1],"no") == 0) axesflag = NO;
      else error->all(FLERR,"Illegal dump image command");
      axeslen = force->numeric(FLERR,arg[iarg+2]);
      axesdiam = force->numeric(FLERR,arg[iarg+3]);
      if (axeslen < 0.0 || axesdiam < 0.0)
        error->all(FLERR,"Illegal dump image command");
      iarg += 4;

    } else if (strcmp(arg[iarg],"shiny") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal dump image command");
      double shiny = force->numeric(FLERR,arg[iarg+1]);
      if (shiny < 0.0 || shiny > 1.0)
        error->all(FLERR,"Illegal dump image command");
      image->shiny = shiny;
      iarg += 2;

    } else if (strcmp(arg[iarg],"ssao") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal dump image command");
      if (strcmp(arg[iarg+1],"yes") == 0) image->ssao = YES;
      else if (strcmp(arg[iarg+1],"no") == 0) image->ssao = NO;
      else error->all(FLERR,"Illegal dump image command");
      int seed = force->inumeric(FLERR,arg[iarg+2]);
      if (seed <= 0) error->all(FLERR,"Illegal dump image command");
      image->seed = seed;
      double ssaoint = force->numeric(FLERR,arg[iarg+3]);
      if (ssaoint < 0.0 || ssaoint > 1.0)
        error->all(FLERR,"Illegal dump image command");
      image->ssaoint = ssaoint;
      iarg += 4;

    } else error->all(FLERR,"Illegal dump image command");
  }

  // allocate image buffer now that image size is known

  image->buffers();

  // communication neede for bonds colored by atoms

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
    bdiamtype = NULL;
    bcolortype = NULL;
  }

  // viewflag = DYNAMIC if any view parameter is dynamic

  viewflag = STATIC;
  if (thetastr || phistr || cflag == DYNAMIC ||
      upxstr || upystr || upzstr || zoomstr || perspstr) viewflag = DYNAMIC;

  box_bounds();
  if (cflag == STATIC) box_center();
  if (viewflag == STATIC) view_params();

  // local data

  maxbufcopy = 0;
  chooseghost = NULL;
  bufcopy = NULL;
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
}

/* ---------------------------------------------------------------------- */

void DumpImage::init_style()
{
  if (multifile == 0 && !multifile_override)
    error->all(FLERR,"Dump image requires one snapshot per file");
  if (sort_flag) error->all(FLERR,"Dump image cannot perform sorting");

  DumpCustom::init_style();

  // check variables

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
  if (perspstr) {
    perspvar = input->variable->find(perspstr);
    if (perspvar < 0)
      error->all(FLERR,"Variable name for dump image persp does not exist");
    if (!input->variable->equalstyle(perspvar))
      error->all(FLERR,"Variable for dump image persp is invalid style");
  }

  // set up type -> element mapping

  if (atomflag && acolor == ELEMENT) {
    for (int i = 1; i <= ntypes; i++) {
      colorelement[i] = image->element2color(typenames[i]);
      if (colorelement[i] == NULL)
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

  // pack buf with color & diameter

  pack(NULL);

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
    if (flag) error->all(FLERR,"Invalid color map min/max values");
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
      fp = NULL;
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
    double theta = input->variable->compute_equal(thetavar);
    if (theta < 0.0 || theta > 180.0)
      error->all(FLERR,"Invalid dump image theta value");
    theta *= MY_PI/180.0;
    image->theta = theta;
  }

  if (phistr) {
    double phi = input->variable->compute_equal(phivar);
    phi *= MY_PI/180.0;
    image->phi = phi;
  }

  // up vector

  if (upxstr) image->up[0] = input->variable->compute_equal(upxvar);
  if (upystr) image->up[1] = input->variable->compute_equal(upyvar);
  if (upzstr) image->up[2] = input->variable->compute_equal(upzvar);

  // zoom and perspective

  if (zoomstr) image->zoom = input->variable->compute_equal(zoomvar);
  if (image->zoom <= 0.0) error->all(FLERR,"Invalid dump image zoom value");
  if (perspstr) image->persp = input->variable->compute_equal(perspvar);
  if (image->persp < 0.0) error->all(FLERR,"Invalid dump image persp value");

  // remainder of view setup is internal to Image class

  image->view_params(boxxlo,boxxhi,boxylo,boxyhi,boxzlo,boxzhi);
}

/* ----------------------------------------------------------------------
   create image for atoms on this proc
   every pixel has depth
------------------------------------------------------------------------- */

void DumpImage::create_image()
{
  int i,j,m,n,itype,atom1,atom2,imol,iatom,btype;
  tagint tagprev;
  double diameter,delx,dely,delz;
  double *color,*color1,*color2;
  double xmid[3];

  // render my atoms

  if (atomflag) {
    double **x = atom->x;

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
      }

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

      image->draw_sphere(x[j],color,diameter);
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
    int nall = atom->nlocal + atom->nghost;
    int newton_bond = force->newton_bond;
    int molecular = atom->molecular;
    Molecule **onemols = atom->avec->onemols;

    // communicate choose flag for ghost atoms to know if they are selected
    // if bcolor/bdiam = ATOM, setup bufcopy to comm atom color/diam attributes

    if (nall > maxbufcopy) {
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

    comm->forward_comm_dump(this);

    for (i = 0; i < nchoose; i++) {
      atom1 = clist[i];
      if (molecular == 1) n = num_bond[atom1];
      else {
        if (molindex[atom1] < 0) continue;
        imol = molindex[atom1];
        iatom = molatom[atom1];
        n = onemols[imol]->num_bond[iatom];
      }

      for (m = 0; m < n; m++) {
        if (molecular == 1) {
          btype = bond_type[atom1][m];
          atom2 = atom->map(bond_atom[atom1][m]);
        } else {
          tagprev = tag[i] - iatom - 1;
          btype = atom->map(onemols[imol]->bond_type[atom1][m]);
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

  // render outline of simulation box, orthogonal or triclinic

  if (boxflag) {
    double diameter = MIN(boxxhi-boxxlo,boxyhi-boxylo);
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
    double diameter = MIN(boxxhi-boxxlo,boxyhi-boxylo);
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

int DumpImage::pack_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
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

  return comm_forward;
}

/* ---------------------------------------------------------------------- */

void DumpImage::unpack_comm(int n, int first, double *buf)
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
    force->bounds(arg[1],atom->ntypes,nlo,nhi);

    // ptrs = list of ncount colornames separated by '/'

    int ncount = 1;
    char *nextptr;
    char *ptr = arg[2];
    while (nextptr = strchr(ptr,'/')) {
      ptr = nextptr + 1;
      ncount++;
    }
    char **ptrs = new char*[ncount+1];
    ncount = 0;
    ptrs[ncount++] = strtok(arg[2],"/");
    while (ptrs[ncount++] = strtok(NULL,"/"));
    ncount--;

    // assign each of ncount colors in round-robin fashion to types

    int m = 0;
    for (int i = nlo; i <= nhi; i++) {
      colortype[i] = image->color2rgb(ptrs[m%ncount]);
      if (colortype[i] == NULL)
        error->all(FLERR,"Invalid color in dump_modify command");
      m++;
    }

    delete [] ptrs;
    return 3;
  }

  if (strcmp(arg[0],"adiam") == 0) {
    if (narg < 3) error->all(FLERR,"Illegal dump_modify command");
    int nlo,nhi;
    force->bounds(arg[1],atom->ntypes,nlo,nhi);
    double diam = force->numeric(FLERR,arg[2]);
    if (diam <= 0.0) error->all(FLERR,"Illegal dump_modify command");
    for (int i = nlo; i <= nhi; i++) diamtype[i] = diam;
    return 3;
  }

  if (strcmp(arg[0],"amap") == 0) {
    if (narg < 6) error->all(FLERR,"Illegal dump_modify command");
    if (strlen(arg[3]) != 2) error->all(FLERR,"Illegal dump_modify command");
    int factor = 2;
    if (arg[3][0] == 's') factor = 1;
    int nentry = force->inumeric(FLERR,arg[5]);
    if (nentry < 1) error->all(FLERR,"Illegal dump_modify command");
    int n = 6 + factor*nentry;
    if (narg < n) error->all(FLERR,"Illegal dump_modify command");
    int flag = image->map_reset(0,n-1,&arg[1]);
    if (flag) error->all(FLERR,"Illegal dump_modify command");
    return n;
  }

  if (strcmp(arg[0],"bcolor") == 0) {
    if (narg < 3) error->all(FLERR,"Illegal dump_modify command");
    if (atom->nbondtypes == 0)
      error->all(FLERR,"Dump modify bcolor not allowed with no bond types");
    int nlo,nhi;
    force->bounds(arg[1],atom->nbondtypes,nlo,nhi);

    // ptrs = list of ncount colornames separated by '/'

    int ncount = 1;
    char *nextptr;
    char *ptr = arg[2];
    while (nextptr = strchr(ptr,'/')) {
      ptr = nextptr + 1;
      ncount++;
    }
    char **ptrs = new char*[ncount+1];
    ncount = 0;
    ptrs[ncount++] = strtok(arg[2],"/");
    while (ptrs[ncount++] = strtok(NULL,"/"));
    ncount--;

    // assign each of ncount colors in round-robin fashion to types

    int m = 0;
    for (int i = nlo; i <= nhi; i++) {
      bcolortype[i] = image->color2rgb(ptrs[m%ncount]);
      if (bcolortype[i] == NULL)
        error->all(FLERR,"Invalid color in dump_modify command");
      m++;
    }

    delete [] ptrs;
    return 3;
  }

  if (strcmp(arg[0],"bdiam") == 0) {
    if (narg < 3) error->all(FLERR,"Illegal dump_modify command");
    if (atom->nbondtypes == 0)
      error->all(FLERR,"Dump modify bdiam not allowed with no bond types");
    int nlo,nhi;
    force->bounds(arg[1],atom->ntypes,nlo,nhi);
    double diam = force->numeric(FLERR,arg[2]);
    if (diam <= 0.0) error->all(FLERR,"Illegal dump_modify command");
    for (int i = nlo; i <= nhi; i++) bdiamtype[i] = diam;
    return 3;
  }

  if (strcmp(arg[0],"backcolor") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    double *color = image->color2rgb(arg[1]);
    if (color == NULL) error->all(FLERR,"Invalid color in dump_modify command");
    image->background[0] = static_cast<int> (color[0]*255.0);
    image->background[1] = static_cast<int> (color[1]*255.0);
    image->background[2] = static_cast<int> (color[2]*255.0);
    return 2;
  }

  if (strcmp(arg[0],"boxcolor") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    image->boxcolor = image->color2rgb(arg[1]);
    if (image->boxcolor == NULL)
      error->all(FLERR,"Invalid color in dump_modify command");
    return 2;
  }

  if (strcmp(arg[0],"color") == 0) {
    if (narg < 5) error->all(FLERR,"Illegal dump_modify command");
    int flag = image->addcolor(arg[1],force->numeric(FLERR,arg[2]),force->numeric(FLERR,arg[3]),force->numeric(FLERR,arg[4]));
    if (flag) error->all(FLERR,"Illegal dump_modify command");
    return 5;
  }

  return 0;
}
