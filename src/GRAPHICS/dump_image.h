/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef DUMP_CLASS
// clang-format off
DumpStyle(image,DumpImage);
// clang-format on
#else

#ifndef LMP_DUMP_IMAGE_H
#define LMP_DUMP_IMAGE_H

#include "dump_custom.h"

namespace LAMMPS_NS {

// forward declarations
class AtomVecBody;
class AtomVecEllipsoid;
class AtomVecLine;
class AtomVecTri;
class Compute;
class Fix;
class Grid2d;
class Grid3d;
class Image;
class LAMMPS;
class Region;

class DumpImage : public DumpCustom {
 public:
  DumpImage(LAMMPS *, int, char **);
  ~DumpImage() override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

  void *extract(const char *, int &) override;

 protected:
  int filetype;
  enum { PPM, JPG, PNG, TGA };    // file type constants

  int atomflag;         // 0/1 for draw atoms
  int acolor, adiam;    // what determines color/diam of atoms
  double adiamvalue;    // atom diameter value

  int lineflag;                   // 0/1 for draw atoms as lines
  int lcolor, ldiam;              // what determines color/diam of lines
  double ldiamvalue;              // line diameter value
  int triflag;                    // 0/1 for draw atoms as triangles
  int tcolor, tstyle;             // what determines color/style of tris
  double tdiamvalue;              // tri edge diameter value
  int ellipsoidflag;              // 0/1 for draw atoms as ellipsoids
  int ecolor, estyle;             // what determines color/style of ellipsoid
  int elevel;                     // mesh refinement level 1, 2, 3, or 4
  double ediamvalue;              // ellipsoid edge diameter value
  int bodyflag;                   // 0/1 for draw atoms as bodies
  int bodycolor;                  // what determines color of bodies
  double bodyflag1, bodyflag2;    // user-specified params for drawing bodies

  int bondflag;         // NO/YES/AUTO for drawing bonds
  int bcolor, bdiam;    // what determines color/diam of bonds
  double bdiamvalue;    // bond diameter value
  double bondcutoff;    // autobond cutoff

  int extraflag;                                    // 0/1 for any of line/tri/body flag set
  char *thetastr, *phistr;                          // variables for view theta,phi
  int thetavar, phivar;                             // index to theta,phi vars
  int cflag;                                        // static/dynamic box center
  double cx, cy, cz;                                // fractional box center
  char *cxstr, *cystr, *czstr;                      // variables for box center
  int cxvar, cyvar, czvar;                          // index to box center vars
  char *upxstr, *upystr, *upzstr;                   // view up vector variables
  int upxvar, upyvar, upzvar;                       // index to up vector vars
  char *zoomstr;                                    // view zoom variable name
  int zoomvar;                                      // index to zoom variable
  int boxflag, axesflag;                            // 0/1 for draw box and axes
  double boxdiam, axeslen, axesdiam;                // params for drawing box and axes
  double boxopacity, axesopacity, subboxopacity;    // opacity for box, subbox, axes
  int subboxflag;
  double subboxdiam;

  int viewflag;    // overall view is static or dynamic

  double *diamtype, *diamelement, *bdiamtype;          // per-type diameters
  double **colortype, **colorelement, **bcolortype;    // per-type colors
  double *aopacity, *bopacity;                         // per-type opacity

  int gridflag;    // 0/1 for draw grid cells
  Grid2d *grid2d;
  Grid3d *grid3d;
  char *id_grid_compute, *id_grid_fix;
  Compute *grid_compute;
  Fix *grid_fix;
  int grid_igrid, grid_idata, grid_index;
  int nxgrid, nygrid, nzgrid;
  int nxlo_in, nxhi_in, nylo_in, nyhi_in, nzlo_in, nzhi_in;
  double *gbuf;
  int ngrid, maxgrid;
  double gcorners[8][3];

  AtomVecLine *avec_line;    // ptrs to atom style (sub)classes
  AtomVecTri *avec_tri;
  AtomVecEllipsoid *avec_ellipsoid;
  AtomVecBody *avec_body;

  struct ObjInfo {
    ObjInfo() = delete;
    ObjInfo(const std::string &_id, Compute *_cptr, Fix *_fptr, int _colorstyle, double _flag1,
            double _flag2, double *_rgb, double _opacity = 1.0) :
        id(_id), cptr(_cptr), fptr(_fptr), colorstyle(_colorstyle), flag1(_flag1), flag2(_flag2),
        rgb(_rgb), opacity(_opacity)
    {
    }

    std::string id;
    Compute *cptr;
    Fix *fptr;
    int colorstyle;
    double flag1;
    double flag2;
    double *rgb;
    double opacity;
  };
  std::vector<ObjInfo> objects;

  Image *image;    // class that renders each image

  struct RegionInfo {
    RegionInfo() = delete;
    RegionInfo(const std::string &_id, Region *_ptr, double *_color, int _style,
               double _diameter = 0.5, double _opacity = 1.0, int _npoints = 0) :
        ptr(_ptr), id(_id), style(_style), color(_color), diameter(_diameter), opacity(_opacity),
        npoints(_npoints)
    {
    }

    Region *ptr;
    std::string id;
    int style;
    double *color;
    double diameter;
    double opacity;
    int npoints;
  };

  std::vector<RegionInfo> regions;

  int *chooseghost;    // extended choose array for comm
  double **bufcopy;    // buffer for communicating bond/atom info
  int maxbufcopy;

  void init_style() override;
  int modify_param(int, char **) override;
  void write() override;

  void box_center();
  void view_params();
  void box_bounds();

  void create_image();
  void grid_cell_corners_2d(int, int);
  void grid_cell_corners_3d(int, int, int);
};

}    // namespace LAMMPS_NS

#endif
#endif
