/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

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

class DumpImage : public DumpCustom {
 public:
  int multifile_override;    // used by write_dump command

  DumpImage(class LAMMPS *, int, char **);
  ~DumpImage() override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

 protected:
  int filetype;
  enum { PPM, JPG, PNG };

  int atomflag;         // 0/1 for draw atoms
  int acolor, adiam;    // what determines color/diam of atoms
  double adiamvalue;    // atom diameter value

  int lineflag;                   // 0/1 for draw atoms as lines
  int lcolor, ldiam;              // what determines color/diam of lines
  double ldiamvalue;              // line diameter value
  int triflag;                    // 0/1 for draw atoms as triangles
  int tcolor, tstyle;             // what determines color/style of tris
  double tdiamvalue;              // tri edge diameter value
  int bodyflag;                   // 0/1 for draw atoms as bodies
  int bodycolor;                  // what determines color of bodies
  double bodyflag1, bodyflag2;    // user-specified params for drawing bodies
  int fixflag;                    // 0/1 to draw what fix provides
  int fixcolor;                   // what determines color of fix objects
  double fixflag1, fixflag2;      // user-specified params for fix objects

  int bondflag;         // 0/1 for draw bonds
  int bcolor, bdiam;    // what determines color/diam of bonds
  double bdiamvalue;    // bond diameter value

  int extraflag;                        // 0/1 for any of line/tri/body flag set
  char *thetastr, *phistr;              // variables for view theta,phi
  int thetavar, phivar;                 // index to theta,phi vars
  int cflag;                            // static/dynamic box center
  double cx, cy, cz;                    // fractional box center
  char *cxstr, *cystr, *czstr;          // variables for box center
  int cxvar, cyvar, czvar;              // index to box center vars
  char *upxstr, *upystr, *upzstr;       // view up vector variables
  int upxvar, upyvar, upzvar;           // index to up vector vars
  char *zoomstr;                        // view zoom variable name
  int zoomvar;                          // index to zoom variable
  int boxflag, axesflag;                // 0/1 for draw box and axes
  double boxdiam, axeslen, axesdiam;    // params for drawing box and axes
  int subboxflag;
  double subboxdiam;

  int viewflag;    // overall view is static or dynamic

  double *diamtype, *diamelement, *bdiamtype;          // per-type diameters
  double **colortype, **colorelement, **bcolortype;    // per-type colors

  class AtomVecLine *avec_line;    // ptrs to atom style (sub)classes
  class AtomVecTri *avec_tri;
  class AtomVecBody *avec_body;

  class Fix *fixptr;    // ptr to Fix that provides image data

  class Image *image;    // class that renders each image

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
};

}    // namespace LAMMPS_NS

#endif
#endif
