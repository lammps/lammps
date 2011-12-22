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

#ifdef DUMP_CLASS

DumpStyle(image,DumpImage)

#else

#ifndef LMP_DUMP_IMAGE_H
#define LMP_DUMP_IMAGE_H

#include "math.h"
#include "dump_custom.h"

namespace LAMMPS_NS {

class DumpImage : public DumpCustom {
 public:
  DumpImage(class LAMMPS *, int, char**);
  ~DumpImage();
  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);

 private:
  int filetype;
  int acolor,adiam;                // what determines color/diam of atoms
  double adiamvalue;               // atom diameter value
  int atomflag,bondflag;           // 0/1 for draw atoms,bonds
  int bcolor,bdiam;                // what determines color/diam of bonds
  double bdiamvalue;               // bond diameter value
  char *thetastr,*phistr;          // variables for view theta,phi
  int thetavar,phivar;             // index to theta,phi vars
  int cflag;                       // static/dynamic box center
  double cx,cy,cz;                 // fractional box center
  char *cxstr,*cystr,*czstr;       // variables for box center
  int cxvar,cyvar,czvar;           // index to box center vars
  char *upxstr,*upystr,*upzstr;    // view up vector variables
  int upxvar,upyvar,upzvar;        // index to up vector vars
  char *zoomstr,*perspstr;         // view zoom and perspective variables
  int zoomvar,perspvar;            // index to zoom,persp vars
  int boxflag,axesflag;            // 0/1 for draw box and axes
  double boxdiam,axeslen,axesdiam; // params for drawing box and axes

  int viewflag;                    // overall view is static or dynamic

  double *diamtype,*diamelement,*bdiamtype;         // per-type diameters
  double **colortype,**colorelement,**bcolortype;   // per-type colors

  class Image *image;              // class that renders each image

  double **bufcopy;                // buffer for communicating bond/atom info
  int maxbufcopy;

  void init_style();
  int modify_param(int, char **);
  void write();

  void box_center();
  void view_params();
  void box_bounds();
  void color_minmax();

  void create_image();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Invalid dump image filename

UNDOCUMENTED

E: Cannot dump JPG file

UNDOCUMENTED

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Dump image bond not allowed with no bond types

UNDOCUMENTED

E: Invalid dump image theta value

UNDOCUMENTED

E: Dump image persp option is not yet supported

UNDOCUMENTED

E: Dump image requires one snapshot per file

UNDOCUMENTED

E: Dump image cannot perform sorting

UNDOCUMENTED

E: Variable name for dump image theta does not exist

UNDOCUMENTED

E: Variable for dump image theta is invalid style

UNDOCUMENTED

E: Variable name for dump image phi does not exist

UNDOCUMENTED

E: Variable for dump image phi is invalid style

UNDOCUMENTED

E: Variable name for dump image center does not exist

UNDOCUMENTED

E: Variable for dump image center is invalid style

UNDOCUMENTED

E: Variable name for dump image zoom does not exist

UNDOCUMENTED

E: Variable for dump image zoom is invalid style

UNDOCUMENTED

E: Variable name for dump image persp does not exist

UNDOCUMENTED

E: Variable for dump image persp is invalid style

UNDOCUMENTED

E: Invalid dump image element name

UNDOCUMENTED

E: Invalid dump image zoom value

UNDOCUMENTED

E: Invalid dump image persp value

UNDOCUMENTED

E: Invalid dump image up vector

UNDOCUMENTED

E: Invalid dump image color range

UNDOCUMENTED

E: Invalid color in dump_modify command

UNDOCUMENTED

E: Illega dump_modify command

UNDOCUMENTED

E: Invalid color map in dump_modify command

UNDOCUMENTED

E: Dump modify bcolor not allowed with no bond types

UNDOCUMENTED

E: Dump modify bdiam not allowed with no bond types

UNDOCUMENTED

*/
