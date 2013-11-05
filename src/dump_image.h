/* -*- c++ -*- ----------------------------------------------------------
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

#include "dump_custom.h"

namespace LAMMPS_NS {

class DumpImage : public DumpCustom {
 public:
  int multifile_override;          // used by write_dump command

  DumpImage(class LAMMPS *, int, char**);
  virtual ~DumpImage();
  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);

 protected:
  int filetype;
  enum{PPM,JPG,PNG};

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

  int *chooseghost;                // extended choose array for comm
  double **bufcopy;                // buffer for communicating bond/atom info
  int maxbufcopy;

  virtual void init_style();
  int modify_param(int, char **);
  void write();

  void box_center();
  void view_params();
  void box_bounds();

  void create_image();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Invalid dump image filename

The file produced by dump image cannot be binary and must
be for a single processor.

E: Support for writing images in JPEG format not included

LAMMPS was not built with the -DLAMMPS_JPEG switch in the Makefile.

E: Support for writing images in PNG format not included

LAMMPS was not built with the -DLAMMPS_PNG switch in the Makefile.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Dump image bond not allowed with no bond types

Self-explanatory.

E: Invalid dump image theta value

Theta must be between 0.0 and 180.0 inclusive.

E: Dump image persp option is not yet supported

Self-explanatory.

E: Dump image requires one snapshot per file

Use a "*" in the filename.

E: Dump image cannot perform sorting

Self-explanatory.

E: Variable name for dump image theta does not exist

Self-explanatory.

E: Variable for dump image theta is invalid style

Must be an equal-style variable.

E: Variable name for dump image phi does not exist

Self-explanatory.

E: Variable for dump image phi is invalid style

Must be an equal-style variable.

E: Variable name for dump image center does not exist

Self-explanatory.

E: Variable for dump image center is invalid style

Must be an equal-style variable.

E: Variable name for dump image zoom does not exist

Self-explanatory.

E: Variable for dump image zoom is invalid style

Must be an equal-style variable.

E: Variable name for dump image persp does not exist

Self-explanatory.

E: Variable for dump image persp is invalid style

Must be an equal-style variable.

E: Invalid dump image element name

The specified element name was not in the standard list of elements.
See the dump_modify doc page.

E: Invalid dump image zoom value

Zoom value must be > 0.0.

E: Invalid dump image persp value

Persp value must be >= 0.0.

E: Invalid color in dump_modify command

The specified color name was not in the list of recognized colors.
See the dump_modify doc page.

E: Dump modify bcolor not allowed with no bond types

Self-explanatory.

E: Dump modify bdiam not allowed with no bond types

Self-explanatory.

*/
