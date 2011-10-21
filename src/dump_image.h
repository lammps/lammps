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
  int acolor,adiam;
  double adiamvalue;
  int atomflag,bondflag;
  int bcolor,bdiam;
  double bdiamvalue;
  int width,height;
  double theta,phi;
  char *thetastr,*phistr;
  int thetavar,phivar;
  int cflag;
  double cx,cy,cz;
  char *cxstr,*cystr,*czstr;
  int cxvar,cyvar,czvar;
  double up[3];
  char *upxstr,*upystr,*upzstr;
  int upxvar,upyvar,upzvar;
  double zoom,persp;
  char *zoomstr,*perspstr;
  int zoomvar,perspvar;
  int boxflag,axesflag;
  double boxdiam,axeslen,axesdiam;
  double shiny;
  int ssao,seed;
  double ssaoint;

  int npixels,viewflag;
  double *depthBuffer,*surfaceBuffer;
  double *depthcopy,*surfacecopy;
  char *imageBuffer,*rgbcopy,*writeBuffer;

  double **bufcopy;
  int maxbufcopy;

  // constant view params

  double FOV;
  double ambientColor[3];

  double keyLightTheta;
  double keyLightPhi;
  double keyLightColor[3];

  double fillLightTheta;
  double fillLightPhi;
  double fillLightColor[3];

  double backLightTheta;
  double backLightPhi;
  double backLightColor[3];

  double specularHardness;
  double specularIntensity;

  double SSAORadius;
  int SSAOSamples;
  double SSAOJitter;

  // dynamic view params

  double xctr,yctr,zctr,zdist;
  double tanPerPixel;
  double camDir[3],camUp[3],camRight[4],camPos[3];
  double keyLightDir[3],fillLightDir[3],backLightDir[3];
  double keyHalfDir[3];

  // dump_modify values

  int ncolors;
  char **username;
  double **userrgb;

  double *diamtype,*diamelement,*bdiamtype;
  double **colortype,**colorelement,**bcolortype;

  double *boxcolor;
  int background[3];

  // color map

  int mstyle,mrange;               // 2-letter style/range of color map
  int mlo,mhi;                     // bounds = NUMERIC or MINVALUE or MAXVALUE
  double mlovalue,mhivalue;        // user bounds if NUMERIC
  double locurrent,hicurrent;      // current bounds for this snapshot
  double mbinsize,mbinsizeinv;     // bin size for sequential color map

  struct MapEntry {
    int single,lo,hi;              // NUMERIC or MINVALUE or MAXVALUE
    double svalue,lvalue,hvalue;   // actual value
    double *color;                 // RGB values
  };

  MapEntry *mentry;
  int nentry;
  double interpolate[3];

  class RanMars *random;

  void init_style();
  int modify_param(int, char **);
  void write();

  void box_center();
  void view_params();
  void box_bounds();
  void color_minmax();

  void create_image();

  // rasterizing methods

  void draw_sphere(double *, double *, double);
  void draw_cylinder(double *, double *, double *, double, int);
  void draw_pixel(int, int, double, double *, double*);
  void compute_SSAO();
  void write_JPG();
  void write_PPM();

  double *value2color(double);
  double *color2rgb(char *);
  double *element2color(char *);
  double element2diam(char *);

  // inlined functions

  inline double saturate(double v) {
    if (v < 0.0) return 0.0;
    else if (v > 1.0) return 1.0;
    else return v;
  }

  inline double distance(double* a, double* b) {
    return sqrt((a[0] - b[0]) * (a[0] - b[0]) + 
		(a[1] - b[1]) * (a[1] - b[1]) + 
		(a[2] - b[2]) * (a[2] - b[2]));
  }
};

}

#endif
#endif
