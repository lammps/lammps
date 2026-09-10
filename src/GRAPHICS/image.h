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

#ifndef LMP_IMAGE_H
#define LMP_IMAGE_H

#include "pointers.h"

#include <array>
#include <cmath>
#include <unordered_map>

namespace LAMMPS_NS {

class Image : protected Pointers {
 public:
  // indices of the colormaps managed by this class.  The order must match
  // how DumpImage allocates and addresses them (amap/gmap/bmap).
  enum { ATOM_MAP = 0, GRID_MAP = 1, BOND_MAP = 2 };

  int width, height;          // size of image
  double theta, phi;          // view image from theta,phi
  double xctr, yctr, zctr;    // center of image in user coords
  double up[3];               // up direction in image
  double zoom;                // zoom factor
  double shiny;               // shininess of objects
  double gamma;               // gamma correction of rendered objects, 1.0 = off
  int fsaa;                   // antialiasing on or off
  int ssao;                   // SSAO on or off
  int seed;                   // RN seed for SSAO
  double ssaoint;             // strength of shading from 0 to 1
  int ssaosamples;            // SSAO samples per pixel; 0 = derived from ssaoint
  int depthcue;               // depth cueing on or off
  double depthcueint;         // strength of depth cueing from 0 to 1
  double *depthcuecolor;      // fog color; fade toward background color if null
  int depthcuestartflag;      // 1 if fading starts at a box fraction, 0 at nearest object
  double depthcuestart;       // start of fading as box fraction along the view direction
  int defocus;                // background defocus on or off
  double defocusint;          // strength of the defocus blur from 0 to 1
  int defocusstartflag;       // 1 if blurring starts at a box fraction, 0 at nearest object
  double defocusstart;        // start of blurring as box fraction along the view direction
  int outline;                // outline drawing on or off
  int outlinewidth;           // width of outlines in pixels
  double *outlinecolor;       // color of the outlines
  double *boxcolor;           // color to draw box outline with
  int background[3];          // RGB values of background
  int background2[3];         // RGB values of second background color for gradient (off if < 0.0)

  double ambientColor[3];    // light color settings (adjustable by caller)
  double keyLightColor[3];
  double fillLightColor[3];
  double backLightColor[3];

  int specularflag;            // 1 if the specular exponent is set explicitly
  int nospecular;              // 1 = disable the specular highlight entirely
  double specularHardness;     // exponent of the specular highlight
  double specularIntensity;    // strength of the specular highlight

  double metallic;         // 0.0 = dielectric ("plastic"), 1.0 = conductor ("metal")
  int finishMirror;        // 1 = mirror the surroundings, 0 = soft light from above
  double finishBand;       // extra brightness of the horizon band, 0.0 = off
  double finishWidth;      // exponent setting the width of the horizon band

  Image(class LAMMPS *, int);
  ~Image() override;
  void buffers();
  void clear();
  void merge();
  void write_JPG(FILE *);
  void write_PNG(FILE *);
  void write_TGA(FILE *, bool compressed = true);
  void write_PPM(FILE *);
  void view_params(double, double, double, double, double, double);

  void draw_sphere(const double *, const double *, double, double opacity = 1.0);
  void draw_cube(const double *, const double *, double, double opacity = 1.0);
  void draw_cylinder(const double *, const double *, const double *, double, int,
                     double opacity = 1.0);
  void draw_triangle(const double *, const double *, const double *, const double *,
                     double opacity = 1.0);
  void draw_trinorm(const double *, const double *, const double *, const double *, const double *,
                    const double *, const double *, const double *, const double *,
                    double opacity = 1.0);
  void draw_box(double (*)[3], double, double opacity = 1.0);
  void draw_axes(double (*)[3], double, double opacity = 1.0);
  void draw_pixmap(const double *, int, int, const unsigned char *, double *, double scale = 1.0,
                   double opacity = 1.0);
  void draw_pixmap(int, int, int, int, const unsigned char *, double *, double scale = 1.0,
                   double opacity = 1.0, double depth = 0.0);

  int map_dynamic(int);
  int map_reset(int, int, char **);
  int map_minmax(int, double, double);
  int map_info(int, double &, double &, bool &);
  double *map_value2color(int, double);

  int addcolor(const std::string &, double, double, double);
  double *element2color(const std::string &);
  double element2diam(const std::string &) const;
  double *color2rgb(const std::string &);
  std::string rgb2color(const double *) const;

 private:
  int me, nprocs;
  int npixels;

  class ColorMap **maps;
  int nmap;

  std::unordered_map<std::string, std::array<double, 3>> rgbcolors;

  struct elementInfo {
    double rgb[3];
    double diam;
  };
  std::unordered_map<std::string, elementInfo> elementdata;

  double *depthBuffer, *surfaceBuffer;
  double *depthcopy, *surfacecopy;
  unsigned char *imageBuffer, *rgbcopy, *writeBuffer;

  // MPI_Gatherv

  int *recvcounts, *displs;

  // constant view params

  double FOV;

  double keyLightTheta;
  double keyLightPhi;
  double fillLightTheta;
  double fillLightPhi;
  double backLightTheta;
  double backLightPhi;

  double SSAORadius;
  int SSAOSamples;
  double SSAOJitter;

  // dynamic view params

  double zdist;
  double tanPerPixel;
  double boxbounds[6];    // box bounds from the last view_params() call
  double camDir[3], camUp[3], camRight[4], camPos[3];
  double keyLightDir[3], fillLightDir[3], backLightDir[3];
  double keyHalfDir[3];

  // internal methods

  void draw_pixel(int, int, double, const double *, const double *);
  void setup_lights();
  void compute_SSAO();
  void compute_outline();
  void compute_depthcue();
  void compute_defocus();
  bool depth_minmax(double &, double &) const;
  void box_depth_minmax(double &, double &) const;

  // inline functions

  inline double saturate(double v)
  {
    if (v < 0.0)
      return 0.0;
    else if (v > 1.0)
      return 1.0;
    else
      return v;
  }

  inline double distance(double *a, double *b)
  {
    return sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) +
                (a[2] - b[2]) * (a[2] - b[2]));
  }
};

// ColorMap class

class ColorMap : protected Pointers {
 public:
  int dynamic;    // 0/1 if lo/hi bounds are static/dynamic

  ColorMap(class LAMMPS *, class Image *);
  ~ColorMap() override;
  int reset(int, char **);
  int minmax(double, double);
  int info(double &, double &, bool &);
  double *value2color(double);

 private:
  class Image *image;              // caller with color2rgb() method
  int mstyle, mrange;              // 2-letter style/range of color map
  int mlo, mhi;                    // bounds = NUMERIC or MINVALUE or MAXVALUE
  double mlovalue, mhivalue;       // user bounds if NUMERIC
  double locurrent, hicurrent;     // current bounds for this snapshot
  double mbinsize, mbinsizeinv;    // bin size for sequential color map
  double interpolate[3];           // local storage for returned RGB color

  struct MapEntry {
    int single, lo, hi;               // NUMERIC or MINVALUE or MAXVALUE
    double svalue, lvalue, hvalue;    // actual value
    double *color;                    // RGB values
  };

  MapEntry *mentry;
  int nentry;
};
}    // namespace LAMMPS_NS
#endif
