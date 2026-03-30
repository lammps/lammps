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

#ifndef LMP_IMAGE_OBJECTS_H
#define LMP_IMAGE_OBJECTS_H

#include "graphics.h"

#include <array>
#include <vector>

namespace LAMMPS_NS {
class Image;
class Region;
namespace ImageObjects {
  constexpr int RESOLUTION = 36;       // default resolution for cylindrical objects
  constexpr int DEF_ELEVEL = 3;        // default refinement level for ellipsoids
  constexpr int DEF_PLEVEL = 6;        // default refinement level for planes

  // custom data types for positions and triangles based on std::array
  using vec3 = std::array<double, 3>;
  using triangle = std::array<vec3, 3>;

  // some basic math operations for positions/vectors
  inline vec3 operator+(const vec3 &a, const vec3 &b)
  {
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
  }

  inline vec3 operator-(const vec3 &a, const vec3 &b)
  {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
  }

  inline vec3 operator*(double s, const vec3 &v)
  {
    return {s * v[0], s * v[1], s * v[2]};
  }

  inline vec3 operator*(const vec3 &v, double s)
  {
    return s * v;
  }

  class ArrowObj {
   public:
    // build an arrow template with length 1 in (1.0, 0.0, 0.0) direction as list of triangles
    ArrowObj(double _tipl = 0.2, double _tipw = 0.1, double radius = 0.1, int res = RESOLUTION);

    // draw custom arrow from unit template using center, direction, and length
    void draw(Image *, const double *, const double *, double, const double *, double, double);
    // draw custom arrow from unit template using two end points
    void draw(Image *, const double *, const double *, const double *, double, double);

   private:
    double tiplength;
    double tipwidth;
    double diameter;
    std::vector<triangle> triangles;
    std::vector<triangle> normals;
    int resolution;
  };

  class ConeObj {
   public:
    // build a truncated cone in (1.0, 0.0, 0.0) direction centered at (0.0, 0.0, 0.0)
    // with given length and top / bottom diameter as list of triangles.
    // flag is bitmap deciding whether top / bottom or side is shown.
    ConeObj(double, double, double, int flag = Graphics::CONE_ALL, int res = RESOLUTION);

    // draw triangle mesh for region. flag 1 is triangles, flag 2 is wireframe, flag 3 both
    void draw(Image *, int, const vec3 &, const vec3 &, const double *, Region *, double, double);

    // draw triangle mesh for fix
    void draw(Image *img, const vec3 &, const vec3 &, const double *, double);

   private:
    std::vector<triangle> triangles;
    std::vector<triangle> normals;
  };

  class EllipsoidObj {
   public:
    // construct (spherical) triangle mesh by refining the triangles of an icosahedron
    EllipsoidObj(int level = DEF_ELEVEL);

    // draw ellipsoid from triangle mesh for ellipsoid and superellipsoid particles
    void draw(Image *, int, const double *, const double *, const double *, const double *, double,
              double opacity = 1.0, const double *block = nullptr);

    // draw ellipsoid from triangle mesh for ellipsoid regions
    void draw(Image *, int, const double *, const double *, const double *, Region *, double,
              double opacity = 1.0);

   private:
    std::vector<triangle> triangles;
    void refine();
  };

  class PlaneObj {
   public:
    // build a plane template with four triangles extending well outside the box
    PlaneObj(int level = DEF_PLEVEL);

    // draw plane for region after transforming it.
    void draw(Image *, int, const double *, const double *, const double *, const double *,
              const double *, double, Region *reg, double, double opacity = 1.0);

   private:
    std::vector<triangle> triangles;
    void refine();
  };
}    // namespace ImageObjects
}    // namespace LAMMPS_NS

#endif
