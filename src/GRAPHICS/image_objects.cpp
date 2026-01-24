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

#include "image_objects.h"

#include "image.h"
#include "math_const.h"
#include "math_extra.h"
#include "region.h"

#include <array>
#include <cmath>
#include <utility>
#include <vector>

using namespace LAMMPS_NS;
using namespace ImageObjects;

namespace {

using LAMMPS_NS::MathConst::MY_2PI;
constexpr double RADOVERLAP = 0.01;
constexpr double SMALL = 1.0e-10;

// helper functions for generating and transforming triangle meshes

// dot product of two vectors
inline double vec3dot(const vec3 &a, const vec3 &b)
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// cross product of two vectors
inline vec3 vec3cross(const vec3 &a, const vec3 &b)
{
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

// length of vector
inline double vec3len(const vec3 &v)
{
  return sqrt(vec3dot(v, v));
}

// return normalized vector
inline vec3 vec3norm(const vec3 &v)
{
  double n = vec3len(v);
  return (n > 0.0) ? (1.0 / n) * v : vec3{0.0, 0.0, 0.0};
}

// scale factor to move a position to the surface of an ellipsoid with given shape parameters
inline double radscale(const double *shape, const vec3 &pos)
{
  return sqrt(1.0 /
              (pos[0] / shape[0] * pos[0] / shape[0] + pos[1] / shape[1] * pos[1] / shape[1] +
               pos[2] / shape[2] * pos[2] / shape[2]));
}

// re-orient list of triangles to point along "dir", then scale and translate it.
std::vector<triangle> transform(const std::vector<triangle> &triangles, const vec3 &dir,
                                const vec3 &offs, double len, double width)
{
  // customized vector
  std::vector<triangle> newtriangles;

  // normalize direction vector
  vec3 u = vec3norm(dir);

  // vector is too short. can't draw anything. return empty list
  if (vec3len(u) < SMALL) return newtriangles;

  // construct orthonormal basis around direction vector
  vec3 a = (std::fabs(u[0]) < 0.9) ? vec3{1.0, 0.0, 0.0} : vec3{0.0, 1.0, 0.0};
  vec3 v = vec3norm(vec3cross(u, a));
  vec3 w = vec3cross(u, v);

  // now process the template triangles and return the transformed list
  newtriangles.reserve(triangles.size());
  for (const auto &tri : triangles) {
    vec3 p1 = (len * tri[0][0] * u) + (width * tri[0][1] * v) + (width * tri[0][2] * w) + offs;
    vec3 p2 = (len * tri[1][0] * u) + (width * tri[1][1] * v) + (width * tri[1][2] * w) + offs;
    vec3 p3 = (len * tri[2][0] * u) + (width * tri[2][1] * v) + (width * tri[2][2] * w) + offs;
    newtriangles.push_back({p1, p2, p3});
  }
  return newtriangles;
}
}    // namespace

// construct an arrow from primitives, mostly triangles and a cylinder, and draw them

// construct arrow template by placing sets of triangles with two corners on a circle at "mid"
// and the third corner in the center either at "mid" or at "tip". A third set of triangles
// it at "bot". The resolution parameter determines how many triangles per set (36 by default)
// "bot" to "tip" is 1.0, "tiplength is "mid" to "tip". "diameter" is the width at "bot"
// "tipwidth" is the additional width at "mid".
//
//          |\
// |--------| \
// |--------| /
//          |/
// ^        ^  ^
// bot    mid tip

ArrowObj::ArrowObj(double _tipl, double _tipw, double radius, int res) :
    tiplength(_tipl), tipwidth(_tipw), diameter(2.0 * radius), resolution(res)
{
  triangles.clear();

  // we want at least 2 iterations.
  if (res < 2) return;

  vec3 tip{0.5, 0.0, 0.0};
  vec3 mid{0.5 - tiplength, 0.0, 0.0};
  vec3 bot{-0.5, 0.0, 0.0};

  // construct list of triangles for the tip of the arrow. p1, p2 are the points on the "rim".

  const double radinc = MY_2PI / resolution;
  vec3 p1{0.5 - tiplength, 0.0, 0.0};
  vec3 p2{0.5 - tiplength, 0.0, 0.0};
  for (int i = 0; i < resolution; ++i) {
    p1[1] = (radius + tipwidth) * sin(radinc * i - RADOVERLAP);
    p1[2] = (radius + tipwidth) * cos(radinc * i - RADOVERLAP);
    p2[1] = (radius + tipwidth) * sin(radinc * (i + 1));
    p2[2] = (radius + tipwidth) * cos(radinc * (i + 1));
    triangles.emplace_back(triangle{p2, tip, p1});
    triangles.emplace_back(triangle{p2, mid, p1});
  }

  // construct list of triangles for the cap at the bottom

  p1[0] = -0.5;
  p2[0] = -0.5;
  for (int i = 0; i < resolution; ++i) {
    p1[1] = radius * sin(radinc * i - RADOVERLAP);
    p1[2] = radius * cos(radinc * i - RADOVERLAP);
    p2[1] = radius * sin(radinc * (i + 1));
    p2[2] = radius * cos(radinc * (i + 1));
    triangles.emplace_back(triangle{p2, bot, p1});
  }
}

// draw custom arrow from unit template using center, direction, and length
void ArrowObj::draw(Image *img, const double *color, const double *center, double length,
                    const double *data, double scale, double opacity)
{
  // nothing to draw
  if (!triangles.size()) return;

  // transform the template into the arrow object we want to draw

  vec3 dir{data[0], data[1], data[2]};
  double lscale = vec3len(dir) * length;
  double wscale = scale / diameter;

  auto arrow = transform(triangles, dir, {center[0], center[1], center[2]}, lscale, wscale);

  // nothing to draw
  if (!arrow.size()) return;

  // draw tip and bottom from list of triangles
  for (const auto &tri : arrow)
    img->draw_triangle(tri[0].data(), tri[1].data(), tri[2].data(), color, opacity);

  // infer cylinder end points for body from list of triangles
  // (middle corner of all triangles in the the second and last set of triangles)
  if ((int) arrow.size() > resolution + 2)
    img->draw_cylinder(arrow[1][1].data(), arrow[arrow.size() - 1][1].data(), color, scale, 0,
                       opacity);
}

// draw custom arrow from unit template using center, direction, and length
void ArrowObj::draw(Image *img, const double *color, const double *bottom, const double *tip,
                    double scale, double opacity)
{
  // nothing to draw
  if (!triangles.size()) return;

  // transform the template into the arrow object we want to draw

  vec3 dir{vec3{tip[0], tip[1], tip[2]} - vec3{bottom[0], bottom[1], bottom[2]}};
  vec3 center{0.5 * dir + vec3{bottom[0], bottom[1], bottom[2]}};
  double lscale = vec3len(dir);
  double wscale = scale / diameter;

  auto arrow = transform(triangles, dir, {center[0], center[1], center[2]}, lscale, wscale);

  // nothing to draw
  if (!arrow.size()) return;

  // draw tip and bottom from list of triangles
  for (const auto &tri : arrow)
    img->draw_triangle(tri[0].data(), tri[1].data(), tri[2].data(), color, opacity);

  // infer cylinder end points for body from list of triangles
  // (middle corner of all triangles in the the second and last set of triangles)
  if ((int) arrow.size() > resolution + 2)
    img->draw_cylinder(arrow[1][1].data(), arrow[arrow.size() - 1][1].data(), color, scale, 0,
                       opacity);
}

// construct a truncated cone from triangles and draw them

// we have two circles and place triangles that connect from bottom to top and back
// where the top of the triangle alternates direction. the caps on either end use the
// same circle coordinates but the tip is the center of the object.
// as an optimization we skip triangles where the bottom is on the circle when the
// diameter on either side of the cone is zero.
// a cylinder is just a special case of a cone with both radii of the same value.
//
// |\
// | \
// |  |  _ center
// |  |
// | /
// |/
// ^  ^
//bot top

ConeObj::ConeObj(double length, double topwidth, double botwidth, int flag, int resolution)
{
  triangles.clear();

  // we want at least 2 iterations.
  if (resolution < 2) return;

  // store settings for cone

  bool dotop = (flag & Graphics::CONE_TOP) > 0;
  bool dobot = (flag & Graphics::CONE_BOT) > 0;
  bool doside = (flag & Graphics::CONE_SIDE) > 0;

  vec3 top{0.5 * length, 0.0, 0.0};
  vec3 bot{-0.5 * length, 0.0, 0.0};

  // construct list of triangles

  const double radinc = MY_2PI / resolution;
  vec3 p1top{top};
  vec3 p2top{top};
  vec3 p1bot{bot};
  vec3 p2bot{bot};

  for (int i = 0; i < resolution; ++i) {
    if (topwidth > 0.0) {
      p1top[1] = topwidth * sin(radinc * i - RADOVERLAP);
      p1top[2] = topwidth * cos(radinc * i - RADOVERLAP);
      p2top[1] = topwidth * sin(radinc * (i + 1));
      p2top[2] = topwidth * cos(radinc * (i + 1));
      // cap on top
      if (dotop) triangles.emplace_back(triangle{p1top, top, p2top});
    }
    if (botwidth > 0.0) {
      p1bot[1] = botwidth * sin(radinc * i - RADOVERLAP);
      p1bot[2] = botwidth * cos(radinc * i - RADOVERLAP);
      p2bot[1] = botwidth * sin(radinc * (i + 1));
      p2bot[2] = botwidth * cos(radinc * (i + 1));
      // cap at bottom
      if (dobot) triangles.emplace_back(triangle{p1bot, bot, p2bot});
    }
    // side
    if (doside) {
      if (topwidth > 0.0) triangles.emplace_back(triangle{p1top, p1bot, p2top});
      if (botwidth > 0.0) triangles.emplace_back(triangle{p1bot, p2bot, p2top});
    }
  }
}

// draw triangle mesh for region. flag 1 is triangles, flag 2 is wireframe, flag 3 both

void ConeObj::draw(Image *img, int flag, const vec3 &dir, const vec3 &mid, const double *color,
                   Region *reg, double diameter, double opacity)
{
  // nothing to draw
  if (!triangles.size()) return;

  // rotate to selected axis and translate from origin to original center
  // no need of scaling here since length and width was already applied during construction
  auto cone = transform(triangles, dir, mid, 1.0, 1.0);

  // nothing to draw
  if (!cone.size()) return;

  int n = 0;
  for (auto &tri : cone) {
    // apply region rotation and translation
    reg->forward_transform(tri[0][0], tri[0][1], tri[0][2]);
    reg->forward_transform(tri[1][0], tri[1][1], tri[1][2]);
    reg->forward_transform(tri[2][0], tri[2][1], tri[2][2]);

    // draw triangle
    if (flag & 1) img->draw_triangle(tri[0].data(), tri[1].data(), tri[2].data(), color, opacity);

    // draw wireframe
    if (flag & 2) {
      // draw bottom rim and straight lines from bottom to top
      img->draw_cylinder(tri[0].data(), tri[1].data(), color, diameter, 3, opacity);
      // only draw top rim by picking coordinates from every other triangle
      ++n;
      if (n & 1) img->draw_cylinder(tri[0].data(), tri[2].data(), color, diameter, 3, opacity);
    }
  }
}

// draw triangle mesh for fix.

void ConeObj::draw(Image *img, const vec3 &bot, const vec3 &top, const double *color,
                   double opacity)
{
  // nothing to draw
  if (!triangles.size()) return;

  vec3 mid{0.5 * (top + bot)};
  vec3 dir{top - bot};
  double length = vec3len(dir);
  dir = vec3norm(dir);

  // rotate to selected axis and translate from origin to original center
  // no need of scaling here since length and width was already applied during construction
  auto cone = transform(triangles, dir, mid, length, 1.0);

  // nothing to draw
  if (!cone.size()) return;

  for (auto &tri : cone) {
    // draw triangle
    img->draw_triangle(tri[0].data(), tri[1].data(), tri[2].data(), color, opacity);
  }
}

/****************************************************************************
 * Refine triangle mesh by replacing each triangle with four triangles.
 * Compute the new positions so they are located on a sphere with radius 1.
 *    /\            /\
 *   /  \          /_ \
 *  /    \   -->  /\  /\
 * /______\      /__\/__\
 ***************************************************************************/

void EllipsoidObj::refine()
{
  std::vector<triangle> newlist;
  for (const auto &tri : triangles) {
    vec3 posa = vec3norm(tri[0] + tri[2]);
    vec3 posb = vec3norm(tri[0] + tri[1]);
    vec3 posc = vec3norm(tri[1] + tri[2]);
    newlist.push_back({tri[0], posb, posa});
    newlist.push_back({posb, tri[1], posc});
    newlist.push_back({posa, posb, posc});
    newlist.push_back({posa, posc, tri[2]});
  }
  triangles = std::move(newlist);
}

// Construct and draw an ellipsoid from primitives, triangles and cylinders.
// Build a triangle mesh by refinining the triangles of an octahedron

EllipsoidObj::EllipsoidObj(int level)
{
  // Define edges of an octahedron to approximate a sphere of radius 1 around the origin.
  constexpr vec3 OCT1 = {-1.0, 0.0, 0.0};
  constexpr vec3 OCT2 = {1.0, 0.0, 0.0};
  constexpr vec3 OCT3 = {0.0, -1.0, 0.0};
  constexpr vec3 OCT4 = {0.0, 1.0, 0.0};
  constexpr vec3 OCT5 = {0.0, 0.0, -1.0};
  constexpr vec3 OCT6 = {0.0, 0.0, 1.0};

  // define level 1 octahedron triangle mesh, normals pointing away from the center.
  triangles = {{OCT5, OCT4, OCT1}, {OCT2, OCT4, OCT5}, {OCT6, OCT4, OCT2}, {OCT1, OCT4, OCT6},
               {OCT1, OCT3, OCT5}, {OCT5, OCT3, OCT2}, {OCT2, OCT3, OCT6}, {OCT6, OCT3, OCT1}};

  // refine the list of triangles to the desired level
  for (int i = 1; i < level; ++i) refine();
}

// draw method for drawing ellipsoids from a region which has its own transformation function
void EllipsoidObj::draw(Image *img, int flag, const double *color, const double *center,
                        const double *shape, Region *reg, double diameter, double opacity)
{
  // select between triangles or cylinders
  bool doframe = false;
  bool dotri = false;
  if (flag == 1) dotri = true;
  if (flag == 2) doframe = true;
  if (diameter <= 0.0) doframe = false;
  if (!dotri && !doframe) return;    // nothing to do

  // optimization: just draw a sphere if a filled surface is requested and the object is a sphere
  if (dotri && (shape[0] == shape[1]) && (shape[0] == shape[2])) {
    img->draw_sphere(center, color, 2.0 * shape[0], opacity);
    return;
  }

  // nothing to draw
  if (!triangles.size()) return;

  // draw triangles

  const vec3 offs{center[0], center[1], center[2]};
  for (auto tri : triangles) {

    // set shape and move
    tri[0] = tri[0] * radscale(shape, tri[0]) + offs;
    reg->forward_transform(tri[0][0], tri[0][1], tri[0][2]);
    tri[1] = tri[1] * radscale(shape, tri[1]) + offs;
    reg->forward_transform(tri[1][0], tri[1][1], tri[1][2]);
    tri[2] = tri[2] * radscale(shape, tri[2]) + offs;
    reg->forward_transform(tri[2][0], tri[2][1], tri[2][2]);

    if (dotri) img->draw_triangle(tri[0].data(), tri[1].data(), tri[2].data(), color, opacity);
    if (doframe) {
      img->draw_cylinder(tri[0].data(), tri[1].data(), color, diameter, 3, opacity);
      img->draw_cylinder(tri[0].data(), tri[2].data(), color, diameter, 3, opacity);
      img->draw_cylinder(tri[1].data(), tri[2].data(), color, diameter, 3, opacity);
    }
  }
}

// draw method for drawing ellipsoids from per-atom data which has a quaternion
// and the shape list to define the orientation and stretch
void EllipsoidObj::draw(Image *img, int flag, const double *color, const double *center,
                        const double *shape, const double *quat, double diameter, double opacity)
{
  // select between triangles or cylinders or both
  bool doframe = true;
  bool dotri = true;
  if (flag == 1) doframe = false;
  if (flag == 2) dotri = false;
  if (diameter <= 0.0) doframe = false;
  if (!dotri && !doframe) return;    // nothing to do

  double p[3][3];
  vec3 e1, e2, e3;
  const vec3 offs{center[0], center[1], center[2]};

  // optimization: just draw a sphere if a filled surface is requested and the object is a sphere
  if (dotri && (shape[0] == shape[1]) && (shape[0] == shape[2])) {
    img->draw_sphere(center, color, 2.0 * shape[0], opacity);
    return;
  }

  // nothing to draw
  if (!triangles.size()) return;

  // get rotation matrix for body frame to box frame
  MathExtra::quat_to_mat(quat, p);

  // draw triangles and edges as requested, work on copy of triangle since we modify it
  for (auto tri : triangles) {

    if (dotri) {
      // set shape by shifting each corner to the surface
      for (int i = 0; i < 3; ++i) {
        auto &t = tri[i];
        t = radscale(shape, t) * t;
      }

      // rotate
      MathExtra::matvec(p, tri[0].data(), e1.data());
      MathExtra::matvec(p, tri[1].data(), e2.data());
      MathExtra::matvec(p, tri[2].data(), e3.data());

      // translate
      e1 = e1 + offs;
      e2 = e2 + offs;
      e3 = e3 + offs;

      img->draw_triangle(e1.data(), e2.data(), e3.data(), color, opacity);
    }

    if (doframe) {
      // set shape
      for (int i = 0; i < 3; ++i) {
        auto &t = tri[i];
        if (dotri) {
          // shift the cylinder positions inward by their diameter when using cylinders and
          // triangles together for a smoother surface to avoid increasing the final size
          double shapeplus[3] = {shape[0] - diameter, shape[1] - diameter, shape[1] - diameter};
          t = radscale(shapeplus, t) * t;
        } else {
          t = radscale(shape, t) * t;
        }
      }

      // rotate
      MathExtra::matvec(p, tri[0].data(), e1.data());
      MathExtra::matvec(p, tri[1].data(), e2.data());
      MathExtra::matvec(p, tri[2].data(), e3.data());

      // translate
      e1 = e1 + offs;
      e2 = e2 + offs;
      e3 = e3 + offs;
      img->draw_cylinder(e1.data(), e2.data(), color, diameter, 3, opacity);
      img->draw_cylinder(e2.data(), e3.data(), color, diameter, 3, opacity);
      img->draw_cylinder(e3.data(), e1.data(), color, diameter, 3, opacity);
    }
  }
}

/***********************************************************************
 * refine triangle mesh by replacing each triangle with four triangles
 *    /\            /\
 *   /  \          /__\
 *  /    \   -->  /\  /\
 * /______\      /__\/__\
***********************************************************************/

void PlaneObj::refine()
{
  std::vector<triangle> newlist;
  for (const auto &tri : triangles) {
    vec3 posa = 0.5 * (tri[0] + tri[2]);
    vec3 posb = 0.5 * (tri[0] + tri[1]);
    vec3 posc = 0.5 * (tri[1] + tri[2]);
    newlist.push_back({tri[0], posb, posa});
    newlist.push_back({posb, tri[1], posc});
    newlist.push_back({posa, posb, posc});
    newlist.push_back({posa, posc, tri[2]});
  }
  triangles = std::move(newlist);
}

// construct a plane from many triangles (so we can truncate it to the box dimensions)

PlaneObj::PlaneObj(int level)
{
  // define edges and center of a square
  constexpr vec3 SQ1 = {0.0, 1.0, 1.0};
  constexpr vec3 SQ2 = {0.0, 1.0, -1.0};
  constexpr vec3 SQ3 = {0.0, -1.0, -1.0};
  constexpr vec3 SQ4 = {0.0, -1.0, 1.0};
  constexpr vec3 CEN = {0.0, 0.0, 0.0};

  // define unit plane with norm (1.0,0.0,0.0) from four triangles
  triangles = {{SQ2, CEN, SQ1}, {SQ3, CEN, SQ2}, {SQ4, CEN, SQ3}, {SQ1, CEN, SQ4}};

  // refine the list of triangles to the desired level
  for (int i = 1; i < level; ++i) refine();
}

// draw method for drawing planes from a region which has its own transformation function

void PlaneObj::draw(Image *img, int flag, const double *color, const double *center,
                    const double *norm, const double *boxlo, const double *boxhi, double scale,
                    Region *reg, double diameter, double opacity)
{
  // select between triangles or cylinders
  bool doframe = false;
  bool dotri = false;
  if (flag == 1) dotri = true;
  if (flag == 2) doframe = true;
  if (diameter <= 0.0) doframe = false;
  if (!dotri && !doframe) return;    // nothing to do

  // nothing to draw

  if (!triangles.size()) return;

  // draw triangles after scaling and shifting the mesh

  const vec3 dir{norm[0], norm[1], norm[2]};
  const vec3 offs{center[0], center[1], center[2]};
  auto plane = transform(triangles, dir, offs, scale, scale);

  for (auto tri : plane) {

    // rotate and translate

    reg->forward_transform(tri[0][0], tri[0][1], tri[0][2]);
    reg->forward_transform(tri[1][0], tri[1][1], tri[1][2]);
    reg->forward_transform(tri[2][0], tri[2][1], tri[2][2]);

    // skip drawing triangle if all corners are outside the box in one direction

    int n = 0;
    for (int i = 0; i < 3; ++i) {
      if (((tri[0][i] < boxlo[i]) || (tri[0][i] > boxhi[i])) &&
          ((tri[1][i] < boxlo[i]) || (tri[1][i] > boxhi[i])) &&
          ((tri[2][i] < boxlo[i]) || (tri[2][i] > boxhi[i])))
        ++n;
    }
    if (n) continue;

    if (dotri) img->draw_triangle(tri[0].data(), tri[1].data(), tri[2].data(), color, opacity);
    if (doframe) {
      img->draw_cylinder(tri[0].data(), tri[1].data(), color, diameter, 3, opacity);
      img->draw_cylinder(tri[0].data(), tri[2].data(), color, diameter, 3, opacity);
      img->draw_cylinder(tri[1].data(), tri[2].data(), color, diameter, 3, opacity);
    }
  }
}
