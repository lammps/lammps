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

// scale factor to move a position to the surface of a superellipsoid with given parameters
inline double superscale(const double *shape, const double *block, const vec3 &pos)
{
  double a = pow(fabs(pos[0] / shape[0]), block[1]) + pow(fabs(pos[1] / shape[1]), block[1]);
  double b = pow(fabs(pos[2] / shape[2]), block[0]);
  return pow(pow(a, block[0] / block[1]) + b, -1.0 / block[0]);
}

// compute surface normal for a superellipsoid at a point on the unit sphere.
// The superellipsoid surface is: ((|x/a|^e + |y/b|^e)^(n/e) + |z/c|^n)^(1/n) = 1
// where n = block[0], e = block[1], a = shape[0], b = shape[1], c = shape[2].
// The normal is proportional to the gradient of the implicit function:
//   nx = u^(n/e - 1) * |x/a|^(e-1) * sign(x) / a
//   ny = u^(n/e - 1) * |y/b|^(e-1) * sign(y) / b
//   nz = |z/c|^(n-1) * sign(z) / c
// where u = |x/a|^e + |y/b|^e.

inline vec3 supernormal(const double *shape, const double *block, const vec3 &pos)
{
  const double n = block[0], e = block[1];
  const double xa = fabs(pos[0] / shape[0]);
  const double yb = fabs(pos[1] / shape[1]);
  const double zc = fabs(pos[2] / shape[2]);

  const double u = pow(xa, e) + pow(yb, e);
  const double ufactor = (u > 0.0) ? pow(u, n / e - 1.0) : 0.0;

  double nx = (xa > 0.0) ? ufactor * pow(xa, e - 1.0) * copysign(1.0, pos[0]) / shape[0] : 0.0;
  double ny = (yb > 0.0) ? ufactor * pow(yb, e - 1.0) * copysign(1.0, pos[1]) / shape[1] : 0.0;
  double nz = (zc > 0.0) ? pow(zc, n - 1.0) * copysign(1.0, pos[2]) / shape[2] : 0.0;

  return vec3norm({nx, ny, nz});
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

// re-orient list of per-vertex normals to point along "dir" (rotation only, no scaling/translation)
std::vector<triangle> transform_normals(const std::vector<triangle> &norms, const vec3 &dir)
{
  std::vector<triangle> newnormals;

  // normalize direction vector
  vec3 u = vec3norm(dir);

  // vector is too short. can't draw anything. return empty list
  if (vec3len(u) < SMALL) return newnormals;

  // construct orthonormal basis around direction vector
  vec3 a = (std::fabs(u[0]) < 0.9) ? vec3{1.0, 0.0, 0.0} : vec3{0.0, 1.0, 0.0};
  vec3 v = vec3norm(vec3cross(u, a));
  vec3 w = vec3cross(u, v);

  // now process the template normals and return the rotated list
  newnormals.reserve(norms.size());
  for (const auto &n : norms) {
    vec3 n1 = vec3norm((n[0][0] * u) + (n[0][1] * v) + (n[0][2] * w));
    vec3 n2 = vec3norm((n[1][0] * u) + (n[1][1] * v) + (n[1][2] * w));
    vec3 n3 = vec3norm((n[2][0] * u) + (n[2][1] * v) + (n[2][2] * w));
    newnormals.push_back({n1, n2, n3});
  }
  return newnormals;
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
  normals.clear();

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
    // tip cone: radial normals at rim,
    //           use mid-point rim normal with a little tilt for a mostly "pointy" tip
    vec3 n1 = vec3norm({0.0, p2[1], p2[2]});
    vec3 n2 = vec3norm({0.25, 0.5 * (p1[1] + p2[1]), 0.5 * (p1[2] + p2[2])});
    vec3 n3 = vec3norm({0.0, p1[1], p1[2]});
    triangles.emplace_back(triangle{p2, tip, p1});
    normals.emplace_back(triangle{n1, n2, n3});
    // flat disc at base of tip: normal points in -x direction
    vec3 nflat = {-1.0, 0.0, 0.0};
    triangles.emplace_back(triangle{p2, mid, p1});
    normals.emplace_back(triangle{nflat, nflat, nflat});
  }

  // construct list of triangles for the cap at the bottom

  p1[0] = -0.5;
  p2[0] = -0.5;
  for (int i = 0; i < resolution; ++i) {
    p1[1] = radius * sin(radinc * i - RADOVERLAP);
    p1[2] = radius * cos(radinc * i - RADOVERLAP);
    p2[1] = radius * sin(radinc * (i + 1));
    p2[2] = radius * cos(radinc * (i + 1));
    // bottom cap: normal points in -x direction
    vec3 nflat = {-1.0, 0.0, 0.0};
    triangles.emplace_back(triangle{p2, bot, p1});
    normals.emplace_back(triangle{nflat, nflat, nflat});
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
  auto norms = transform_normals(normals, dir);

  // nothing to draw
  if (!arrow.size()) return;

  // draw tip and bottom from list of triangles with per-vertex normals
  for (size_t i = 0; i < arrow.size(); ++i)
    img->draw_trinorm(arrow[i][0].data(), arrow[i][1].data(), arrow[i][2].data(),
                      norms[i][0].data(), norms[i][1].data(), norms[i][2].data(), color, color,
                      color, opacity);

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
  auto norms = transform_normals(normals, dir);

  // nothing to draw
  if (!arrow.size()) return;

  // draw tip and bottom from list of triangles with per-vertex normals
  for (size_t i = 0; i < arrow.size(); ++i)
    img->draw_trinorm(arrow[i][0].data(), arrow[i][1].data(), arrow[i][2].data(),
                      norms[i][0].data(), norms[i][1].data(), norms[i][2].data(), color, color,
                      color, opacity);

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
  normals.clear();

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

  // cap normal directions
  vec3 ntop = {1.0, 0.0, 0.0};
  vec3 nbot = {-1.0, 0.0, 0.0};

  for (int i = 0; i < resolution; ++i) {
    if (topwidth > 0.0) {
      p1top[1] = topwidth * sin(radinc * i - RADOVERLAP);
      p1top[2] = topwidth * cos(radinc * i - RADOVERLAP);
      p2top[1] = topwidth * sin(radinc * (i + 1));
      p2top[2] = topwidth * cos(radinc * (i + 1));
      // cap on top
      if (dotop) {
        triangles.emplace_back(triangle{p1top, top, p2top});
        normals.emplace_back(triangle{ntop, ntop, ntop});
      }
    }
    if (botwidth > 0.0) {
      p1bot[1] = botwidth * sin(radinc * i - RADOVERLAP);
      p1bot[2] = botwidth * cos(radinc * i - RADOVERLAP);
      p2bot[1] = botwidth * sin(radinc * (i + 1));
      p2bot[2] = botwidth * cos(radinc * (i + 1));
      // cap at bottom
      if (dobot) {
        triangles.emplace_back(triangle{p1bot, bot, p2bot});
        normals.emplace_back(triangle{nbot, nbot, nbot});
      }
    }
    // side: use radial normals for smooth shading
    if (doside) {
      vec3 n1top = vec3norm({0.0, p1top[1], p1top[2]});
      vec3 n2top = vec3norm({0.0, p2top[1], p2top[2]});
      vec3 n1bot = vec3norm({0.0, p1bot[1], p1bot[2]});
      vec3 n2bot = vec3norm({0.0, p2bot[1], p2bot[2]});
      if (topwidth > 0.0) {
        triangles.emplace_back(triangle{p1top, p1bot, p2top});
        normals.emplace_back(triangle{n1top, n1bot, n2top});
      }
      if (botwidth > 0.0) {
        triangles.emplace_back(triangle{p1bot, p2bot, p2top});
        normals.emplace_back(triangle{n1bot, n2bot, n2top});
      }
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
  auto norms = transform_normals(normals, dir);

  // nothing to draw
  if (!cone.size()) return;

  // get the offset from forward_transform to extract rotation-only for normals
  double ox = 0.0, oy = 0.0, oz = 0.0;
  reg->forward_transform(ox, oy, oz);

  int n = 0;
  for (size_t k = 0; k < cone.size(); ++k) {
    auto &tri = cone[k];
    // apply region rotation and translation
    reg->forward_transform(tri[0][0], tri[0][1], tri[0][2]);
    reg->forward_transform(tri[1][0], tri[1][1], tri[1][2]);
    reg->forward_transform(tri[2][0], tri[2][1], tri[2][2]);

    // apply region rotation to normals (subtract translation offset)
    auto rn = norms[k];
    reg->forward_transform(rn[0][0], rn[0][1], rn[0][2]);
    reg->forward_transform(rn[1][0], rn[1][1], rn[1][2]);
    reg->forward_transform(rn[2][0], rn[2][1], rn[2][2]);
    rn[0] = vec3norm({rn[0][0] - ox, rn[0][1] - oy, rn[0][2] - oz});
    rn[1] = vec3norm({rn[1][0] - ox, rn[1][1] - oy, rn[1][2] - oz});
    rn[2] = vec3norm({rn[2][0] - ox, rn[2][1] - oy, rn[2][2] - oz});

    // draw triangle with per-vertex normals
    if (flag & 1)
      img->draw_trinorm(tri[0].data(), tri[1].data(), tri[2].data(), rn[0].data(), rn[1].data(),
                        rn[2].data(), color, color, color, opacity);

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
  auto norms = transform_normals(normals, dir);

  // nothing to draw
  if (!cone.size()) return;

  for (size_t k = 0; k < cone.size(); ++k) {
    // draw triangle with per-vertex normals
    img->draw_trinorm(cone[k][0].data(), cone[k][1].data(), cone[k][2].data(), norms[k][0].data(),
                      norms[k][1].data(), norms[k][2].data(), color, color, color, opacity);
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
// Build a triangle mesh by refining the triangles of an icosahedron

EllipsoidObj::EllipsoidObj(int level)
{
  // Define vertices of an icosahedron to approximate a sphere of radius 1 around the origin.
  // A and B are the normalized coordinates derived from the golden ratio:
  // phi = (1 + sqrt(5)) / 2; A = 1 / sqrt(1 + phi^2); B = phi / sqrt(1 + phi^2)
  constexpr double A = 0.5257311121191336;
  constexpr double B = 0.8506508083520399;
  // clang-format off
  constexpr vec3 ICO00 = { -A,   B, 0.0};
  constexpr vec3 ICO01 = {  A,   B, 0.0};
  constexpr vec3 ICO02 = { -A,  -B, 0.0};
  constexpr vec3 ICO03 = {  A,  -B, 0.0};
  constexpr vec3 ICO04 = {0.0,  -A,   B};
  constexpr vec3 ICO05 = {0.0,   A,   B};
  constexpr vec3 ICO06 = {0.0,  -A,  -B};
  constexpr vec3 ICO07 = {0.0,   A,  -B};
  constexpr vec3 ICO08 = {  B, 0.0,  -A};
  constexpr vec3 ICO09 = {  B, 0.0,   A};
  constexpr vec3 ICO10 = { -B, 0.0,  -A};
  constexpr vec3 ICO11 = { -B, 0.0,   A};
  // clang-format on

  // define level 1 icosahedron triangle mesh, normals pointing away from the center.
  triangles = {
      {ICO00, ICO05, ICO11}, {ICO00, ICO01, ICO05}, {ICO00, ICO07, ICO01}, {ICO00, ICO10, ICO07},
      {ICO00, ICO11, ICO10}, {ICO01, ICO09, ICO05}, {ICO05, ICO04, ICO11}, {ICO11, ICO02, ICO10},
      {ICO10, ICO06, ICO07}, {ICO07, ICO08, ICO01}, {ICO03, ICO04, ICO09}, {ICO03, ICO02, ICO04},
      {ICO03, ICO06, ICO02}, {ICO03, ICO08, ICO06}, {ICO03, ICO09, ICO08}, {ICO04, ICO05, ICO09},
      {ICO02, ICO11, ICO04}, {ICO06, ICO10, ICO02}, {ICO08, ICO07, ICO06}, {ICO09, ICO01, ICO08}};

  // refine the list of triangles to the desired level
  for (int i = 1; i < level; ++i) refine();

  // Rotate the sphere mesh so that the Cartesian axes point through (or near)
  // triangle face centers rather than through vertices or edges.  This improves
  // the visual quality of ellipsoids and superellipsoids by making them appear
  // less "pointy" along the three principal axes, especially with lower
  // triangle counts.  The default orientation has the principal axes pass
  // through corners or edges.  We prefer smooth geometries for simulations
  // anyway and thus the rotation allows to get a better approximation from the
  // triangulation with lower refinement levels and thus require less
  // computational effort for creating an acceptable representation.  The
  // rotation is constructed by finding the face centers closest to the +x and
  // +y axes and building an orthonormal basis from them via Gram-Schmidt. This
  // yields a rotation around a tilted axis that moves vertices off all three
  // coordinate axes simultaneously.

  if (!triangles.empty()) {

    // Find the face centers (normalized to unit sphere) closest to +x and +y
    const vec3 ax = {1.0, 0.0, 0.0};
    const vec3 ay = {0.0, 1.0, 0.0};
    double best_dx = -1.0, best_dy = -1.0;
    vec3 cx = ax, cy = ay;

    for (const auto &tri : triangles) {
      vec3 c = vec3norm(tri[0] + tri[1] + tri[2]);
      double dx = vec3dot(c, ax);
      double dy = vec3dot(c, ay);
      if (dx > best_dx) { best_dx = dx; cx = c; }
      if (dy > best_dy) { best_dy = dy; cy = c; }
    }

    // Build orthonormal frame {e1, e2, e3} from the two face center directions
    // using Gram-Schmidt orthogonalization.  The resulting rotation matrix has
    // e1, e2, e3 as its rows so that e1 maps to +x, e2 maps to +y, and
    // e3 = e1 x e2 maps to +z.
    vec3 e1 = cx;
    vec3 e2 = vec3norm(cy - vec3dot(cy, e1) * e1);
    vec3 e3 = vec3cross(e1, e2);

    // clang-format off
    double R[3][3] = {{e1[0], e1[1], e1[2]},
                      {e2[0], e2[1], e2[2]},
                      {e3[0], e3[1], e3[2]}};
    // clang-format on

    // Apply rotation to all triangle vertices
    for (auto &tri : triangles) {
      for (auto &v : tri) {
        vec3 rv;
        rv[0] = R[0][0] * v[0] + R[0][1] * v[1] + R[0][2] * v[2];
        rv[1] = R[1][0] * v[0] + R[1][1] * v[1] + R[1][2] * v[2];
        rv[2] = R[2][0] * v[0] + R[2][1] * v[1] + R[2][2] * v[2];
        v = rv;
      }
    }
  }
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

  // get the offset from forward_transform to extract rotation-only for normals
  double ox = 0.0, oy = 0.0, oz = 0.0;
  reg->forward_transform(ox, oy, oz);

  // draw triangles

  const vec3 offs{center[0], center[1], center[2]};
  for (auto tri : triangles) {

    // compute ellipsoid surface normals from gradient of x^2/a^2 + y^2/b^2 + z^2/c^2
    const double sa = shape[0] * shape[0], sb = shape[1] * shape[1], sc = shape[2] * shape[2];
    vec3 n1 = vec3norm({tri[0][0] / sa, tri[0][1] / sb, tri[0][2] / sc});
    vec3 n2 = vec3norm({tri[1][0] / sa, tri[1][1] / sb, tri[1][2] / sc});
    vec3 n3 = vec3norm({tri[2][0] / sa, tri[2][1] / sb, tri[2][2] / sc});

    // set shape and move
    tri[0] = tri[0] * radscale(shape, tri[0]) + offs;
    reg->forward_transform(tri[0][0], tri[0][1], tri[0][2]);
    tri[1] = tri[1] * radscale(shape, tri[1]) + offs;
    reg->forward_transform(tri[1][0], tri[1][1], tri[1][2]);
    tri[2] = tri[2] * radscale(shape, tri[2]) + offs;
    reg->forward_transform(tri[2][0], tri[2][1], tri[2][2]);

    if (dotri) {
      // apply region rotation to normals (subtract translation offset)
      reg->forward_transform(n1[0], n1[1], n1[2]);
      reg->forward_transform(n2[0], n2[1], n2[2]);
      reg->forward_transform(n3[0], n3[1], n3[2]);
      n1 = vec3norm({n1[0] - ox, n1[1] - oy, n1[2] - oz});
      n2 = vec3norm({n2[0] - ox, n2[1] - oy, n2[2] - oz});
      n3 = vec3norm({n3[0] - ox, n3[1] - oy, n3[2] - oz});
      img->draw_trinorm(tri[0].data(), tri[1].data(), tri[2].data(), n1.data(), n2.data(),
                        n3.data(), color, color, color, opacity);
    }
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
                        const double *shape, const double *quat, double diameter, double opacity,
                        const double *block)
{
  // select between triangles or cylinders or both
  bool doframe = false;
  bool dotri = true;
  if (flag == 1) doframe = false;
  if (flag == 2) {
    dotri = false;
    doframe = true;
  }
  if (diameter <= 0.0) doframe = false;
  if (!dotri && !doframe) return;    // nothing to do

  double p[3][3];
  vec3 e1, e2, e3;
  const vec3 offs{center[0], center[1], center[2]};

  // optimization: just draw a sphere if a filled surface is requested and the object is a sphere
  // note: this does not apply to superellipsoids
  if (dotri && !block && (shape[0] == shape[1]) && (shape[0] == shape[2])) {
    img->draw_sphere(center, color, 2.0 * shape[0], opacity);
    return;
  }

  // nothing to draw
  if (!triangles.size()) return;

  // get rotation matrix for body frame to box frame
  MathExtra::quat_to_mat(quat, p);

  // draw triangles and edges as requested, work on copy of triangle since we modify it
  for (auto tri : triangles) {

    // compute surface normals from unit sphere coordinates before scaling
    vec3 n1, n2, n3;
    if (dotri) {
      if (block) {
        // compute superellipsoid surface normals from gradient of implicit function
        n1 = supernormal(shape, block, tri[0]);
        n2 = supernormal(shape, block, tri[1]);
        n3 = supernormal(shape, block, tri[2]);
      } else {
        // compute ellipsoid surface normals from gradient of x^2/a^2 + y^2/b^2 + z^2/c^2
        const double sa = shape[0] * shape[0], sb = shape[1] * shape[1], sc = shape[2] * shape[2];
        n1 = vec3norm({tri[0][0] / sa, tri[0][1] / sb, tri[0][2] / sc});
        n2 = vec3norm({tri[1][0] / sa, tri[1][1] / sb, tri[1][2] / sc});
        n3 = vec3norm({tri[2][0] / sa, tri[2][1] / sb, tri[2][2] / sc});
      }
    }

    // set shape by shifting each corner to the surface
    if (block) {
      for (int i = 0; i < 3; ++i) {
        auto &t = tri[i];
        t = superscale(shape, block, t) * t;
      }
    } else {
      for (int i = 0; i < 3; ++i) {
        auto &t = tri[i];
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

    if (dotri) {
      // rotate normals (no translation or scaling)
      vec3 rn1, rn2, rn3;
      MathExtra::matvec(p, n1.data(), rn1.data());
      MathExtra::matvec(p, n2.data(), rn2.data());
      MathExtra::matvec(p, n3.data(), rn3.data());

      img->draw_trinorm(e1.data(), e2.data(), e3.data(), rn1.data(), rn2.data(), rn3.data(), color,
                        color, color, opacity);
    }

    if (doframe) {
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
