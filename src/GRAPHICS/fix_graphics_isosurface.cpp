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

#include "fix_graphics_isosurface.h"

#include "arg_info.h"
#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "graphics.h"
#include "input.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "respa.h"
#include "safe_pointers.h"
#include "update.h"
#include "variable.h"

#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>

using namespace LAMMPS_NS;

namespace {
constexpr double BIG = 1.0e200;

// for choosing the mesh resolution
enum { NONE, SURFMIN, SURFLOW, SURFMED, SURFHIGH, SURFMAX };
constexpr double GRIDPOINTS[] = {0, 8.0, 16.0, 64.0, 128.0, 256.0};

// for choosing property
enum { NUMBER, MASS, CHARGE, VARIABLE, COMPUTE, FIX, CUSTOM };

// for truncating gaussian spreading
constexpr double CUTVAL = 9.210340371976182;    // = -std::log(0.0001);
// extra grid points at the sides of the subbox grid
constexpr int GRIDEXTRA = 4;

// custom data types for positions and triangles based on std::array
using vec3 = std::array<double, 3>;
using triangle = struct {
  std::array<vec3, 3> triangle;
  std::array<vec3, 3> normals;
  std::array<int, 3> type;
};
using gridcell = struct {
  vec3 pos[8];
  double iso[8];
  vec3 grad[8];
  int type[8];
};

inline vec3 operator-(const vec3 &a, const vec3 &b)
{
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

// compute the gradient at a grid point using central differences (clamped at boundaries)
void compute_gradient(double ***isogrid, int cx, int cy, int cz, int nx, int ny, int nz, vec3 &grad)
{
  int xm = (cx > 0) ? cx - 1 : 0;
  int xp = (cx < nx - 1) ? cx + 1 : nx - 1;
  int ym = (cy > 0) ? cy - 1 : 0;
  int yp = (cy < ny - 1) ? cy + 1 : ny - 1;
  int zm = (cz > 0) ? cz - 1 : 0;
  int zp = (cz < nz - 1) ? cz + 1 : nz - 1;
  grad[0] = isogrid[xp][cy][cz] - isogrid[xm][cy][cz];
  grad[1] = isogrid[cx][yp][cz] - isogrid[cx][ym][cz];
  grad[2] = isogrid[cx][cy][zp] - isogrid[cx][cy][zm];
}

// get vertex position and interpolated normal for a grid cell edge
// by interpolating between the two corners based on their iso values
void get_vertex(const gridcell &g, int c1, int c2, vec3 &vert, vec3 &norm, int &type)
{
  const double diff = g.iso[c2] - g.iso[c1];
  const double fraction = (fabs(diff) > 0.0) ? -g.iso[c1] / diff : 0.0;

  // assign type closest to vertex, if available
  if (g.type[c1] && !g.type[c2]) {
    type = g.type[c1];
  } else if (!g.type[c1] && g.type[c2]) {
    type = g.type[c2];
  } else if (g.type[c1] && g.type[c2]) {
    type = (fraction < 0.5) ? g.type[c1] : g.type[c2];
  } else {
    type = 0;
  }

  for (int d = 0; d < 3; ++d) {
    vert[d] = g.pos[c1][d] + fraction * (g.pos[c2][d] - g.pos[c1][d]);
    norm[d] = g.grad[c1][d] + fraction * (g.grad[c2][d] - g.grad[c1][d]);
  }
  double len = sqrt(norm[0] * norm[0] + norm[1] * norm[1] + norm[2] * norm[2]);
  if (len > 0.0) {
    norm[0] /= len;
    norm[1] /= len;
    norm[2] /= len;
  }
}

// spread out atom data across the grid
void distribute(double ***pgrid, double ***dgrid, int ***tgrid, const double *pos, int type,
                double val, const double *sublo, double delta, int nx, int ny, int nz, int nrange,
                double rcutsq, double sigma)
{
  // locate the primary grid cell of the atom
  const int ix = GRIDEXTRA + static_cast<int>(floor((pos[0] - sublo[0]) / delta));
  const int iy = GRIDEXTRA + static_cast<int>(floor((pos[1] - sublo[1]) / delta));
  const int iz = GRIDEXTRA + static_cast<int>(floor((pos[2] - sublo[2]) / delta));

  // compute relative position inside the cell
  const double dx = (pos[0] - (sublo[0] + (ix - GRIDEXTRA) * delta));
  const double dy = (pos[1] - (sublo[1] + (iy - GRIDEXTRA) * delta));
  const double dz = (pos[2] - (sublo[2] + (iz - GRIDEXTRA) * delta));

  // loop over possible grid positions for spreading the data
#if defined(_OPENMP)
#pragma omp parallel for
#endif
  for (int jx = -nrange; jx <= nrange; ++jx) {
    int kx = ix + jx;
    // skip if outside the local grid
    if ((kx < 0) || (kx >= nx)) continue;
    for (int jy = -nrange; jy <= nrange; ++jy) {
      int ky = iy + jy;
      // skip if outside the local grid
      if ((ky < 0) || (ky >= ny)) continue;
      for (int jz = -nrange; jz <= nrange; ++jz) {
        int kz = iz + jz;
        // skip if outside the local grid
        if ((kz < 0) || (kz >= nz)) continue;

        // compute squared distance of grid point to atom and apply cutoff and set value
        double xpos = jx * delta - dx;
        double ypos = jy * delta - dy;
        double zpos = jz * delta - dz;
        double distsq = xpos * xpos + ypos * ypos + zpos * zpos;
        if (distsq < rcutsq) {
          pgrid[kx][ky][kz] += val * exp(-distsq / sigma);
          if (distsq < dgrid[kx][ky][kz]) {
            dgrid[kx][ky][kz] = distsq;
            tgrid[kx][ky][kz] = type;
          }
        }
      }
    }
  }
}

// for identifying the edges that contain the isosurface
constexpr int EDGETABLE[256] = {
    0x0,   0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c, 0x80c, 0x905, 0xa0f, 0xb06, 0xc0a,
    0xd03, 0xe09, 0xf00, 0x190, 0x99,  0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c, 0x99c, 0x895,
    0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90, 0x230, 0x339, 0x33,  0x13a, 0x636, 0x73f, 0x435,
    0x53c, 0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30, 0x3a0, 0x2a9, 0x1a3, 0xaa,
    0x7a6, 0x6af, 0x5a5, 0x4ac, 0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0, 0x460,
    0x569, 0x663, 0x76a, 0x66,  0x16f, 0x265, 0x36c, 0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963,
    0xa69, 0xb60, 0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff,  0x3f5, 0x2fc, 0xdfc, 0xcf5, 0xfff,
    0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0, 0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55,  0x15c,
    0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950, 0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6,
    0x2cf, 0x1c5, 0xcc,  0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0, 0x8c0, 0x9c9,
    0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc, 0xcc,  0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9,
    0x7c0, 0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c, 0x15c, 0x55,  0x35f, 0x256,
    0x55a, 0x453, 0x759, 0x650, 0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc, 0x2fc,
    0x3f5, 0xff,  0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0, 0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f,
    0xd65, 0xc6c, 0x36c, 0x265, 0x16f, 0x66,  0x76a, 0x663, 0x569, 0x460, 0xca0, 0xda9, 0xea3,
    0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac, 0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa,  0x1a3, 0x2a9, 0x3a0,
    0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c, 0x53c, 0x435, 0x73f, 0x636, 0x13a,
    0x33,  0x339, 0x230, 0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c, 0x69c, 0x795,
    0x49f, 0x596, 0x29a, 0x393, 0x99,  0x190, 0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905,
    0x80c, 0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0};

// to map the "vertex" index from above to a list of triangles
// of the interpolated edge vertices.
constexpr int TRITABLE[256][16] = {
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1},
    {3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1},
    {3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1},
    {3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1},
    {9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1},
    {9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
    {2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1},
    {8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1},
    {9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
    {4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1},
    {3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1},
    {1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1},
    {4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1},
    {4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
    {5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1},
    {2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1},
    {9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1},
    {0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1},
    {2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1},
    {10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1},
    {5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1},
    {5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1},
    {9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1},
    {0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1},
    {1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1},
    {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1},
    {8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1},
    {2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1},
    {7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1},
    {2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1},
    {11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1},
    {5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1},
    {11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1},
    {11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
    {1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1},
    {9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1},
    {5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1},
    {2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1},
    {5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1},
    {6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1},
    {3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1},
    {6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1},
    {5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1},
    {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1},
    {6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1},
    {8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1},
    {7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1},
    {3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1},
    {5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1},
    {0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1},
    {9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1},
    {8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1},
    {5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1},
    {0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1},
    {6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1},
    {10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1},
    {10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1},
    {8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1},
    {1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1},
    {0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1},
    {10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1},
    {3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1},
    {6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1},
    {9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1},
    {8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1},
    {3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1},
    {6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1},
    {0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1},
    {10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1},
    {10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1},
    {2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1},
    {7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1},
    {7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1},
    {2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1},
    {1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1},
    {11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1},
    {8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1},
    {0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1},
    {7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
    {10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
    {2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1},
    {6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1},
    {7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1},
    {2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1},
    {1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1},
    {10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1},
    {10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1},
    {0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1},
    {7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1},
    {6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1},
    {8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1},
    {9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1},
    {6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1},
    {4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1},
    {10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1},
    {8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1},
    {0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1},
    {1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1},
    {8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1},
    {10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1},
    {4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1},
    {10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1},
    {5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
    {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1},
    {9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1},
    {6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1},
    {7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1},
    {3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1},
    {7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1},
    {9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1},
    {3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1},
    {6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1},
    {9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1},
    {1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1},
    {4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1},
    {7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1},
    {6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1},
    {3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1},
    {0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1},
    {6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1},
    {0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1},
    {11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1},
    {6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1},
    {5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1},
    {9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1},
    {1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1},
    {1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1},
    {10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1},
    {0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1},
    {5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1},
    {10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1},
    {11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1},
    {9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1},
    {7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1},
    {2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1},
    {8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1},
    {9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1},
    {9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1},
    {1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1},
    {9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1},
    {9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1},
    {5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1},
    {0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1},
    {10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1},
    {2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1},
    {0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1},
    {0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1},
    {9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1},
    {5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1},
    {3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1},
    {5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1},
    {8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1},
    {0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1},
    {9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1},
    {0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1},
    {1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1},
    {3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1},
    {4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1},
    {9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1},
    {11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1},
    {11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1},
    {2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1},
    {9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1},
    {3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1},
    {1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1},
    {4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1},
    {4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1},
    {0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1},
    {3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1},
    {3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1},
    {0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1},
    {9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1},
    {1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}};
}    // namespace

/* ---------------------------------------------------------------------- */

FixGraphicsIsosurface::FixGraphicsIsosurface(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), pdata(nullptr), pstr(nullptr), pcomp(nullptr), pfix(nullptr),
    imgobjs(nullptr), imgparms(nullptr)
{
  if (narg < 6) utils::missing_cmd_args(FLERR, "fix graphics/isosurface", error);

  // fix settings
  global_freq = nevery;
  dynamic_group_allow = 1;
  comm_forward = 1;
  nlevels_respa = 1;
  nmax = -1;

  if (domain->triclinic)
    error->all(FLERR, "Fix graphics/isosurface is currently not compatible with triclinic cells");
  if (domain->dimension == 2)
    error->all(FLERR, "Fix graphics/isosurface is currently not compatible with 2d systems");

  // defaults
  numobjs = 0;
  quality = SURFLOW;
  pflag = NUMBER;
  binary = 0;
  pad = 0;
  pvar = -1;
  pindex = -1;

  // parse mandatory args

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  if (nevery <= 0) error->all(FLERR, 3, "Illegal fix graphics/isosurface nevery value {}", nevery);
  iso = utils::numeric(FLERR, arg[4], false, lmp);
  rad = 2.0 * utils::numeric(FLERR, arg[5], false, lmp);
  if (rad <= 0.0) error->all(FLERR, 5, "Illegal fix graphics/isosurface radius value {}", rad);

  // parse optional args

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "quality") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix graphics/isosurface quality", error);
      if (strcmp(arg[iarg + 1], "min") == 0) {
        quality = SURFMIN;
      } else if (strcmp(arg[iarg + 1], "low") == 0) {
        quality = SURFLOW;
      } else if (strcmp(arg[iarg + 1], "med") == 0) {
        quality = SURFMED;
      } else if (strcmp(arg[iarg + 1], "high") == 0) {
        quality = SURFHIGH;
      } else if (strcmp(arg[iarg + 1], "max") == 0) {
        quality = SURFMAX;
      } else {
        error->all(FLERR, iarg + 1, "Unknown fix graphics/isosurface quality setting {}",
                   arg[iarg + 1]);
      }
      iarg += 2;
    } else if (strcmp(arg[iarg], "property") == 0) {
      if (iarg + 2 > narg)
        utils::missing_cmd_args(FLERR, "fix graphics/isosurface property", error);
      if (strcmp(arg[iarg + 1], "none") == 0) {
        pflag = ArgInfo::NONE;
      } else if (strcmp(arg[iarg + 1], "mass") == 0) {
        pflag = ArgInfo::MASS;
      } else {
        ArgInfo argi(arg[iarg + 1]);
        pflag = argi.get_type();
        pindex = argi.get_index1();
        delete[] pstr;
        pstr = utils::strdup(argi.get_name());

        // check validity of property argument
        if ((pflag == ArgInfo::UNKNOWN) || (pflag == ArgInfo::NONE) || (argi.get_dim() > 1))
          error->all(FLERR, iarg + 1, "Invalid fix graphics/isosurface property {}", arg[iarg + 1]);

        if (pflag == ArgInfo::COMPUTE) {
          pcomp = modify->get_compute_by_id(pstr);
          if (!pcomp)
            error->all(FLERR, iarg + 1, "Compute ID {} for fix graphics/isosurface does not exist",
                       pstr);
          if (pcomp->peratom_flag == 0)
            error->all(FLERR, iarg + 1,
                       "Fix graphics/isosurface compute {} does not calculate per-atom values",
                       pstr);
          if (pindex == 0 && pcomp->size_peratom_cols != 0)
            error->all(FLERR, iarg + 1,
                       "Fix graphics/isosurface compute {} does not calculate a per-atom vector",
                       pstr);
          if (pindex && pcomp->size_peratom_cols == 0)
            error->all(FLERR, iarg + 1,
                       "Fix graphics/isosurface compute {} does not calculate a per-atom array",
                       pstr);
          if (pindex && pindex > pcomp->size_peratom_cols)
            error->all(FLERR, iarg + 1,
                       "Fix graphics/isosurface compute {} array is accessed out-of-range{}", pstr,
                       utils::errorurl(20));

        } else if (pflag == ArgInfo::FIX) {
          pfix = modify->get_fix_by_id(pstr);
          if (!pfix)
            error->all(FLERR, iarg + 1, "Fix ID {} for fix graphics/isosurface does not exist",
                       pstr);
          if (pfix->peratom_flag == 0)
            error->all(FLERR, iarg + 1,
                       "Fix graphics/isosurface fix {} does not calculate per-atom values", pstr);
          if (pindex == 0 && pfix->size_peratom_cols != 0)
            error->all(FLERR, iarg + 1,
                       "Fix graphics/isosurface fix {} does not calculate a per-atom vector", pstr);
          if (pindex && pfix->size_peratom_cols == 0)
            error->all(FLERR, iarg + 1,
                       "Fix graphics/isosurface fix {} does not calculate a per-atom array", pstr);
          if (pindex && pindex > pfix->size_peratom_cols)
            error->all(FLERR, iarg + 1,
                       "Fix graphics/isosurface fix {} array is accessed out-of-range{}", pstr,
                       utils::errorurl(20));
          if (nevery % pfix->peratom_freq)
            error->all(FLERR, iarg + 1,
                       "Fix {} for fix graphics/isosurface not computed at compatible time{}", pstr,
                       utils::errorurl(7));

        } else if (pflag == ArgInfo::VARIABLE) {
          pvar = input->variable->find(pstr);
          if (pvar < 0)
            error->all(FLERR, iarg + 1,
                       "Variable name {} for fix graphics/isosurface does not exist", pstr);
          if (input->variable->atomstyle(pvar) == 0)
            error->all(FLERR, iarg + 1,
                       "Fix graphics/isosurface variable {} is not atom-style variable", pstr);
        }
      }
      iarg += 2;
    } else if (strcmp(arg[iarg], "filename") == 0) {
      if (iarg + 2 > narg)
        utils::missing_cmd_args(FLERR, "fix graphics/isosurface filename", error);
      filename = arg[iarg + 1];
      if (filename.find('*') == std::string::npos)
        error->all(FLERR, iarg + 1, "Output to STL file requires file name with '*'");
      iarg += 2;
    } else if (strcmp(arg[iarg], "binary") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix graphics/isosurface binary", error);
      binary = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "pad") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix graphics/isosurface pad", error);
      pad = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else {
      error->all(FLERR, iarg, "Unknown fix graphics/isosurface keyword {}", arg[iarg]);
    }
  }
}

/* ---------------------------------------------------------------------- */

FixGraphicsIsosurface::~FixGraphicsIsosurface()
{
  delete[] pstr;
  memory->destroy(pdata);
  memory->destroy(imgobjs);
  memory->destroy(imgparms);
}

/* ---------------------------------------------------------------------- */

int FixGraphicsIsosurface::setmask()
{
  int mask = 0;
  mask |= FixConst::POST_FORCE;
  mask |= FixConst::POST_FORCE_RESPA;
  mask |= FixConst::END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixGraphicsIsosurface::init()
{
  if (utils::strmatch(update->integrate_style, "^respa"))
    nlevels_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels;

  // check validity of all compute, fix, or variable and update

  if (pflag == ArgInfo::COMPUTE) {
    pcomp = modify->get_compute_by_id(pstr);
    if (!pcomp)
      error->all(FLERR, Error::NOLASTLINE, "Compute ID {} for fix ave/atom does not exist", pstr);

  } else if (pflag == ArgInfo::FIX) {
    pfix = modify->get_fix_by_id(pstr);
    if (!pfix)
      error->all(FLERR, Error::NOLASTLINE, "Fix ID {} for fix ave/atom does not exist", pstr);

  } else if (pflag == ArgInfo::VARIABLE) {
    pvar = input->variable->find(pstr);
    if (pvar < 0)
      error->all(FLERR, Error::NOLASTLINE, "Variable name {} for fix ave/atom does not exist",
                 pstr);
  }

  // request updates of computes, if needed
  if ((pflag != ArgInfo::NONE) && (pflag != ArgInfo::MASS)) {
    bigint nextstep = (update->ntimestep / nevery) * nevery + nevery;
    if ((nextstep - nevery) == update->ntimestep) nextstep = update->ntimestep;
    modify->addstep_compute(nextstep);
  }
}

/* ---------------------------------------------------------------------- */

void FixGraphicsIsosurface::setup(int vflag)
{
  post_force(vflag);
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixGraphicsIsosurface::post_force(int /*vflag*/)
{
  // manage updates of computes, if needed
  if ((pflag != ArgInfo::NONE) && (pflag != ArgInfo::MASS)) modify->clearstep_compute();

  // ensure data storage is sufficient
  if (nmax < atom->nmax) {
    nmax = MAX(1, atom->nmax);
    memory->destroy(pdata);
    memory->create(pdata, nmax, "fix_graphics/isosurface:pdata");
  }

  // fill per-atom data storage with requested data for local atoms
  const int nlocal = atom->nlocal;
  const int *const type = atom->type;
  const double *const mass = atom->mass;
  const double *const rmass = atom->rmass;

  if (pflag == ArgInfo::NONE) {
    for (int i = 0; i < nlocal; ++i) pdata[i] = 1.0;
  } else if (pflag == ArgInfo::MASS) {
    if (rmass) {
      for (int i = 0; i < nlocal; ++i) pdata[i] = rmass[i];
    } else {
      for (int i = 0; i < nlocal; ++i) pdata[i] = mass[type[i]];
    }
  } else if (pflag == ArgInfo::COMPUTE) {
    if (!(pcomp->invoked_flag & Compute::INVOKED_PERATOM)) {
      pcomp->compute_peratom();
      pcomp->invoked_flag |= Compute::INVOKED_PERATOM;
    }
    if (pindex == 0) {
      double *compute_vector = pcomp->vector_atom;
      for (int i = 0; i < nlocal; ++i) pdata[i] = compute_vector[i];
    } else {
      int jm1 = pindex - 1;
      double **compute_array = pcomp->array_atom;
      for (int i = 0; i < nlocal; ++i) pdata[i] = compute_array[i][jm1];
    }
  } else if (pflag == ArgInfo::FIX) {
    if (pindex == 0) {
      double *fix_vector = pfix->vector_atom;
      for (int i = 0; i < nlocal; ++i) pdata[i] = fix_vector[i];
    } else {
      int jm1 = pindex - 1;
      double **fix_array = pfix->array_atom;
      for (int i = 0; i < nlocal; ++i) pdata[i] = fix_array[i][jm1];
    }
  } else if (pflag == ArgInfo::VARIABLE) {
    input->variable->compute_atom(pvar, igroup, pdata, 1, 0);
  }

  // set data for ghost atoms
  comm->forward_comm(this);

  // request updates of computes, if needed
  if ((pflag != ArgInfo::NONE) && (pflag != ArgInfo::MASS)) {
    bigint nextstep = (update->ntimestep / nevery) * nevery + nevery;
    if ((nextstep - nevery) == update->ntimestep) nextstep = update->ntimestep;
    modify->addstep_compute(nextstep);
  }
}

/* ---------------------------------------------------------------------- */

int FixGraphicsIsosurface::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/,
                                             int * /*pbc*/)
{
  int m = 0;
  for (int i = 0; i < n; ++i) buf[m++] = pdata[list[i]];

  return m;
}

/* ---------------------------------------------------------------------- */

void FixGraphicsIsosurface::unpack_forward_comm(int n, int first, double *buf)
{
  int m = 0;
  int last = first + n;

  for (int i = first; i < last; ++i) pdata[i] = buf[m++];
}

/* ---------------------------------------------------------------------- */

void FixGraphicsIsosurface::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa - 1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixGraphicsIsosurface::end_of_step()
{
  // determine grid dimensions

  const double *const sublo = domain->sublo;
  const double *const subhi = domain->subhi;

  // get grid spacing and sub-domain lengths
  double delta = cbrt(domain->xprd * domain->yprd * domain->zprd) / GRIDPOINTS[quality];
  double sublen[3] = {subhi[0] - sublo[0], subhi[1] - sublo[1], subhi[2] - sublo[2]};

  // get grid dims for subdomains, add extra points on each side
  int nx = 2 * GRIDEXTRA + static_cast<int>(ceil(sublen[0] / delta));
  int ny = 2 * GRIDEXTRA + static_cast<int>(ceil(sublen[1] / delta));
  int nz = 2 * GRIDEXTRA + static_cast<int>(ceil(sublen[2] / delta));

  // determine cutoff for spreading and number of gridpoints for spreading
  const double rcutsq = rad * CUTVAL;
  const int nrange = static_cast<int>(ceil(sqrt(rcutsq) / delta));

  // allocate grids and zero it out the property grid
  double ***isogrid, ***distgrid;
  int ***typegrid;
  memory->create(isogrid, nx, ny, nz, "fix graphics/isosurface:isogrid");
  memory->create(distgrid, nx, ny, nz, "fix graphics/isosurface:distgrid");
  memory->create(typegrid, nx, ny, nz, "fix graphics/isosurface:typegrid");
  memset(isogrid[0][0], 0, sizeof(double) * nx * ny * nz);
  memset(typegrid[0][0], 0, sizeof(int) * nx * ny * nz);
  for (int ix = 0; ix < nx; ++ix) {
    for (int iy = 0; iy < ny; ++iy) {
      for (int iz = 0; iz < nz; ++iz) distgrid[ix][iy][iz] = BIG;
    }
  }

  // loop over local and ghost atoms, match group and check if within grid range
  // per-atom data was copied to pdata array and sent to ghost atoms in post_force()
  const int nall = atom->nlocal + atom->nghost;
  const int *const mask = atom->mask;
  const int *const type = atom->type;
  const double *const *const x = atom->x;
  for (int i = 0; i < nall; ++i) {
    if (mask[i] & groupbit)
      distribute(isogrid, distgrid, typegrid, x[i], type[i], pdata[i], sublo, delta, nx, ny, nz,
                 nrange, rcutsq, rad);
  }

  // subtract the isovalue from grid data so that the isosurface would be drawn for a value of < 0.0
  for (int ix = 0; ix < nx; ++ix) {
    for (int iy = 0; iy < ny; ++iy) {
      for (int iz = 0; iz < nz; ++iz) isogrid[ix][iy][iz] -= iso;
    }
  }

  // now construct list of triangles for the isosurface using the marching cubes algorithm
  // evaluate variable if necessary, wrap with clear/add

  gridcell g;
  std::vector<triangle> triangles;

  for (int ix = 0; ix < nx - 1; ++ix) {
    for (int iy = 0; iy < ny - 1; ++iy) {
      for (int iz = 0; iz < nz - 1; ++iz) {
        // clang-format off

        // get lower edge coordinates of the grid cell
        double gx = sublo[0] + (ix - GRIDEXTRA) * delta;
        double gy = sublo[1] + (iy - GRIDEXTRA) * delta;
        double gz = sublo[2] + (iz - GRIDEXTRA) * delta;

        // store position and isovalue for each corner of the grid cell
        g.iso[0] = isogrid[ix  ][iy  ][iz  ];
        g.iso[1] = isogrid[ix+1][iy  ][iz  ];
        g.iso[3] = isogrid[ix  ][iy+1][iz  ];
        g.iso[2] = isogrid[ix+1][iy+1][iz  ];
        g.iso[4] = isogrid[ix  ][iy  ][iz+1];
        g.iso[5] = isogrid[ix+1][iy  ][iz+1];
        g.iso[7] = isogrid[ix  ][iy+1][iz+1];
        g.iso[6] = isogrid[ix+1][iy+1][iz+1];
        g.pos[0] = {gx,       gy,       gz};
        g.pos[1] = {gx+delta, gy,        gz};
        g.pos[3] = {gx,       gy+delta, gz};
        g.pos[2] = {gx+delta, gy+delta, gz};
        g.pos[4] = {gx,       gy,       gz+delta};
        g.pos[5] = {gx+delta, gy,       gz+delta};
        g.pos[7] = {gx,       gy+delta, gz+delta};
        g.pos[6] = {gx+delta, gy+delta, gz+delta};
        g.type[0] = typegrid[ix  ][iy  ][iz  ];
        g.type[1] = typegrid[ix+1][iy  ][iz  ];
        g.type[3] = typegrid[ix  ][iy+1][iz  ];
        g.type[2] = typegrid[ix+1][iy+1][iz  ];
        g.type[4] = typegrid[ix  ][iy  ][iz+1];
        g.type[5] = typegrid[ix+1][iy  ][iz+1];
        g.type[7] = typegrid[ix  ][iy+1][iz+1];
        g.type[6] = typegrid[ix+1][iy+1][iz+1];
        // clang-format on

        // compute gradient at each corner for smooth normals
        compute_gradient(isogrid, ix, iy, iz, nx, ny, nz, g.grad[0]);
        compute_gradient(isogrid, ix + 1, iy, iz, nx, ny, nz, g.grad[1]);
        compute_gradient(isogrid, ix + 1, iy + 1, iz, nx, ny, nz, g.grad[2]);
        compute_gradient(isogrid, ix, iy + 1, iz, nx, ny, nz, g.grad[3]);
        compute_gradient(isogrid, ix, iy, iz + 1, nx, ny, nz, g.grad[4]);
        compute_gradient(isogrid, ix + 1, iy, iz + 1, nx, ny, nz, g.grad[5]);
        compute_gradient(isogrid, ix + 1, iy + 1, iz + 1, nx, ny, nz, g.grad[6]);
        compute_gradient(isogrid, ix, iy + 1, iz + 1, nx, ny, nz, g.grad[7]);

        // determine edge table index
        int idx = 0;
        if (g.iso[0] < 0.0) idx |= 1;
        if (g.iso[1] < 0.0) idx |= 2;
        if (g.iso[2] < 0.0) idx |= 4;
        if (g.iso[3] < 0.0) idx |= 8;
        if (g.iso[4] < 0.0) idx |= 16;
        if (g.iso[5] < 0.0) idx |= 32;
        if (g.iso[6] < 0.0) idx |= 64;
        if (g.iso[7] < 0.0) idx |= 128;

        // gridcube is not crossed by isosurface
        if (EDGETABLE[idx] == 0) continue;

        // compute the possible 12 triangle vertices and normals
        std::array<vec3, 12> vertices;
        std::array<vec3, 12> vnormals;
        std::array<int, 12> vtypes;
        if (EDGETABLE[idx] & 1) get_vertex(g, 0, 1, vertices[0], vnormals[0], vtypes[0]);
        if (EDGETABLE[idx] & 2) get_vertex(g, 1, 2, vertices[1], vnormals[1], vtypes[1]);
        if (EDGETABLE[idx] & 4) get_vertex(g, 2, 3, vertices[2], vnormals[2], vtypes[2]);
        if (EDGETABLE[idx] & 8) get_vertex(g, 3, 0, vertices[3], vnormals[3], vtypes[3]);
        if (EDGETABLE[idx] & 16) get_vertex(g, 4, 5, vertices[4], vnormals[4], vtypes[4]);
        if (EDGETABLE[idx] & 32) get_vertex(g, 5, 6, vertices[5], vnormals[5], vtypes[5]);
        if (EDGETABLE[idx] & 64) get_vertex(g, 6, 7, vertices[6], vnormals[6], vtypes[6]);
        if (EDGETABLE[idx] & 128) get_vertex(g, 7, 4, vertices[7], vnormals[7], vtypes[7]);
        if (EDGETABLE[idx] & 256) get_vertex(g, 0, 4, vertices[8], vnormals[8], vtypes[8]);
        if (EDGETABLE[idx] & 512) get_vertex(g, 1, 5, vertices[9], vnormals[9], vtypes[9]);
        if (EDGETABLE[idx] & 1024) get_vertex(g, 2, 6, vertices[10], vnormals[10], vtypes[10]);
        if (EDGETABLE[idx] & 2048) get_vertex(g, 3, 7, vertices[11], vnormals[11], vtypes[11]);

        // compute the triangles for this grid cell and add them to the list
        for (int i = 0; TRITABLE[idx][i] != -1; i += 3)
          triangles.emplace_back(
              triangle{{vertices[TRITABLE[idx][i]], vertices[TRITABLE[idx][i + 1]],
                        vertices[TRITABLE[idx][i + 2]]},
                       {vnormals[TRITABLE[idx][i]], vnormals[TRITABLE[idx][i + 1]],
                        vnormals[TRITABLE[idx][i + 2]]},
                       {vtypes[TRITABLE[idx][i]], vtypes[TRITABLE[idx][i + 1]],
                        vtypes[TRITABLE[idx][i + 2]]}});
      }
    }
  }
  // we don't need the iso value grid anymore.
  memory->destroy(isogrid);

  // allocate and assign list of graphics objects
  numobjs = triangles.size();
  memory->destroy(imgobjs);
  memory->destroy(imgparms);
  memory->create(imgobjs, numobjs, "fix_graphics:imgobjs");
  memory->create(imgparms, numobjs, 21, "fix_graphics:imgparms");

  int n = 0;
  for (const auto &tri : triangles) {
    // skip if any part of the triangle is outside the subdomain
    bool addme = true;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        if ((tri.triangle[i][j] < sublo[j]) || (tri.triangle[i][j] > subhi[j])) addme = false;
    if (addme) {
      imgobjs[n] = Graphics::TRINORM;
      imgparms[n][0] = tri.type[0];
      imgparms[n][1] = tri.type[1];
      imgparms[n][2] = tri.type[2];
      imgparms[n][3] = tri.triangle[0][0];
      imgparms[n][4] = tri.triangle[0][1];
      imgparms[n][5] = tri.triangle[0][2];
      imgparms[n][6] = tri.triangle[1][0];
      imgparms[n][7] = tri.triangle[1][1];
      imgparms[n][8] = tri.triangle[1][2];
      imgparms[n][9] = tri.triangle[2][0];
      imgparms[n][10] = tri.triangle[2][1];
      imgparms[n][11] = tri.triangle[2][2];
      imgparms[n][12] = tri.normals[0][0];
      imgparms[n][13] = tri.normals[0][1];
      imgparms[n][14] = tri.normals[0][2];
      imgparms[n][15] = tri.normals[1][0];
      imgparms[n][16] = tri.normals[1][1];
      imgparms[n][17] = tri.normals[1][2];
      imgparms[n][18] = tri.normals[2][0];
      imgparms[n][19] = tri.normals[2][1];
      imgparms[n][20] = tri.normals[2][2];
      ++n;
    }
  }
  numobjs = n;

  // write grid to STL format file, if requested
  if (filename.size() > 0) {
    int maxobjs = 0;
    uint32_t allobjs = 0U;
    MPI_Allreduce(&numobjs, &maxobjs, 1, MPI_INT, MPI_SUM, world);
    allobjs = maxobjs;
    MPI_Allreduce(&numobjs, &maxobjs, 1, MPI_INT, MPI_MAX, world);

    // convert to single precision data (for binary output) and compute normals
    auto *stldata = new float[12 * maxobjs];    // 3 floats for normal and 9 floats for corners
    double normal[3];
    vec3 d1, d2;
    int n = 0;
    for (const auto &tri : triangles) {
      // skip if any part of the triangle is outside the subdomain
      bool addme = true;
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          if ((tri.triangle[i][j] < sublo[j]) || (tri.triangle[i][j] > subhi[j])) addme = false;
      if (addme) {
        d1 = tri.triangle[1] - tri.triangle[0];
        d2 = tri.triangle[2] - tri.triangle[0];
        MathExtra::cross3(d1.data(), d2.data(), normal);
        MathExtra::norm3(normal);
        stldata[12 * n] = normal[0];
        stldata[12 * n + 1] = normal[1];
        stldata[12 * n + 2] = normal[2];
        stldata[12 * n + 3] = tri.triangle[0][0];
        stldata[12 * n + 4] = tri.triangle[0][1];
        stldata[12 * n + 5] = tri.triangle[0][2];
        stldata[12 * n + 6] = tri.triangle[1][0];
        stldata[12 * n + 7] = tri.triangle[1][1];
        stldata[12 * n + 8] = tri.triangle[1][2];
        stldata[12 * n + 9] = tri.triangle[2][0];
        stldata[12 * n + 10] = tri.triangle[2][1];
        stldata[12 * n + 11] = tri.triangle[2][2];
        ++n;
      }
    }

    SafeFilePtr fp;
    if (comm->me == 0) {    // only MPI rank 0 writes to the file
      auto *filecurrent = utils::strdup(utils::star_subst(filename, update->ntimestep, pad));
      if (platform::has_compress_extension(filename)) {
        if (binary)
          error->one(FLERR, Error::NOLASTLINE, "Connot use compression with binary output: {}",
                     filename);
        fp.set_pclose();
        fp = platform::compressed_write(filecurrent);
      } else if (binary) {
        fp = fopen(filecurrent, "wb");
      } else {
        fp = fopen(filecurrent, "w");
      }
      if (fp == nullptr)
        error->one(FLERR, Error::NOLASTLINE, "Cannot open STL output file {} for writing: {}",
                   filecurrent, utils::getsyserror());

      auto title = fmt::format("STL isosurface from fix {} graphics/isosurface on step {}", id,
                               update->ntimestep);
      title.resize(80, '\0');

      // write out data from rank 0
      if (binary) {
        uint16_t attributes = 0;
        fwrite(title.c_str(), 1, 80, fp);
        fwrite(&allobjs, sizeof(uint32_t), 1, fp);
        for (int i = 0; i < numobjs; ++i) {
          fwrite(stldata + 12 * i, sizeof(float), 12, fp);
          fwrite(&attributes, sizeof(uint16_t), 1, fp);
        }
      } else {
        fprintf(fp, "solid %s\n", title.c_str());
        for (int i = 0; i < numobjs; ++i) {
          utils::print(fp, "  facet normal {:e} {:e} {:e}\n", stldata[12 * i], stldata[12 * i + 1],
                       stldata[12 * i + 2]);
          fputs("    outer loop\n", fp);
          utils::print(fp, "      vertex {:e} {:e} {:e}\n", stldata[12 * i + 3],
                       stldata[12 * i + 4], stldata[12 * i + 5]);
          utils::print(fp, "      vertex {:e} {:e} {:e}\n", stldata[12 * i + 6],
                       stldata[12 * i + 7], stldata[12 * i + 8]);
          utils::print(fp, "      vertex {:e} {:e} {:e}\n", stldata[12 * i + 9],
                       stldata[12 * i + 10], stldata[12 * i + 11]);
          fputs("    endloop\n  endfacet\n", fp);
        }
      }

      // receive data from other processes
      MPI_Status status;
      for (int j = 1; j < comm->nprocs; ++j) {
        MPI_Recv(stldata, 12 * maxobjs, MPI_FLOAT, MPI_ANY_SOURCE, 0, world, &status);
        int numtriangles = 0;
        MPI_Get_count(&status, MPI_FLOAT, &numtriangles);
        numtriangles /= 12;

        // write out received data
        if (binary) {
          uint16_t attributes = 0;
          for (int i = 0; i < numtriangles; ++i) {
            fwrite(stldata + 12 * i, sizeof(float), 12, fp);
            fwrite(&attributes, sizeof(uint16_t), 1, fp);
          }
        } else {
          for (int i = 0; i < numtriangles; ++i) {
            utils::print(fp, "  facet normal {:e} {:e} {:e}\n", stldata[12 * i],
                         stldata[12 * i + 1], stldata[12 * i + 2]);
            fputs("    outer loop\n", fp);
            utils::print(fp, "    vertex {:e} {:e} {:e}\n", stldata[12 * i + 3],
                         stldata[12 * i + 4], stldata[12 * i + 5]);
            utils::print(fp, "    vertex {:e} {:e} {:e}\n", stldata[12 * i + 6],
                         stldata[12 * i + 7], stldata[12 * i + 8]);
            utils::print(fp, "    vertex {:e} {:e} {:e}\n", stldata[12 * i + 9],
                         stldata[12 * i + 10], stldata[12 * i + 11]);
            fputs("    endloop\n  endfacet\n", fp);
          }
        }
      }
      if (!binary) fprintf(fp, "endsolid %s\n", title.c_str());
      delete[] filecurrent;
    } else {
      MPI_Send(stldata, 12 * numobjs, MPI_FLOAT, 0, 0, world);
    }
    delete[] stldata;
  }
}

/* ----------------------------------------------------------------------
   provide graphics information to dump image
------------------------------------------------------------------------- */

int FixGraphicsIsosurface::image(int *&objs, double **&parms)
{
  objs = imgobjs;
  parms = imgparms;
  return numobjs;
}
