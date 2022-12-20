// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   ----------------------------------------------------------------------- */

#include "manifold_thylakoid_shared.h"
#include <cmath>

using namespace LAMMPS_NS;
using namespace user_manifold;


thyla_part::thyla_part( int type, double *args, double xlo, double ylo, double zlo,
                        double xhi, double yhi, double zhi )
  : type(type), xlo(xlo), xhi(xhi),
    ylo(ylo), yhi(yhi), zlo(zlo), zhi(zhi)
{
  switch(type) {
    case THYLA_TYPE_PLANE: // a*(x-x0) + b*(y-y0) + c*(z-z0) = 0
      params[0] = args[0]; // a
      params[1] = args[1]; // b
      params[2] = args[2]; // c
      params[3] = args[3]; // x0
      params[4] = args[4]; // y0
      params[5] = args[5]; // z0
      break;
    case THYLA_TYPE_SPHERE: // (x-x0)^2 + (y-y0)^2 + (z-z0)^2 - R^2 = 0
      params[0] = args[0]; // R
      params[1] = args[1]; // x0
      params[2] = args[2]; // y0
      params[3] = args[3]; // z0
      break;
    case THYLA_TYPE_CYL: // a*(x-x0)^2 + b*(y-y0)^2 + c*(z-z0)^2 - R^2 = 0
      params[0] = args[0]; // a
      params[1] = args[1]; // b
      params[2] = args[2]; // c
      params[3] = args[3]; // x0
      params[4] = args[4]; // y0
      params[5] = args[5]; // z0
      params[6] = args[6]; // R
      if ((args[0] != 0.0) && (args[1] != 0.0) && (args[2] != 0.0)) {
        err_flag = -1;
        return;
      }
      // The others should be 1.
      if ( (args[0] != 1.0) && (args[0] != 0.0) &&
           (args[1] != 1.0) && (args[1] != 0.0) &&
           (args[2] != 1.0) && (args[2] != 0.0)) {
        err_flag = -1;
      }
      break;
    case THYLA_TYPE_CYL_TO_PLANE:
      /*
       * Funky bit that connects a cylinder to a plane.
       * It is what you get by rotating the equation
       * r(x) = R0 + R*( 1 - sqrt( 1 - ( ( X0 - x ) /R )^2 ) ) around the x-axis.
       * I kid you not.
       *
       * The shape is symmetric so you have to set whether it is the "left" or
       * "right" end by truncating it properly. It is designed to run from
       * X0 to X0 + R if "right" or X0 - R to X0 if "left".
       *
       * As params it takes X0, R0, R, and a sign that determines whether it is
       * "left" (args[3] < 0) or "right" (args[3] > 0).
       *
       * The args are: X0, R0, R, x0, y0, z0, sign.
       */
      params[0] = args[0];
      params[1] = args[1];
      params[2] = args[2];
      params[3] = args[3];
      params[4] = args[4];
      params[5] = args[5];
      params[6] = args[6];
      break;
    default:
      err_flag = -1;
  }
  x0 = (type == THYLA_TYPE_SPHERE) ? params[1] : params[3];
  y0 = (type == THYLA_TYPE_SPHERE) ? params[2] : params[4];
  z0 = (type == THYLA_TYPE_SPHERE) ? params[3] : params[5];
}

double thyla_part::g(const double *x)
{
  switch(type) {
    case THYLA_TYPE_PLANE:{ // a*(x-x0) + b*(y-y0) + c*(z-z0) = 0
      double a  = params[0];
      double b  = params[1];
      double c  = params[2];
      double dx = x[0] - params[3];
      double dy = x[1] - params[4];
      double dz = x[2] - params[5];
      return a*dx + b*dy + c*dz;
      break;
    }
    case THYLA_TYPE_SPHERE:{ // (x-x0)^2 + (y-y0)^2 + (z-z0)^2 - R^2 = 0
      double R2 = params[0]*params[0];
      double dx = x[0] - params[1];
      double dy = x[1] - params[2];
      double dz = x[2] - params[3];
      return dx*dx + dy*dy + dz*dz - R2;

      break;
    }
    case THYLA_TYPE_CYL:{ // a*(x-x0)^2 + b*(y-y0)^2 + c*(z-z0)^2 - R^2 = 0
      double a  = params[0];
      double b  = params[1];
      double c  = params[2];
      double X0 = params[3];
      double Y0 = params[4];
      double Z0 = params[5];
      double R  = params[6];
      double dx = x[0] - X0;
      double dy = x[1] - Y0;
      double dz = x[2] - Z0;

      return a*dx*dx + b*dy*dy + c*dz*dz - R*R;
      break;
    }
    case THYLA_TYPE_CYL_TO_PLANE:{
      double X0 = params[0];
      double R0 = params[1];
      double R  = params[2];

      // Determine the centre of the sphere.
      double dx   = (x[0] - X0);
      double dyz  = sqrt( x[1]*x[1] + x[2]*x[2] );
      double rdyz = dyz - (R0 + R);
      double norm = 1.0 / (2.0*R);
      // Maybe sign is important here...
      double g = norm*(dx*dx + rdyz*rdyz - R*R);
      return g;

    }
    default:
      err_flag = -1;
      return 0; // Mostly to get rid of compiler werrors.
      break;
  }
}



void   thyla_part::n( const double *x, double *n )
{
  switch(type) {
    case THYLA_TYPE_PLANE:{ // a*(x-x0) + b*(y-y0) + c*(z-z0) = 0
      double a  = params[0];
      double b  = params[1];
      double c  = params[2];
      n[0] = a;
      n[1] = b;
      n[2] = c;
      break;
    }
    case THYLA_TYPE_SPHERE:{ // (x-x0)^2 + (y-y0)^2 + (z-z0)^2 - R^2 = 0
      double dx = x[0] - params[1];
      double dy = x[1] - params[2];
      double dz = x[2] - params[3];
      n[0] = 2*dx;
      n[1] = 2*dy;
      n[2] = 2*dz;
      break;
    }
    case THYLA_TYPE_CYL:{ // a*(x-x0)^2 + b*(y-y0)^2 + c*(z-z0)^2 - R^2 = 0
      double a  = params[0];
      double b  = params[1];
      double c  = params[2];
      double X0 = params[3];
      double Y0 = params[4];
      double Z0 = params[5];
      double dx = x[0] - X0;
      double dy = x[1] - Y0;
      double dz = x[2] - Z0;

      n[0] = 2*a*dx;
      n[1] = 2*b*dy;
      n[2] = 2*c*dz;
      break;
    }
    case THYLA_TYPE_CYL_TO_PLANE:{
      double X0 = params[0];
      double R0 = params[1];
      double R  = params[2];
      double s  = (params[6] > 0.0) ? 1.0 : -1.0;

      // Determine the centre of the sphere.
      double dx   = s*(x[0] - X0);
      double ryz  = sqrt( x[1]*x[1] + x[2]*x[2] );
      // Maybe sign is important here...
      // Normalize g and n so that the normal is continuous:
      double norm = 1.0 / (2.0 * R);

      n[0] = s*2*dx*norm;

      double const_part = 1.0 - (R0 + R) / ryz;
      n[1] = 2*x[1]*const_part*norm;
      n[2] = 2*x[2]*const_part*norm;
      break;
    }
    default:
      err_flag = -1;
      break;
  }
}


void thyla_part_geom::mirror( unsigned int axis, thyla_part_geom *m,
                              const thyla_part_geom *o )
{
  // Since dir is already the index of the array this is really simple:
  m->lo[0] = o->lo[0];
  m->lo[1] = o->lo[1];
  m->lo[2] = o->lo[2];
  m->pt[0] = o->pt[0];
  m->pt[1] = o->pt[1];
  m->pt[2] = o->pt[2];
  m->hi[0] = o->hi[0];
  m->hi[1] = o->hi[1];
  m->hi[2] = o->hi[2];

  m->lo[axis] = -o->hi[axis];
  m->hi[axis] = -o->lo[axis];
  m->pt[axis] = -o->pt[axis];
}


