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

/* ----------------------------------------------------------------------
   Contributing authors: Mike Brown (ORNL), brownw@ornl.gov
------------------------------------------------------------------------- */

#ifndef PAIR_PRECISION_H
#define PAIR_PRECISION_H

struct _lgpu_float2 {
  float x; float y;
};

struct _lgpu_float4 {
  float x; float y; float z; float w;
};

struct _lgpu_double2 {
  double x; double y;
};

struct _lgpu_double4 {
  double x; double y; double z; double w;
};

#include <iostream>
inline std::ostream & operator<<(std::ostream &out, const _lgpu_float2 &v) {
  out << v.x << " " << v.y;
  return out;
}
  
inline std::ostream & operator<<(std::ostream &out, const _lgpu_float4 &v) {
  out << v.x << " " << v.y << " " << v.z;
  return out;
}
  
inline std::ostream & operator<<(std::ostream &out, const _lgpu_double2 &v) {
  out << v.x << " " << v.y;
  return out;
}
  
inline std::ostream & operator<<(std::ostream &out, const _lgpu_double4 &v) {
  out << v.x << " " << v.y << " " << v.z;
  return out;
}

// PRECISION - Precision for rsq, energy, force, and torque calculation
// ACC_PRECISION - Precision for accumulation of energies, forces, and torques
#ifdef _SINGLE_DOUBLE
#define OCL_PRECISION_COMPILE "-D_SINGLE_DOUBLE"
#define PRECISION float
#define ACC_PRECISION double
#define numtyp2 _lgpu_float2
#define numtyp4 _lgpu_float4
#define acctyp4 _lgpu_double4
#endif

#ifdef _DOUBLE_DOUBLE
#define OCL_PRECISION_COMPILE "-D_DOUBLE_DOUBLE"
#define PRECISION double
#define ACC_PRECISION double
#define numtyp2 _lgpu_double2
#define numtyp4 _lgpu_double4
#define acctyp4 _lgpu_double4
#endif

#ifndef PRECISION
#define OCL_PRECISION_COMPILE "-D_SINGLE_SINGLE"
#define PRECISION float
#define ACC_PRECISION float
#define numtyp2 _lgpu_float2
#define numtyp4 _lgpu_float4
#define acctyp4 _lgpu_float4
#endif

#define MAX_SHARED_TYPES 8
enum{SPHERE_SPHERE,SPHERE_ELLIPSE,ELLIPSE_SPHERE,ELLIPSE_ELLIPSE};

#endif

