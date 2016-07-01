/***************************************************************************
                                 precision.h
                             -------------------
                            W. Michael Brown (ORNL)

  Data and preprocessor definitions for different precision modes

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : brownw@ornl.gov
 ***************************************************************************/

#ifndef LAL_PRECISION_H
#define LAL_PRECISION_H

#if defined(USE_CUDART)
#include <cuda_runtime.h>
#endif

struct _lgpu_int2 {
  int x; int y;
};

#ifndef int2
#define int2 _lgpu_int2
#endif

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

enum{SPHERE_SPHERE,SPHERE_ELLIPSE,ELLIPSE_SPHERE,ELLIPSE_ELLIPSE};

// OCL_DEFAULT_VENDOR: preprocessor define for hardware
// specific sizes of OpenCL kernel related constants

#ifdef FERMI_OCL
#define OCL_DEFAULT_VENDOR "fermi"
#endif

#ifdef KEPLER_OCL
#define OCL_DEFAULT_VENDOR "kepler"
#endif

#ifdef CYPRESS_OCL
#define OCL_DEFAULT_VENDOR "cypress"
#endif

#ifdef GENERIC_OCL
#define OCL_DEFAULT_VENDOR "generic"
#endif

#ifdef INTEL_OCL
#define OCL_DEFAULT_VENDOR "intel"
#endif

#ifdef PHI_OCL
#define OCL_DEFAULT_VENDOR "phi"
#endif

#ifndef OCL_DEFAULT_VENDOR
#define OCL_DEFAULT_VENDOR "none"
#endif

// default to 32-bit smallint and other ints, 64-bit bigint: same as defined in src/lmptype.h
#if !defined(LAMMPS_SMALLSMALL) && !defined(LAMMPS_BIGBIG) && !defined(LAMMPS_SMALLBIG)
#define LAMMPS_SMALLBIG
#endif

#ifdef LAMMPS_SMALLBIG
typedef int tagint;
#define OCL_INT_TYPE "-DLAMMPS_SMALLBIG"
#endif
#ifdef LAMMPS_BIGBIG
#include "inttypes.h"
typedef int64_t tagint;
#define OCL_INT_TYPE "-DLAMMPS_BIGBIG"
#endif
#ifdef LAMMPS_SMALLSMALL
typedef int tagint;
#define OCL_INT_TYPE "-DLAMMPS_SMALLSMALL"
#endif

#endif // LAL_PRECISION_H
