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

// ---------------------- OPENMP PREPROCESSOR STUFF ------------------
#if defined(_OPENMP)
  #if !defined(LAL_USE_OMP)
  #define LAL_USE_OMP 1
  #endif

  #if !defined(LAL_USE_OMP_SIMD)
    #if (_OPENMP >= 201307)
    #define LAL_USE_OMP_SIMD 1
    #else
    #define LAL_USE_OMP_SIMD 0
    #endif
  #endif
#else
  #if !defined(LAL_USE_OMP)
  #define LAL_USE_OMP 0
  #endif

  #if !defined(LAL_USE_OMP_SIMD)
  #define LAL_USE_OMP_SIMD 0
  #endif
#endif

struct _lgpu_int2 {
  int x; int y;
};

#ifndef USE_HIP
#ifndef int2
#define int2 _lgpu_int2
#endif
#endif

struct _lgpu_float2 {
  float x; float y;
};

struct _lgpu_float3 {
  float x; float y; float z;
};

struct _lgpu_float4 {
  float x; float y; float z; float w;
};

struct _lgpu_double2 {
  double x; double y;
};

struct _lgpu_double3 {
  double x; double y; double z;
};

struct _lgpu_double4 {
  double x; double y; double z; double w;
};

#include <iostream>
inline std::ostream & operator<<(std::ostream &out, const _lgpu_float2 &v) {
  out << v.x << " " << v.y;
  return out;
}

inline std::ostream & operator<<(std::ostream &out, const _lgpu_float3 &v) {
  out << v.x << " " << v.y << " " << v.z;
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

inline std::ostream & operator<<(std::ostream &out, const _lgpu_double3 &v) {
  out << v.x << " " << v.y << " " << v.z;
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
#define numtyp3 _lgpu_float3
#define numtyp4 _lgpu_float4
#define acctyp2 _lgpu_double2
#define acctyp3 _lgpu_double3
#define acctyp4 _lgpu_double4
#endif

#ifdef _DOUBLE_DOUBLE
#define OCL_PRECISION_COMPILE "-D_DOUBLE_DOUBLE"
#define PRECISION double
#define ACC_PRECISION double
#define numtyp2 _lgpu_double2
#define numtyp3 _lgpu_double3
#define numtyp4 _lgpu_double4
#define acctyp2 _lgpu_double2
#define acctyp3 _lgpu_double3
#define acctyp4 _lgpu_double4
#endif

#ifndef PRECISION
#define OCL_PRECISION_COMPILE "-D_SINGLE_SINGLE"
#define PRECISION float
#define ACC_PRECISION float
#define numtyp2 _lgpu_float2
#define numtyp3 _lgpu_float3
#define numtyp4 _lgpu_float4
#define acctyp2 _lgpu_float2
#define acctyp3 _lgpu_float3
#define acctyp4 _lgpu_float4
#endif

enum{SPHERE_SPHERE,SPHERE_ELLIPSE,ELLIPSE_SPHERE,ELLIPSE_ELLIPSE};

// default to 32-bit smallint and other ints, 64-bit bigint:
//   same as defined in src/lmptype.h
#if !defined(LAMMPS_SMALLSMALL) && !defined(LAMMPS_BIGBIG) && \
  !defined(LAMMPS_SMALLBIG)
#define LAMMPS_SMALLBIG
#endif

#ifdef LAMMPS_SMALLBIG
typedef int tagint;
#define OCL_INT_TYPE "-DLAMMPS_SMALLBIG"
#endif
#ifdef LAMMPS_BIGBIG
#include "stdint.h"
typedef int64_t tagint;
#define OCL_INT_TYPE "-DLAMMPS_BIGBIG"
#endif
#ifdef LAMMPS_SMALLSMALL
typedef int tagint;
#define OCL_INT_TYPE "-DLAMMPS_SMALLSMALL"
#endif

#endif // LAL_PRECISION_H
