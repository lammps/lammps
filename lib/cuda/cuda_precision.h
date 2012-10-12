/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator

   Original Version:
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   See the README file in the top-level LAMMPS directory.

   -----------------------------------------------------------------------

   USER-CUDA Package and associated modifications:
   https://sourceforge.net/projects/lammpscuda/

   Christian Trott, christian.trott@tu-ilmenau.de
   Lars Winterfeld, lars.winterfeld@tu-ilmenau.de
   Theoretical Physics II, University of Technology Ilmenau, Germany

   See the README file in the USER-CUDA directory.

   This software is distributed under the GNU General Public License.
------------------------------------------------------------------------- */

#ifndef CUDA_PRECISION_H_
#define CUDA_PRECISION_H_
/* This File gives Type definitions for mixed precision calculation in the cuda part of LAMMPS-CUDA.
 * Predefined behaviour is given by global CUDA_PRECISION (can be overwritten during compilation).
 * ***_FLOAT: type definition of given property
 * ***_F: constant extension in code (1.0 is interpreted as double while 1.0f is interpreted as float, now use: 1.0CUDA_F)
 */

#ifdef CUDA_USE_BINNING
#define CUDA_IF_BINNING(a) a
#else
#define CUDA_IF_BINNING(a)
#endif

//GLOBAL

#ifdef CUDA_PRECISION
#if CUDA_PRECISION == 1
#define CUDA_FLOAT float
#define CUDA_F(x) x##f
#endif
#if CUDA_PRECISION == 2
#define CUDA_FLOAT double
#define CUDA_F(x) x
#endif
#endif

#ifndef CUDA_PRECISION
#define CUDA_FLOAT double
#define CUDA_F(x) x
#define CUDA_PRECISION 2
#endif
//--------------------------------
//-----------FFT-----------------
//--------------------------------

#ifdef FFT_PRECISION_CU
#if FFT_PRECISION_CU == 1
#define FFT_FLOAT float
#define FFT_F(x) x##f
#endif
#if FFT_PRECISION_CU == 2
#define FFT_FLOAT double
#define FFT_F(x) x
#endif
#endif

#ifndef FFT_PRECISION_CU
#define FFT_FLOAT CUDA_FLOAT
#define FFT_F(x) CUDA_F(x)
#define FFT_PRECISION_CU CUDA_PRECISION
#endif

//--------------------------------
//-----------PPPM-----------------
//--------------------------------

#ifndef PPPM_PRECISION
#define PPPM_PRECISION CUDA_PRECISION
#endif

#ifdef PPPM_PRECISION
#if PPPM_PRECISION == 1
#define PPPM_FLOAT float
#ifdef float3
#define PPPM_FLOAT3 float3
#else
struct PPPM_FLOAT3 {
  PPPM_FLOAT x;
  PPPM_FLOAT y;
  PPPM_FLOAT z;
};
#endif
#define PPPM_F(x) x##f
#endif
#if PPPM_PRECISION == 2
#define PPPM_FLOAT double
struct PPPM_FLOAT3 {
  PPPM_FLOAT x;
  PPPM_FLOAT y;
  PPPM_FLOAT z;
};
#define PPPM_F(x) x
#endif
#endif


//--------------------------------
//-----------FORCE-----------------
//--------------------------------


#ifdef F_PRECISION
#if F_PRECISION == 1
#define F_FLOAT float
#define F_F(x) x##f
#endif
#if F_PRECISION == 2
#define F_FLOAT double
#define F_F(x) x
#endif
#endif

#ifndef F_PRECISION
#define F_FLOAT CUDA_FLOAT
#define F_F(x) CUDA_F(x)
#define F_PRECISION CUDA_PRECISION
#endif

#if F_PRECISION == 1
#define _SQRT_ sqrtf
#define _RSQRT_ rsqrtf
#define _EXP_ expf
#else
#define _SQRT_ sqrt
#define _RSQRT_ rsqrt
#define _EXP_ exp
#endif

#if F_PRECISION == 2
struct F_FLOAT2 {
  F_FLOAT x;
  F_FLOAT y;
};
struct F_FLOAT3 {
  F_FLOAT x;
  F_FLOAT y;
  F_FLOAT z;
};
struct F_FLOAT4 {
  F_FLOAT x;
  F_FLOAT y;
  F_FLOAT z;
  F_FLOAT w;
};
#else
#define F_FLOAT2 float2
#define F_FLOAT3 float3
#define F_FLOAT4 float4
#endif
//--------------------------------
//-----------ENERGY-----------------
//--------------------------------

#ifndef ENERGY_PRECISION
#define ENERGY_FLOAT CUDA_FLOAT
#define ENERGY_F(x) CUDA_F(x)
#endif

#ifdef ENERGY_PRECISION
#if ENERGY_PRECISION == 1
#define ENERGY_FLOAT float
#define ENERGY_F(x) x##f
#endif
#if ENERGY_PRECISION == 2
#define ENERGY_FLOAT double
#define ENERGY_F(x) x
#endif
#endif

#ifndef ENERGY_PRECISION
#define ENERGY_FLOAT CUDA_FLOAT
#define ENERGY_F(x) CUDA_F(x)
#define ENERGY_PRECISION CUDA_PRECISION
#endif

//--------------------------------
//-----------POSITIONS------------
//--------------------------------

#ifdef X_PRECISION
#if X_PRECISION == 1
#define X_FLOAT float
#define X_F(x) x##f
#endif
#if X_PRECISION == 2
#define X_FLOAT double
#define X_F(x) x
#endif
#endif

#ifndef X_PRECISION
#define X_FLOAT CUDA_FLOAT
#define X_F(x) CUDA_F(x)
#define X_PRECISION CUDA_PRECISION
#endif

#if X_PRECISION == 2
struct X_FLOAT2 {
  X_FLOAT x;
  X_FLOAT y;
};
struct X_FLOAT3 {
  X_FLOAT x;
  X_FLOAT y;
  X_FLOAT z;
};
struct X_FLOAT4 {
  X_FLOAT x;
  X_FLOAT y;
  X_FLOAT z;
  X_FLOAT w;
};
#else
#define X_FLOAT2 float2
#define X_FLOAT3 float3
#define X_FLOAT4 float4
#endif

//--------------------------------
//-----------velocities-----------
//--------------------------------

#ifdef V_PRECISION
#if V_PRECISION == 1
#define V_FLOAT float
#define V_F(x) x##f
#endif
#if V_PRECISION == 2
#define V_FLOAT double
#define V_F(x) x
#endif
#endif

#ifndef V_PRECISION
#define V_FLOAT CUDA_FLOAT
#define V_F(x) CUDA_F(x)
#define V_PRECISION CUDA_PRECISION
#endif

#if V_PRECISION == 2
struct V_FLOAT4 {
  V_FLOAT x;
  V_FLOAT y;
  V_FLOAT z;
  V_FLOAT w;
};
#else
#define V_FLOAT4 float4
#endif

#ifdef NO_PREC_TIMING
struct timespec_2 {
  unsigned int tv_sec;
  unsigned int tv_nsec;
};

#define timespec timespec_2
#define clock_gettime(a,b)
#endif
#endif /*CUDA_PRECISION_H_*/
