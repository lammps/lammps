/***************************************************************************
                              pair_gpu_texture.h
                             -------------------
                               W. Michael Brown

  Tricks for templating textures

 __________________________________________________________________________
    This file is part of the LAMMPS GPU Library
 __________________________________________________________________________

    begin                : Tue Jun 23 2009
    copyright            : (C) 2009 by W. Michael Brown
    email                : wmbrown@sandia.gov
 ***************************************************************************/

/* -----------------------------------------------------------------------
   Copyright (2009) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.
   ----------------------------------------------------------------------- */

#include "nvc_memory.h"

#ifndef PAIR_GPU_TEXTURE_H
#define PAIR_GPU_TEXTURE_H

#ifdef _SINGLE_DOUBLE
#define GB_GPU_DOUBLE
#endif

#ifdef _DOUBLE_DOUBLE
#define GB_GPU_DOUBLE
#endif

template <class numtyp> class cu_vec_traits;
template <> class cu_vec_traits<float> { public: typedef float2 vec2; };
template <> class cu_vec_traits<double> { public: typedef double2 vec2; };

// ------------------------------- x ------------------------------------

static texture<float, 2, cudaReadModeElementType> x_float_tex;
static texture<int2, 2, cudaReadModeElementType> x_double_tex;
template <class numtyp> inline void x_bind_texture(NVC_ConstMatT &mat)
  { mat.bind_texture(x_float_tex); }

template <> inline void x_bind_texture<double>(NVC_ConstMatD &mat)
  { mat.bind_texture(x_double_tex); } 
template <class numtyp> inline void x_unbind_texture()
  { cudaUnbindTexture(x_float_tex); }
template <> inline void x_unbind_texture<double>() 
  { cudaUnbindTexture(x_double_tex); }
template <class numtyp>
static __inline__ __device__ numtyp _x_(const int i, const int j) {
  return tex2D(x_float_tex,i,j);
}
#ifdef GB_GPU_DOUBLE
template <>
static __inline__ __device__ double _x_<double>(const int i,const int j) {
  int2 t=tex2D(x_double_tex,i,j);
  return __hiloint2double(t.y, t.x);
}
#endif

// ------------------------------- form ------------------------------------

static texture<int, 2, cudaReadModeElementType> form_tex;
inline void form_bind_texture(NVC_ConstMatI &mat)
  { mat.bind_texture(form_tex); }
inline void form_unbind_texture()
  { cudaUnbindTexture(form_tex); }
static __inline__ __device__ int _form_(const int i, const int j) {
  return tex2D(form_tex,i,j);
}

// ------------------------------- lshape ------------------------------------

static texture<float, 1, cudaReadModeElementType> lshape_float_tex;
static texture<int2, 1, cudaReadModeElementType> lshape_double_tex;
static cudaChannelFormatDesc channel_lshape;
template <class numtyp> inline void lshape_bind_texture(NVC_VecT &vec) 
 { vec.bind_texture(lshape_float_tex,channel_lshape); }
template <> inline void lshape_bind_texture<double>(NVC_VecD &vec)
 { vec.bind_texture(lshape_double_tex,channel_lshape); }
template <class numtyp> inline void lshape_unbind_texture() 
  { cudaUnbindTexture(lshape_float_tex); }
template <> inline void lshape_unbind_texture<double>() 
  { cudaUnbindTexture(lshape_double_tex); }
template <class numtyp>
static __inline__ __device__ numtyp _lshape_(const int i) 
  { return tex1Dfetch(lshape_float_tex,i); }
#ifdef GB_GPU_DOUBLE
template <>
static __inline__ __device__ double _lshape_<double>(const int i) {
  int2 t=tex1Dfetch(lshape_double_tex,i);
  return __hiloint2double(t.y, t.x);
}
#endif

// ------------------------------- shape ------------------------------------

static texture<float, 2, cudaReadModeElementType> shape_float_tex;
static texture<int2, 2, cudaReadModeElementType> shape_double_tex;
template <class numtyp> inline void shape_bind_texture(NVC_ConstMatT &mat) 
 { mat.bind_texture(shape_float_tex); }
template <> inline void shape_bind_texture<double>(NVC_ConstMatD &mat) 
 { mat.bind_texture(shape_double_tex); } 
template <class numtyp> inline void shape_unbind_texture() 
  { cudaUnbindTexture(shape_float_tex); }
template <> inline void shape_unbind_texture<double>() 
  { cudaUnbindTexture(shape_double_tex); }
template <class numtyp> 
static __inline__ __device__ numtyp _shape_(const int i, const int j) {
  return tex2D(shape_float_tex,j,i);
}
#ifdef GB_GPU_DOUBLE
template <>
static __inline__ __device__ double _shape_<double>(const int i, const int j) {
  int2 t=tex2D(shape_double_tex,j,i);
  return __hiloint2double(t.y, t.x);
}
#endif

// ------------------------------- well ------------------------------------

static texture<float, 2, cudaReadModeElementType> well_float_tex;
static texture<int2, 2, cudaReadModeElementType> well_double_tex;
template <class numtyp> inline void well_bind_texture(NVC_ConstMatT &mat) 
  { mat.bind_texture(well_float_tex); }
template <> inline void well_bind_texture<double>(NVC_ConstMatD &mat)
  { mat.bind_texture(well_double_tex); }
template <class numtyp> inline void well_unbind_texture() 
  { cudaUnbindTexture(well_float_tex); }
template <> inline void well_unbind_texture<double>() 
  { cudaUnbindTexture(well_double_tex); }
template <class numtyp> 
static __inline__ __device__ numtyp _well_(const int i, const int j) 
  { return tex2D(well_float_tex,j,i); }
#ifdef GB_GPU_DOUBLE
template <>
static __inline__ __device__ double _well_<double>(const int i,const int j) {
  int2 t=tex2D(well_double_tex,j,i);
  return __hiloint2double(t.y, t.x);
}
#endif

// ------------------------------- sigma ------------------------------------

static texture<float, 2, cudaReadModeElementType> sigma_float_tex;
static texture<int2, 2, cudaReadModeElementType> sigma_double_tex;
template <class numtyp> inline void sigma_bind_texture(NVC_ConstMatT &mat)
  { mat.bind_texture(sigma_float_tex); }
template <> inline void sigma_bind_texture<double>(NVC_ConstMatD &mat)
  { mat.bind_texture(sigma_double_tex); }
template <class numtyp> inline void sigma_unbind_texture() 
  { cudaUnbindTexture(sigma_float_tex); }
template <> inline void sigma_unbind_texture<double>() 
  { cudaUnbindTexture(sigma_double_tex); }
template <class numtyp>
static __inline__ __device__ numtyp _sigma_(const int i, const int j) {
  return tex2D(sigma_float_tex,j,i);
}
#ifdef GB_GPU_DOUBLE
template <>
static __inline__ __device__ double _sigma_<double>(const int i,const int j) {
  int2 t=tex2D(sigma_double_tex,j,i);
  return __hiloint2double(t.y, t.x);
}
#endif

// ------------------------------- epsilon ------------------------------------

static texture<float, 2, cudaReadModeElementType> epsilon_float_tex;
static texture<int2, 2, cudaReadModeElementType> epsilon_double_tex;
template <class numtyp> inline void epsilon_bind_texture(NVC_ConstMatT &mat)
  { mat.bind_texture(epsilon_float_tex); }
template <> inline void epsilon_bind_texture<double>(NVC_ConstMatD &mat)
  { mat.bind_texture(epsilon_double_tex); }
template <class numtyp> inline void epsilon_unbind_texture() 
  { cudaUnbindTexture(epsilon_float_tex); }
template <> inline void epsilon_unbind_texture<double>() 
  { cudaUnbindTexture(epsilon_double_tex); }
template <class numtyp>
static __inline__ __device__ numtyp _epsilon_(const int i, const int j) {
  return tex2D(epsilon_float_tex,j,i);
}
#ifdef GB_GPU_DOUBLE
template <>
static __inline__ __device__ double _epsilon_<double>(const int i,const int j) {
  int2 t=tex2D(epsilon_double_tex,j,i);
  return __hiloint2double(t.y, t.x);
}
#endif

// ------------------------------- cutsq ------------------------------------

static texture<float, 2, cudaReadModeElementType> cutsq_float_tex;
static texture<int2, 2, cudaReadModeElementType> cutsq_double_tex;
template <class numtyp> inline void cutsq_bind_texture(NVC_ConstMatT &mat)
  { mat.bind_texture(cutsq_float_tex); }
template <> inline void cutsq_bind_texture<double>(NVC_ConstMatD &mat)
  { mat.bind_texture(cutsq_double_tex); }
template <class numtyp> inline void cutsq_unbind_texture() 
 { cudaUnbindTexture(cutsq_float_tex); }
template <> inline void cutsq_unbind_texture<double>() 
 { cudaUnbindTexture(cutsq_double_tex); }
template <class numtyp>
static __inline__ __device__ numtyp _cutsq_(const int i, const int j) {
  return tex2D(cutsq_float_tex,j,i);
}
#ifdef GB_GPU_DOUBLE
template <>
static __inline__ __device__ double _cutsq_<double>(const int i,const int j) {
  int2 t=tex2D(cutsq_double_tex,j,i);
  return __hiloint2double(t.y, t.x);
}
#endif

// ------------------------------- lj1 ------------------------------------

static texture<float2, 2, cudaReadModeElementType> lj1_float_tex;
static texture<int4, 2, cudaReadModeElementType> lj1_double_tex;
template <class numtyp> inline void lj1_bind_texture(NVC_ConstMatT &mat)
  { mat.bind_texture(lj1_float_tex); }
template <> inline void lj1_bind_texture<double2>(NVC_ConstMatD2 &mat)
  { mat.bind_texture(lj1_double_tex); }
template <class numtyp> inline void lj1_unbind_texture() 
 { cudaUnbindTexture(lj1_float_tex); }
template <> inline void lj1_unbind_texture<double2>() 
 { cudaUnbindTexture(lj1_double_tex); }
template <class numtyp>
static __inline__ __device__ 
typename cu_vec_traits<numtyp>::vec2 _lj1_(const int i, const int j) {
  return tex2D(lj1_float_tex,j,i);
}
#ifdef GB_GPU_DOUBLE
template <>
static __inline__ __device__ double2 _lj1_<double>(const int i,const int j) {
  int4 t=tex2D(lj1_double_tex,j,i);
  double2 ans;
  ans.x=__hiloint2double(t.y, t.x);
  ans.y=__hiloint2double(t.w, t.z);
  return ans;
}
#endif

// ------------------------------- lj3 ------------------------------------

static texture<float2, 2, cudaReadModeElementType> lj3_float_tex;
static texture<int4, 2, cudaReadModeElementType> lj3_double_tex;
template <class numtyp> inline void lj3_bind_texture(NVC_ConstMatT &mat)
  { mat.bind_texture(lj3_float_tex); }
template <> inline void lj3_bind_texture<double2>(NVC_ConstMatD2 &mat)
  { mat.bind_texture(lj3_double_tex); }
template <class numtyp> inline void lj3_unbind_texture() 
 { cudaUnbindTexture(lj3_float_tex); }
template <> inline void lj3_unbind_texture<double2>() 
 { cudaUnbindTexture(lj3_double_tex); }
template <class numtyp>
static __inline__ __device__ 
typename cu_vec_traits<numtyp>::vec2 _lj3_(const int i, const int j) {
  return tex2D(lj3_float_tex,j,i);
}
#ifdef GB_GPU_DOUBLE
template <>
static __inline__ __device__ double2 _lj3_<double>(const int i,const int j) {
  int4 t=tex2D(lj3_double_tex,j,i);
  double2 ans;
  ans.x=__hiloint2double(t.y, t.x);
  ans.y=__hiloint2double(t.w, t.z);
  return ans;
}
#endif

// ------------------------------- offset ------------------------------------

static texture<float, 2, cudaReadModeElementType> offset_float_tex;
static texture<int2, 2, cudaReadModeElementType> offset_double_tex;
template <class numtyp> inline void offset_bind_texture(NVC_ConstMatT &mat)
  { mat.bind_texture(offset_float_tex); }
template <> inline void offset_bind_texture<double>(NVC_ConstMatD &mat)
  { mat.bind_texture(offset_double_tex); }
template <class numtyp> inline void offset_unbind_texture() 
 { cudaUnbindTexture(offset_float_tex); }
template <> inline void offset_unbind_texture<double>() 
 { cudaUnbindTexture(offset_double_tex); }
template <class numtyp>
static __inline__ __device__ numtyp _offset_(const int i, const int j) {
  return tex2D(offset_float_tex,j,i);
}
#ifdef GB_GPU_DOUBLE
template <>
static __inline__ __device__ double _offset_<double>(const int i,const int j) {
  int2 t=tex2D(offset_double_tex,j,i);
  return __hiloint2double(t.y, t.x);
}
#endif

#endif
