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

#include "nvc_traits.h"
#include "nvc_memory.h"

#ifndef PAIR_GPU_TEXTURE_H
#define PAIR_GPU_TEXTURE_H

#ifdef _SINGLE_DOUBLE
#define GB_GPU_DOUBLE
#endif

#ifdef _DOUBLE_DOUBLE
#define GB_GPU_DOUBLE
#endif

// ------------------------------- x ------------------------------------

static texture<float, 2, cudaReadModeElementType> x_float_tex;
static texture<int2, 2, cudaReadModeElementType> x_double_tex;
template <class numtyp> inline textureReference * x_get_texture() {
  const textureReference *ptr;
  cudaGetTextureReference(&ptr,"x_float_tex");
  return const_cast<textureReference *>(ptr);
}
template <> inline textureReference * x_get_texture<double>() {
  const textureReference *ptr;
  cudaGetTextureReference(&ptr,"x_double_tex");
  return const_cast<textureReference *>(ptr);
}
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
inline textureReference * form_get_texture() {
  const textureReference *ptr;
  cudaGetTextureReference(&ptr,"form_tex");
  return const_cast<textureReference *>(ptr);
}
static __inline__ __device__ int _form_(const int i, const int j) {
  return tex2D(form_tex,i,j);
}

// ------------------------------- lshape ------------------------------------

static texture<float, 1, cudaReadModeElementType> lshape_float_tex;
static texture<int2, 1, cudaReadModeElementType> lshape_double_tex;
template <class numtyp> inline textureReference * lshape_get_texture() {
  const textureReference *ptr;
  cudaGetTextureReference(&ptr,"lshape_float_tex");
  return const_cast<textureReference *>(ptr);
}
template <> inline textureReference * lshape_get_texture<double>() {
  const textureReference *ptr;
  cudaGetTextureReference(&ptr,"lshape_double_tex");
  return const_cast<textureReference *>(ptr);
}
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
template <class numtyp> inline textureReference * shape_get_texture() {
  const textureReference *ptr;
  cudaGetTextureReference(&ptr,"shape_float_tex");
  return const_cast<textureReference *>(ptr);
}
template <> inline textureReference * shape_get_texture<double>() {
  const textureReference *ptr;
  cudaGetTextureReference(&ptr,"shape_double_tex");
  return const_cast<textureReference *>(ptr);
}
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
template <class numtyp> inline textureReference * well_get_texture() {
  const textureReference *ptr;
  cudaGetTextureReference(&ptr,"well_float_tex");
  return const_cast<textureReference *>(ptr);
}
template <> inline textureReference * well_get_texture<double>() {
  const textureReference *ptr;
  cudaGetTextureReference(&ptr,"well_double_tex");
  return const_cast<textureReference *>(ptr);
}
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
template <class numtyp> inline textureReference * sigma_get_texture() {
  const textureReference *ptr;
  cudaGetTextureReference(&ptr,"sigma_float_tex");
  return const_cast<textureReference *>(ptr);
}
template <> inline textureReference * sigma_get_texture<double>() {
  const textureReference *ptr;
  cudaGetTextureReference(&ptr,"sigma_double_tex");
  return const_cast<textureReference *>(ptr);
}
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
template <class numtyp> inline textureReference * epsilon_get_texture() {
  const textureReference *ptr;
  cudaGetTextureReference(&ptr,"epsilon_float_tex");
  return const_cast<textureReference *>(ptr);
}
template <> inline textureReference * epsilon_get_texture<double>() {
  const textureReference *ptr;
  cudaGetTextureReference(&ptr,"epsilon_double_tex");
  return const_cast<textureReference *>(ptr);
}
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
template <class numtyp> inline textureReference * cutsq_get_texture() {
  const textureReference *ptr;
  cudaGetTextureReference(&ptr,"cutsq_float_tex");
  return const_cast<textureReference *>(ptr);
}
template <> inline textureReference * cutsq_get_texture<double>() {
  const textureReference *ptr;
  cudaGetTextureReference(&ptr,"cutsq_double_tex");
  return const_cast<textureReference *>(ptr);
}
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
template <class numtyp> inline textureReference * lj1_get_texture() {
  const textureReference *ptr;
  cudaGetTextureReference(&ptr,"lj1_float_tex");
  return const_cast<textureReference *>(ptr);
}
template <> inline textureReference * lj1_get_texture<double>() {
  const textureReference *ptr;
  cudaGetTextureReference(&ptr,"lj1_double_tex");
  return const_cast<textureReference *>(ptr);
}
template <class numtyp>
static __inline__ __device__ 
typename nvc_vec_traits<numtyp>::vec2 _lj1_(const int i, const int j) {
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
template <class numtyp> inline textureReference * lj3_get_texture() {
  const textureReference *ptr;
  cudaGetTextureReference(&ptr,"lj3_float_tex");
  return const_cast<textureReference *>(ptr);
}
template <> inline textureReference * lj3_get_texture<double>() {
  const textureReference *ptr;
  cudaGetTextureReference(&ptr,"lj3_double_tex");
  return const_cast<textureReference *>(ptr);
}
template <class numtyp>
static __inline__ __device__ 
typename nvc_vec_traits<numtyp>::vec2 _lj3_(const int i, const int j) {
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
template <class numtyp> inline textureReference * offset_get_texture() {
  const textureReference *ptr;
  cudaGetTextureReference(&ptr,"offset_float_tex");
  return const_cast<textureReference *>(ptr);
}
template <> inline textureReference * offset_get_texture<double>() {
  const textureReference *ptr;
  cudaGetTextureReference(&ptr,"offset_double_tex");
  return const_cast<textureReference *>(ptr);
}
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
