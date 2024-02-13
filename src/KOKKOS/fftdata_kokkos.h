// clang-format off
/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// data types for 2d/3d FFTs

#ifndef LMP_FFT_DATA_KOKKOS_H
#define LMP_FFT_DATA_KOKKOS_H

#include "kokkos_type.h"
#include "lmpfftsettings.h"

// -------------------------------------------------------------------------

// if a user sets FFTW, it means FFTW3

#ifdef LMP_KOKKOS
# ifdef FFT_KOKKOS_FFTW
#  undef FFT_KOKKOS_FFTW
#  define FFT_KOKKOS_FFTW3
# endif
# ifdef FFT_KOKKOS_FFTW_THREADS
#  if !defined(FFT_KOKKOS_FFTW3)
#   error "Must use -DFFT_KOKKOS_FFTW3 with -DFFT_KOKKOS_FFTW_THREADS"
#  endif
# endif
#endif

// with KOKKOS in CUDA or HIP mode we can only have
//  CUFFT/HIPFFT or KISS, thus undefine all other
//  FFTs here

#ifdef KOKKOS_ENABLE_CUDA
# if defined(FFT_KOKKOS_FFTW)
#  undef FFT_KOKKOS_FFTW
# endif
# if defined(FFT_KOKKOS_FFTW3)
#  undef FFT_KOKKOS_FFTW3
# endif
# if defined(FFT_KOKKOS_MKL)
#  undef FFT_KOKKOS_MKL
# endif
# if !defined(FFT_KOKKOS_CUFFT) && !defined(FFT_KOKKOS_KISS)
#  define FFT_KOKKOS_KISS
# endif
#elif defined(KOKKOS_ENABLE_HIP)
# if defined(FFT_KOKKOS_FFTW)
#  undef FFT_KOKKOS_FFTW
# endif
# if defined(FFT_KOKKOS_FFTW3)
#  undef FFT_KOKKOS_FFTW3
# endif
# if defined(FFT_KOKKOS_MKL)
#  undef FFT_KOKKOS_MKL
# endif
# if !defined(FFT_KOKKOS_HIPFFT) && !defined(FFT_KOKKOS_KISS)
#  define FFT_KOKKOS_KISS
# endif
#else
# if defined(FFT_KOKKOS_CUFFT)
#  error "Must enable CUDA with KOKKOS to use -DFFT_KOKKOS_CUFFT"
# endif
# if defined(FFT_KOKKOS_HIPFFT)
#  error "Must enable HIP with KOKKOS to use -DFFT_KOKKOS_HIPFFT"
# endif
#endif

// set strings for library info output

#if defined(FFT_KOKKOS_CUFFT)
#define LMP_FFT_KOKKOS_LIB "cuFFT"
#elif defined(FFT_KOKKOS_HIPFFT)
#define LMP_FFT_KOKKOS_LIB "hipFFT"
#elif defined(FFT_KOKKOS_FFTW3)
#define LMP_FFT_KOKKOS_LIB "FFTW3"
#elif defined(FFT_KOKKOS_MKL)
#define LMP_FFT_KOKKOS_LIB "MKL FFT"
#else
#define LMP_FFT_KOKKOS_LIB "KISS FFT"
#endif


#if defined(FFT_KOKKOS_MKL)
  #include "mkl_dfti.h"
  #if defined(FFT_SINGLE)
    typedef float _Complex FFT_KOKKOS_DATA;
    #define FFT_KOKKOS_MKL_PREC DFTI_SINGLE
  #else
    typedef double _Complex FFT_KOKKOS_DATA;
    #define FFT_KOKKOS_MKL_PREC DFTI_DOUBLE
  #endif
#elif defined(FFT_KOKKOS_FFTW3)
  #include "fftw3.h"
  #if defined(FFT_SINGLE)
    typedef fftwf_complex FFT_KOKKOS_DATA;
    #define FFTW_API(function)  fftwf_ ## function
  #else
    typedef fftw_complex FFT_KOKKOS_DATA;
    #define FFTW_API(function) fftw_ ## function
  #endif
#elif defined(FFT_KOKKOS_CUFFT)
  #include "cufft.h"
  #if defined(FFT_SINGLE)
    #define cufftExec cufftExecC2C
    #define CUFFT_TYPE CUFFT_C2C
    typedef cufftComplex FFT_KOKKOS_DATA;
  #else
    #define cufftExec cufftExecZ2Z
    #define CUFFT_TYPE CUFFT_Z2Z
    typedef cufftDoubleComplex FFT_KOKKOS_DATA;
  #endif
#elif defined(FFT_KOKKOS_HIPFFT)
  #include <hipfft/hipfft.h>
  #if defined(FFT_SINGLE)
    #define hipfftExec hipfftExecC2C
    #define HIPFFT_TYPE HIPFFT_C2C
    typedef hipfftComplex FFT_KOKKOS_DATA;
  #else
    #define hipfftExec hipfftExecZ2Z
    #define HIPFFT_TYPE HIPFFT_Z2Z
    typedef hipfftDoubleComplex FFT_KOKKOS_DATA;
  #endif
#else
  #if defined(FFT_SINGLE)
    #define kiss_fft_scalar float
  #else
    #define kiss_fft_scalar double
  #endif
  typedef struct {
    kiss_fft_scalar re;
    kiss_fft_scalar im;
  } FFT_KOKKOS_DATA;
  #ifndef FFT_KOKKOS_KISS
  #define FFT_KOKKOS_KISS
  #endif
#endif

// (double[2]*) is not a 1D pointer
#if defined(FFT_KOKKOS_FFTW3)
  typedef FFT_SCALAR* FFT_KOKKOS_DATA_POINTER;
#else
  typedef FFT_KOKKOS_DATA* FFT_KOKKOS_DATA_POINTER;
#endif


template <class DeviceType>
struct FFTArrayTypes;

template <>
struct FFTArrayTypes<LMPDeviceType> {

typedef Kokkos::
  DualView<FFT_SCALAR*, Kokkos::LayoutRight, LMPDeviceType> tdual_FFT_SCALAR_1d;
typedef tdual_FFT_SCALAR_1d::t_dev t_FFT_SCALAR_1d;
typedef tdual_FFT_SCALAR_1d::t_dev_um t_FFT_SCALAR_1d_um;

typedef Kokkos::DualView<FFT_SCALAR**,Kokkos::LayoutRight,LMPDeviceType> tdual_FFT_SCALAR_2d;
typedef tdual_FFT_SCALAR_2d::t_dev t_FFT_SCALAR_2d;

typedef Kokkos::DualView<FFT_SCALAR**[3],Kokkos::LayoutRight,LMPDeviceType> tdual_FFT_SCALAR_2d_3;
typedef tdual_FFT_SCALAR_2d_3::t_dev t_FFT_SCALAR_2d_3;

typedef Kokkos::DualView<FFT_SCALAR***,Kokkos::LayoutRight,LMPDeviceType> tdual_FFT_SCALAR_3d;
typedef tdual_FFT_SCALAR_3d::t_dev t_FFT_SCALAR_3d;

typedef Kokkos::
  DualView<FFT_KOKKOS_DATA*, Kokkos::LayoutRight, LMPDeviceType> tdual_FFT_DATA_1d;
typedef tdual_FFT_DATA_1d::t_dev t_FFT_DATA_1d;
typedef tdual_FFT_DATA_1d::t_dev_um t_FFT_DATA_1d_um;

typedef Kokkos::
  DualView<int*, LMPDeviceType::array_layout, LMPDeviceType> tdual_int_64;
typedef tdual_int_64::t_dev t_int_64;
typedef tdual_int_64::t_dev_um t_int_64_um;

};

#ifdef LMP_KOKKOS_GPU
template <>
struct FFTArrayTypes<LMPHostType> {

//Kspace

typedef Kokkos::
  DualView<FFT_SCALAR*, Kokkos::LayoutRight, LMPDeviceType> tdual_FFT_SCALAR_1d;
typedef tdual_FFT_SCALAR_1d::t_host t_FFT_SCALAR_1d;
typedef tdual_FFT_SCALAR_1d::t_host_um t_FFT_SCALAR_1d_um;

typedef Kokkos::DualView<FFT_SCALAR**,Kokkos::LayoutRight,LMPDeviceType> tdual_FFT_SCALAR_2d;
typedef tdual_FFT_SCALAR_2d::t_host t_FFT_SCALAR_2d;

typedef Kokkos::DualView<FFT_SCALAR**[3],Kokkos::LayoutRight,LMPDeviceType> tdual_FFT_SCALAR_2d_3;
typedef tdual_FFT_SCALAR_2d_3::t_host t_FFT_SCALAR_2d_3;

typedef Kokkos::DualView<FFT_SCALAR***,Kokkos::LayoutRight,LMPDeviceType> tdual_FFT_SCALAR_3d;
typedef tdual_FFT_SCALAR_3d::t_host t_FFT_SCALAR_3d;

typedef Kokkos::
  DualView<FFT_KOKKOS_DATA*, Kokkos::LayoutRight, LMPDeviceType> tdual_FFT_DATA_1d;
typedef tdual_FFT_DATA_1d::t_host t_FFT_DATA_1d;
typedef tdual_FFT_DATA_1d::t_host_um t_FFT_DATA_1d_um;

typedef Kokkos::
  DualView<int*, LMPDeviceType::array_layout, LMPDeviceType> tdual_int_64;
typedef tdual_int_64::t_host t_int_64;
typedef tdual_int_64::t_host_um t_int_64_um;

};
#endif

typedef struct FFTArrayTypes<LMPDeviceType> FFT_DAT;
typedef struct FFTArrayTypes<LMPHostType> FFT_HAT;


#if defined(FFT_KOKKOS_KISS)
#include "kissfft_kokkos.h" // uses t_FFT_DATA_1d, needs to come last
#endif


#endif
