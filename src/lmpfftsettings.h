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

// common FFT library related defines and compilation settings

#ifndef LMP_FFT_SETTINGS_H
#define LMP_FFT_SETTINGS_H

// if a user sets FFTW, it means FFTW3

#ifdef FFT_FFTW
#ifndef FFT_FFTW3
#undef FFT_FFTW
#define FFT_FFTW3
#endif
#endif

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

// set strings for library info output

#if defined(FFT_HEFFTE)
#if defined(FFT_HEFFTE_FFTW)
#define LMP_FFT_LIB "HeFFTe(FFTW3)"
#elif defined(FFT_HEFFTE_MKL)
#define LMP_FFT_LIB "HeFFTe(MKL)"
#else
#define LMP_FFT_LIB "HeFFTe(builtin)"
#endif
#elif defined(FFT_FFTW3)
#define LMP_FFT_LIB "FFTW3"
#elif defined(FFT_MKL)
#define LMP_FFT_LIB "MKL FFT"
#elif defined(FFT_CUFFT)
#define LMP_FFT_LIB "cuFFT"
#elif defined(FFT_HIPFFT)
#define LMP_FFT_LIB "hipFFT"
#else
#define LMP_FFT_LIB "KISS FFT"
#endif

#ifdef LMP_KOKKOS

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
#endif

#ifdef FFT_SINGLE
typedef float FFT_SCALAR;
#define FFT_PRECISION 1
#define LMP_FFT_PREC "single"
#define MPI_FFT_SCALAR MPI_FLOAT
#else

typedef double FFT_SCALAR;
#define FFT_PRECISION 2
#define LMP_FFT_PREC "double"
#define MPI_FFT_SCALAR MPI_DOUBLE
#endif

#endif
