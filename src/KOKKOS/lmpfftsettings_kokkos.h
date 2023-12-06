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

#ifndef LMP_FFT_KOKKOS_SETTINGS_H
#define LMP_FFT_KOKKOS_SETTINGS_H

// if user set FFTW, it means FFTW3

#ifdef FFT_KOKKOS_FFTW
#ifndef FFT_KOKKOS_FFTW3
#define FFT_KOKKOS_FFTW3
#endif
#endif

// set strings for library info output

#if defined(FFT_KOKKOS_FFTW3)
#define LMP_FFT_KOKKOS_LIB "FFTW3"
#elif defined(FFT_KOKKOS_MKL)
#define LMP_FFT_KOKKOS_LIB "MKL FFT"
#elif defined(FFT_KOKKOS_CUFFT)
#define LMP_FFT_KOKKOS_LIB "cuFFT"
#elif defined(FFT_KOKKOS_HIPFFT)
#define LMP_FFT_KOKKOS_LIB "hipFFT"
#else
#define LMP_FFT_KOKKOS_LIB "KISS FFT"
#endif

#ifdef FFT_KOKKOS_SINGLE
typedef float FFT_KOKKOS_SCALAR;
#define FFT_KOKKOS_PRECISION 1
#define LMP_FFT_KOKKOS_PREC "single"
#define MPI_FFT_KOKKOS_SCALAR MPI_FLOAT
#else

typedef double FFT_KOKKOS_SCALAR;
#define FFT_KOKKOS_PRECISION 2
#define LMP_FFT_KOKKOS_PREC "double"
#define MPI_FFT_KOKKOS_SCALAR MPI_DOUBLE
#endif

#endif
