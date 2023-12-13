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

#include "lmpfftsettings.h"

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

#endif
