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

#define MAX(A,B) ((A) > (B) ? (A) : (B))


// data types for 2d/3d FFTs

#ifndef FFT_DATA_KOKKOS_H
#define FFT_DATA_KOKKOS_H

// User-settable FFT precision

// FFT_PRECISION = 1 is single-precision complex (4-byte real, 4-byte imag)
// FFT_PRECISION = 2 is double-precision complex (8-byte real, 8-byte imag)

#ifdef FFT_SINGLE
#define FFT_PRECISION 1
#define MPI_FFT_SCALAR MPI_FLOAT
typedef float FFT_SCALAR;
#else
#define FFT_PRECISION 2
#define MPI_FFT_SCALAR MPI_DOUBLE
typedef double FFT_SCALAR;
#endif

// -------------------------------------------------------------------------

// Data types for single-precision complex

#if FFT_PRECISION == 1
#elif FFT_PRECISION == 2
#else
#error "FFT_PRECISION needs to be either 1 (=single) or 2 (=double)"
#endif

#endif
