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

#ifndef LMP_FFT3D_KOKKOS_H
#define LMP_FFT3D_KOKKOS_H

#include "pointers.h"
#include "kokkos_type.h"
#include "remap_kokkos.h"
#include "fftdata_kokkos.h"

// with KOKKOS in CUDA mode we can only have
// CUFFT or KISSFFT, thus undefine all other
// FFTs here, since they may be valid in fft3d.cpp

#if defined(KOKKOS_ENABLE_CUDA)
# if defined(FFT_FFTW)
#  undef FFT_FFTW
# endif
# if defined(FFT_FFTW3)
#  undef FFT_FFTW3
# endif
# if defined(FFT_MKL)
#  undef FFT_MKL
# endif
# if !defined(FFT_CUFFT) && !defined(FFT_KISSFFT)
#  define FFT_KISSFFT
# endif
#else
# if defined(FFT_CUFFT)
#  error "Must enable CUDA with KOKKOS to use -DFFT_CUFFT"
# endif
#endif

// if user set FFTW, it means FFTW3

#ifdef FFT_FFTW
#define FFT_FFTW3
#endif


#ifdef FFT_FFTW_THREADS
#if !defined(FFT_FFTW3)
#error "Must use -DFFT_FFTW3 with -DFFT_FFTW_THREADS"
#endif
#endif

#if defined(FFT_MKL)
  #include "mkl_dfti.h"
  #if defined(FFT_SINGLE)
    typedef float _Complex FFT_DATA;
    #define FFT_MKL_PREC DFTI_SINGLE
  #else
    typedef double _Complex FFT_DATA;
    #define FFT_MKL_PREC DFTI_DOUBLE
  #endif
#elif defined(FFT_FFTW3)
  #include "fftw3.h"
  #if defined(FFT_SINGLE)
    typedef fftwf_complex FFT_DATA;
    #define FFTW_API(function)  fftwf_ ## function
  #else
    typedef fftw_complex FFT_DATA;
    #define FFTW_API(function) fftw_ ## function
  #endif
#elif defined(FFT_CUFFT)
  #include "cufft.h"
  #if defined(FFT_SINGLE)
    #define cufftExec cufftExecC2C
    #define CUFFT_TYPE CUFFT_C2C
    typedef cufftComplex FFT_DATA;
  #else
    #define cufftExec cufftExecZ2Z
    #define CUFFT_TYPE CUFFT_Z2Z
    typedef cufftDoubleComplex FFT_DATA;
  #endif
#else
  #include "kissfft_kokkos.h"
  #ifndef FFT_KISSFFT
  #define FFT_KISSFFT
  #endif
#endif

namespace LAMMPS_NS {

// -------------------------------------------------------------------------

// plan for how to perform a 3d FFT

template<class DeviceType>
struct fft_plan_3d_kokkos {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  struct remap_plan_3d_kokkos<DeviceType> *pre_plan;       // remap from input -> 1st FFTs
  struct remap_plan_3d_kokkos<DeviceType> *mid1_plan;      // remap from 1st -> 2nd FFTs
  struct remap_plan_3d_kokkos<DeviceType> *mid2_plan;      // remap from 2nd -> 3rd FFTs
  struct remap_plan_3d_kokkos<DeviceType> *post_plan;      // remap from 3rd FFTs -> output
  typename AT::t_FFT_DATA_1d d_copy;                   // memory for remap results (if needed)
  typename AT::t_FFT_DATA_1d d_scratch;                // scratch space for remaps
  int total1,total2,total3;         // # of 1st,2nd,3rd FFTs (times length)
  int length1,length2,length3;      // length of 1st,2nd,3rd FFTs
  int pre_target;                   // where to put remap results
  int mid1_target,mid2_target;
  int scaled;                       // whether to scale FFT results
  int normnum;                      // # of values to rescale
  double norm;                      // normalization factor for rescaling

                                    // system specific 1d FFT info
#if defined(FFT_MKL)
  DFTI_DESCRIPTOR *handle_fast;
  DFTI_DESCRIPTOR *handle_mid;
  DFTI_DESCRIPTOR *handle_slow;
#elif defined(FFT_FFTW3)
  FFTW_API(plan) plan_fast_forward;
  FFTW_API(plan) plan_fast_backward;
  FFTW_API(plan) plan_mid_forward;
  FFTW_API(plan) plan_mid_backward;
  FFTW_API(plan) plan_slow_forward;
  FFTW_API(plan) plan_slow_backward;
#elif defined(FFT_CUFFT)
  cufftHandle plan_fast;
  cufftHandle plan_mid;
  cufftHandle plan_slow;
#else
  kiss_fft_state_kokkos<DeviceType> cfg_fast_forward;
  kiss_fft_state_kokkos<DeviceType> cfg_fast_backward;
  kiss_fft_state_kokkos<DeviceType> cfg_mid_forward;
  kiss_fft_state_kokkos<DeviceType> cfg_mid_backward;
  kiss_fft_state_kokkos<DeviceType> cfg_slow_forward;
  kiss_fft_state_kokkos<DeviceType> cfg_slow_backward;
#endif
};

template<class DeviceType>
class FFT3dKokkos : protected Pointers {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FFT3dKokkos(class LAMMPS *, MPI_Comm,
        int,int,int,int,int,int,int,int,int,int,int,int,int,int,int,
        int,int,int *,int);
  ~FFT3dKokkos();
  void compute(typename AT::t_FFT_SCALAR_1d, typename AT::t_FFT_SCALAR_1d, int);
  void timing1d(typename AT::t_FFT_SCALAR_1d, int, int);

 private:
  struct fft_plan_3d_kokkos<DeviceType> *plan;
  RemapKokkos<DeviceType> *remapKK;

#ifdef FFT_KISSFFT
  KissFFTKokkos<DeviceType> *kissfftKK;
#endif

  void fft_3d_kokkos(typename AT::t_FFT_DATA_1d, typename AT::t_FFT_DATA_1d, int, struct fft_plan_3d_kokkos<DeviceType> *);

  struct fft_plan_3d_kokkos<DeviceType> *fft_3d_create_plan_kokkos(MPI_Comm, int, int, int,
                                         int, int, int, int, int,
                                         int, int, int, int, int, int, int,
                                         int, int, int *, int, int);

  void fft_3d_destroy_plan_kokkos(struct fft_plan_3d_kokkos<DeviceType> *);

  void fft_3d_1d_only_kokkos(typename AT::t_FFT_DATA_1d, int, int, struct fft_plan_3d_kokkos<DeviceType> *);

  void bifactor(int, int *, int *);
};

}

#endif

/* ERROR/WARNING messages:

E: Could not create 3d FFT plan

The FFT setup for the PPPM solver failed, typically due
to lack of memory.  This is an unusual error.  Check the
size of the FFT grid you are requesting.

E: Cannot use the FFTW library with Kokkos CUDA on GPUs 

Kokkos CUDA doesn't support using the FFTW library to calculate FFTs for 
PPPM on GPUs. 

E: Cannot use the cuFFT library with Kokkos CUDA on the host CPUs

Kokkos CUDA doesn't support using the cuFFT library to calculate FFTs 
for PPPM on the host CPUs, use KISS FFT instead. 

*/
