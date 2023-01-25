/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "amoeba_convolution_gpu.h"
#include "comm.h"
#include "fft3d_wrap.h"
#include "remap_wrap.h"
#include "grid3d.h"

using namespace LAMMPS_NS;

// DEBUG

#define DEBUG_AMOEBA 0
#if DEBUG_AMOEBA
char *labels[7] =
  {(char *) "MPOLE_GRID", (char *) "POLAR_GRID",
   (char *) "POLAR_GRIDC", (char *) "DISP_GRID",
   (char *) "INDUCE_GRID", (char *) "INDUCE_GRIDC"};

enum{GRIDBRICK_OUT,GRIDBRICK_IN,FFT,CFFT1,CFFT2};
#endif
// END DEBUG

#define SCALE 0

//#define USE_AMOEBA_FFT
#ifdef USE_AMOEBA_FFT
// External functions from GPU library
int amoeba_setup_fft(const int size, const int numel, const int element_type);
int amoeba_compute_fft1d(void* in, void* out, const int numel, const int mode);
#endif

/* ----------------------------------------------------------------------
   partition an FFT grid across processors
   both for a brick and FFT x pencil decomposition
   nx,nz,nz = global FFT grid size
   order = size of stencil in each dimension that maps atoms to grid
   adapted from PPPM::set_grid_local()
------------------------------------------------------------------------- */

AmoebaConvolutionGPU::AmoebaConvolutionGPU(LAMMPS *lmp, Pair *pair,
                                     int nx_caller, int ny_caller, int nz_caller,
                                     int order_caller, int which_caller) :
  AmoebaConvolution(lmp, pair, nx_caller, ny_caller,  nz_caller, order_caller,
                    which_caller)
{

}

/* ----------------------------------------------------------------------
   perform pre-convolution grid operations for 4d cgrid_brick array
------------------------------------------------------------------------- */

FFT_SCALAR *AmoebaConvolutionGPU::pre_convolution_4d()
{
  int ix,iy,iz,n;

  // reverse comm for 4d brick grid + ghosts

#if DEBUG_AMOEBA
  debug_scalar(GRIDBRICK_OUT,"PRE Convo / PRE Grid3d");
#endif

  gc->reverse_comm(Grid3d::PAIR,amoeba,which,2,sizeof(FFT_SCALAR),
                   gc_buf1,gc_buf2,MPI_FFT_SCALAR);

#if DEBUG_AMOEBA
  debug_scalar(GRIDBRICK_IN,"PRE Convo / POST Grid3d");
  debug_file(GRIDBRICK_IN,"pre.convo.post.grid3d");
#endif
  // copy owned 4d brick grid values to FFT grid

  n = 0;
  for (iz = nzlo_in; iz <= nzhi_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        cfft[n++] = cgrid_brick[iz][iy][ix][0];
        cfft[n++] = cgrid_brick[iz][iy][ix][1];
      }

  // remap FFT grid from brick to x pencil partitioning
  // NOTE: could just setup FFT to start from brick decomp and skip remap

  remap->perform(cfft,cfft,remap_buf);

#if DEBUG_AMOEBA
  debug_scalar(FFT,"PRE Convo / POST Remap");
  debug_file(FFT,"pre.convo.post.remap");
#endif

  double time0,time1;

  MPI_Barrier(world);
  time0 = platform::walltime();

  // perform forward FFT

  #ifdef USE_AMOEBA_FFT
  amoeba_compute_fft1d(cfft,cfft,2*nfft_owned,FFT3d::FORWARD);
  #else
  fft1->compute(cfft,cfft,FFT3d::FORWARD);
  #endif

  time1 = platform::walltime();

  time_fft += time1 - time0;

  if (SCALE) {
    double scale = 1.0/nfft_global;
    for (int i = 0; i < 2*nfft_owned; i++) cfft[i] *= scale;
  }

#if DEBUG_AMOEBA
  debug_scalar(CFFT1,"PRE Convo / POST FFT");
  debug_file(CFFT1,"pre.convo.post.fft");
#endif
  return cfft;
}

/* ----------------------------------------------------------------------
   perform post-convolution grid operations for 4d cgrid_brick array
------------------------------------------------------------------------- */

void *AmoebaConvolutionGPU::post_convolution_4d()
{
  int ix,iy,iz,n;

  // perform backward FFT

#if DEBUG_AMOEBA
  debug_scalar(CFFT1,"POST Convo / PRE FFT");
  debug_file(CFFT1,"post.convo.pre.fft");
#endif

  double time0,time1;

  MPI_Barrier(world);
  time0 = platform::walltime();

  fft2->compute(cfft,cfft,FFT3d::BACKWARD);

  time1 = platform::walltime();

  time_fft += time1 - time0;

#if DEBUG_AMOEBA
  debug_scalar(CFFT2,"POST Convo / POST FFT");
  debug_file(CFFT2,"post.convo.post.fft");
#endif
  // copy 1d complex values into 4d complex grid

  n = 0;
  for (iz = nzlo_in; iz <= nzhi_in; iz++)
    for (iy = nylo_in; iy <= nyhi_in; iy++)
      for (ix = nxlo_in; ix <= nxhi_in; ix++) {
        cgrid_brick[iz][iy][ix][0] = cfft[n++];
        cgrid_brick[iz][iy][ix][1] = cfft[n++];
      }

  // forward comm to populate ghost grid values

#if DEBUG_AMOEBA
  debug_scalar(GRIDBRICK_IN,"POST Convo / PRE grid3d");
  debug_file(GRIDBRICK_IN,"post.convo.pre.grid3d");
#endif
  gc->forward_comm(Grid3d::PAIR,amoeba,which,2,sizeof(FFT_SCALAR),
                   gc_buf1,gc_buf2,MPI_FFT_SCALAR);

  return (void *) cgrid_brick;
}
