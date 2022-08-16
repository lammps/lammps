/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

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
#include "gridcomm.h"

using namespace LAMMPS_NS;

#define SCALE 0

enum {FORWARD,BACKWARD};

// External functions from GPU library

int amoeba_setup_fft(const int size, const int element_type);
int amoeba_compute_fft1d(void* in, void* out, const int mode);

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
  debug_scalar(GRIDBRICK_OUT,"PRE Convo / PRE GridComm");
#endif

  gc->reverse_comm(GridComm::PAIR,amoeba,2,sizeof(FFT_SCALAR),which,
                   gc_buf1,gc_buf2,MPI_FFT_SCALAR);

#if DEBUG_AMOEBA
  debug_scalar(GRIDBRICK_IN,"PRE Convo / POST GridComm");
  debug_file(GRIDBRICK_IN,"pre.convo.post.gridcomm");
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

  // perform forward FFT

  fft1->compute(cfft,cfft,FFT3d::FORWARD);

  //amoeba_compute_fft1d(cfft,cfft,FORWARD);

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
  fft2->compute(cfft,cfft,FFT3d::BACKWARD);

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
  debug_scalar(GRIDBRICK_IN,"POST Convo / PRE gridcomm");
  debug_file(GRIDBRICK_IN,"post.convo.pre.gridcomm");
#endif
  gc->forward_comm(GridComm::PAIR,amoeba,2,sizeof(FFT_SCALAR),which,
                   gc_buf1,gc_buf2,MPI_FFT_SCALAR);

  return (void *) cgrid_brick;
}
