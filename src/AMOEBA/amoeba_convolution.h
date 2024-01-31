/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_AMOEBA_CONVOLUTION_H
#define LMP_AMOEBA_CONVOLUTION_H

#include "pointers.h"
#include "lmpfftsettings.h"

namespace LAMMPS_NS {

class AmoebaConvolution : protected Pointers {
 public:
  int nx, ny, nz;
  int order;
  int nfft_owned;    // owned grid points in FFT decomp
  int nxlo_in, nxhi_in, nylo_in, nyhi_in, nzlo_in, nzhi_in;
  int nxlo_out, nxhi_out, nylo_out, nyhi_out, nzlo_out, nzhi_out;
  int nxlo_fft, nxhi_fft, nylo_fft, nyhi_fft, nzlo_fft, nzhi_fft;
  bigint nfft_global;          // nx * ny * nz
  FFT_SCALAR *grid_brick_start;    // lower left corner of (c)grid_brick data

  AmoebaConvolution(class LAMMPS *, class Pair *, int, int, int, int, int);
  ~AmoebaConvolution();
  void reset_grid();
  void *zero();
  FFT_SCALAR *pre_convolution();
  void *post_convolution();

  double time_fft;

 protected:
  int which;                   // caller name for convolution being performed
  int flag3d;                  // 1 if using 3d grid_brick, 0 for 4d cgrid_brick
  int nbrick_owned;            // owned grid points in brick decomp
  int nbrick_ghosts;           // owned + ghost brick grid points
  int ngrid_either;            // max of nbrick_owned or nfft_owned

  class Pair *amoeba;
  class FFT3d *fft1, *fft2;
  class Grid3d *gc;
  class Remap *remap;

  FFT_SCALAR ***grid_brick;      // 3d real brick grid with ghosts
  FFT_SCALAR ****cgrid_brick;    // 4d complex brick grid with ghosts

  FFT_SCALAR *grid_fft;    // 3d FFT grid as 1d vector
  FFT_SCALAR *cfft;        // 3d complex FFT grid as 1d vector

  FFT_SCALAR *gc_buf1, *gc_buf2;    // buffers for GridComm
  FFT_SCALAR *remap_buf;            // buffer for Remap

  void allocate_grid();
  void deallocate_grid();
  void *zero_3d();
  void *zero_4d();
  FFT_SCALAR *pre_convolution_3d();
  virtual FFT_SCALAR *pre_convolution_4d();
  void *post_convolution_3d();
  virtual void *post_convolution_4d();
  void procs2grid2d(int, int, int, int &, int &);

  // DEBUG

  void debug_scalar(int, const char *);
  void debug_file(int, const char *);
};
}    // namespace LAMMPS_NS
#endif
