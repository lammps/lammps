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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(xrd/fft,ComputeXRDFFT);
// clang-format on
#else

#ifndef LMP_COMPUTE_XRD_FFT_H
#define LMP_COMPUTE_XRD_FFT_H

#include "compute_xrd.h"
#include "lmpfftsettings.h"

namespace LAMMPS_NS {

class ComputeXRDFFT : public ComputeXRD {
 public:
  ComputeXRDFFT(class LAMMPS *, int, char **);
  ~ComputeXRDFFT() override;
  void init() override;
  void compute_array() override;
  double memory_usage() override;

 protected:
  int nprocs;
  int setup_done;

  // Kaiser-Bessel spreading kernel and FFT mesh

  int order;             // stencil width in grid points (odd, >= 3)
  double oversample;     // FFT mesh oversampling factor
  int nlower, nupper;    // stencil extent, = -(order-1)/2 and (order-1)/2
  int nmesh[3];          // global FFT mesh dimensions
  double kb_beta[3];     // Kaiser-Bessel shape parameter, per dimension

  // the FFT mesh spans the "diffraction cell", the exact period of the phase
  // factor of node index d.  mesh_vec[d] projects a Cartesian coordinate onto
  // that index and scales it to (unfolded) mesh units, so it is the d-th
  // reciprocal lattice vector times the mesh dimension.  for an orthogonal cell
  // only component d is nonzero and this is a plain multiplication.

  double mesh_vec[3][3];

  // z-slab decomposition of the mesh, also used as the FFT3d in/out layout

  int nzlo_fft, nzhi_fft, nslab;
  int nfft;             // local mesh points, = nmesh[0]*nmesh[1]*nslab
  MPI_Comm fft_comm;    // ranks owning a slab; MPI_COMM_NULL on the others

  FFT_SCALAR *density_all;     // full mesh, spread into by every rank
  FFT_SCALAR *density_slab;    // this rank's z-slab after the reduction
  FFT_SCALAR *work1;           // complex FFT buffer, 2*nfft
  int *recvcounts;             // MPI_Reduce_scatter block sizes

  // atomic scattering factors are shared by all atom types of the same element,
  // so the transform is done once per distinct element rather than per type

  int nslot;
  int *slot_of_type;     // ntypes -> slot
  int *ztype_of_slot;    // slot -> row of ASFXRD

  // reciprocal lattice nodes whose FFT mode is owned by this rank

  int nown;
  int *own_row;         // row index in the output array
  int *own_idx;         // offset into the local FFT mesh
  double *own_deconv;   // 1 / (product of the three window transforms)
  double *own_lp;       // Lorentz-polarization factor, 1.0 when LP is off
  double *own_asf;      // nslot*nown scattering factors

  double *Fre, *Fim;    // per-row structure factor accumulators
  double *Iloc, *Iall;  // per-row intensity, before and after reduction

  class FFT3d *fft1;

  void set_grid();
  void setup_mesh();
  void refresh_scaling();
  void deallocate();
  void spread(int);
  double kb_window(int, int, double);
  static int factorable(int);
  static double bessel_i0(double);
};

}    // namespace LAMMPS_NS

#endif
#endif
