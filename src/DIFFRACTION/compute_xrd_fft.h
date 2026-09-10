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
  double kb_norm[3];     // 1/I0(beta), normalizes the window to a peak of one

  // the FFT mesh spans the "diffraction cell", the exact period of the phase
  // factor of node index d.  mesh_vec[d] projects a Cartesian coordinate onto
  // that index and scales it to (unfolded) mesh units, so it is the d-th
  // reciprocal lattice vector times the mesh dimension.  for an orthogonal cell
  // only component d is nonzero and this is a plain multiplication.

  double mesh_vec[3][3];

  // decomposition of the mesh into bricks, one per rank, also used as the
  // FFT3d in/out layout.  every rank's contribution has to reach the rank that
  // owns the mesh it lands on, so a decomposition into slabs would funnel the
  // whole reduction into as many ranks as the mesh has planes; with many more
  // ranks than that the owners cannot absorb it.  The grid falls back to slabs
  // when there are no more ranks than planes, which keeps the reduction a
  // single MPI call in that case.

  int pgrid[3];         // ranks along each dimension of the mesh
  int pme[3];           // position of this rank in that grid
  int fftlo[3];         // first mesh point of this rank's brick
  int ffthi[3];         // last mesh point, or one less than fftlo when empty
  int fftn[3];          // brick dimensions
  int is_slab;          // 1 when the grid is 1 x 1 x nprocs
  int nfft;             // local mesh points, = fftn[0]*fftn[1]*fftn[2]
  MPI_Comm fft_comm;    // ranks owning a brick; MPI_COMM_NULL on the others

  // the part of the mesh this rank spreads into.  the map from a Cartesian
  // coordinate to a mesh index is affine and then folded, so the image of this
  // rank's atoms is one interval per dimension, contiguous before the fold and
  // wrapping at most once after it.  foot_lo is where that interval starts and
  // foot_n how long it is, in mesh points; foot_n reaches nmesh[d] when the
  // atoms wrap the whole dimension, which puts a full copy of the mesh back on
  // every rank as before.

  int foot_lo[3], foot_n[3];
  bigint nfoot;               // foot_n[0]*foot_n[1]*foot_n[2]
  int *foot_all;              // 6 numbers per rank, the table above gathered
  int all_full;               // 1 when every rank holds the whole mesh
  int *recvcounts;            // MPI_Reduce_scatter block sizes, used when it does

  FFT_SCALAR *density_own;     // this rank's footprint, spread into
  FFT_SCALAR *density_slab;    // this rank's brick of the mesh after the exchange
  FFT_SCALAR *work1;           // complex FFT buffer, 2*nfft
  FFT_SCALAR *sendbuf;         // packed pieces of the footprint, by destination
  FFT_SCALAR *recvbuf;         // incoming piece from one source
  bigint maxsend, maxrecv;     // allocated lengths of the two buffers above

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
  void set_footprint();
  void set_pgrid();

  // a derived style may keep the mesh and the atoms somewhere else, so the
  // parts of the calculation that touch either of them are replaceable

  virtual void allocate_mesh();
  virtual void grow_density_own();
  virtual int minmax_u(double *, double *);
  virtual void refresh_scaling();
  virtual void deallocate();
  virtual void bucket_atoms();
  virtual void fold_reduce(int);
  static int segments(int, int, int, int, int, int (*)[3]);
  int qgrid_lo(int q, int d) const {
    int c = (d == 0) ? (q % pgrid[0]) : ((d == 1) ? ((q/pgrid[0]) % pgrid[1])
                                                  : (q/(pgrid[0]*pgrid[1])));
    return c*nmesh[d]/pgrid[d];
  }
  int qgrid_hi(int q, int d) const {
    int c = (d == 0) ? (q % pgrid[0]) : ((d == 1) ? ((q/pgrid[0]) % pgrid[1])
                                                  : (q/(pgrid[0]*pgrid[1])));
    return (c+1)*nmesh[d]/pgrid[d] - 1;
  }
  bigint brick_count(int);
  void pack_brick(int, FFT_SCALAR *);
  void unpack_brick(const FFT_SCALAR *, int);
  virtual void spread(int);
  double kb_window(int, int, int);
  static int factorable(int);
  double bessel_i0(double);
  void build_window();

  // 1/k^2 for the series in bessel_i0().  the division belongs to the
  // dependent chain of that loop and costs several times the rest of an
  // iteration, so it is done once here instead.

  static constexpr int BESSEL_TERMS = 128;
  double inv_ksq[BESSEL_TERMS];

  // Chebyshev expansion of the spreading window.  Evaluating the window
  // directly costs a Bessel function per stencil point per dimension of every
  // atom, which is several times the cost of the spreading it feeds.  As a
  // function of the offset of an atom from its nearest mesh point the window
  // is entire, being a power series in that offset squared, so a Chebyshev
  // expansion converges to round-off within a few terms.  Sixteen reaches it
  // for every stencil width.

  static constexpr int KB_NCHEB = 16;
  double *kb_cheb;    // 3*order*KB_NCHEB coefficients, [dim][stencil point][i]
  double *kb_tcheb;   // the KB_NCHEB Chebyshev polynomials at one offset

  // scratch for spread(), sized once rather than per element per invocation

  double *rho0, *rho1, *rho2;
  int *mx, *my, *mz;

  // the group atoms grouped by element, so that spreading one element does not
  // walk over the atoms of all the others.  the order within an element stays
  // the order of the atoms themselves, so the sum is the same one as before.

  int *slot_atoms;    // group atoms, ordered by element
  int *slot_start;    // where each element starts in slot_atoms, nslot+1 long
  int maxbucket;      // allocated length of slot_atoms

  // scratch for fold_reduce(), likewise

  int *sstart, *scount, *destlist;
  MPI_Request *requests;
};

}    // namespace LAMMPS_NS

#endif
#endif
