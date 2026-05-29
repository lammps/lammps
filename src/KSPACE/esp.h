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

#ifdef KSPACE_CLASS
// clang-format off
KSpaceStyle(esp,ESP);
// clang-format on
#else

#ifndef LMP_ESP_H
#define LMP_ESP_H

#include "kspace.h"
#include "lmpfftsettings.h"    // IWYU pragma: export
#include "math_const.h"

#include <cmath>

using namespace LAMMPS_NS;
using MathConst::MY_PI;

namespace LAMMPS_NS {

class ESP : public KSpace {
 public:
  ESP(class LAMMPS *);
  ~ESP() override;
  void settings(int, char **) override;
  void init() override;
  void NewFunction();
  void setup() override;
  void reset_grid() override;
  void compute(int, int) override;
  int timing_1d(int, double &) override;
  int timing_3d(int, double &) override;
  double memory_usage() override;
  void build_table(double);
  int estimate_order(double);
  void compute_group_group(int, int, int) override;

  // grid communication

  void pack_forward_grid(int, void *, int, int *) override;
  void unpack_forward_grid(int, void *, int, int *) override;
  void pack_reverse_grid(int, void *, int, int *) override;
  void unpack_reverse_grid(int, void *, int, int *) override;

 protected:
  int me, nprocs;
  int nfactors;
  int *factors;
  double cutoff;
  double volume;
  double delxinv, delyinv, delzinv, delvolinv;
  double h_x, h_y, h_z;
  double shift, shiftone, shiftatom_lo, shiftatom_hi;
  int peratom_allocate_flag;

  int nxlo_in, nylo_in, nzlo_in, nxhi_in, nyhi_in, nzhi_in;
  int nxlo_out, nylo_out, nzlo_out, nxhi_out, nyhi_out, nzhi_out;
  int nxlo_ghost, nxhi_ghost, nylo_ghost, nyhi_ghost, nzlo_ghost, nzhi_ghost;
  int nxlo_fft, nylo_fft, nzlo_fft, nxhi_fft, nyhi_fft, nzhi_fft;
  int nlower, nupper;
  int ngrid, nfft_brick, nfft, nfft_both;

  FFT_SCALAR ***density_brick;
  FFT_SCALAR ***vdx_brick, ***vdy_brick, ***vdz_brick;
  FFT_SCALAR ***u_brick;
  FFT_SCALAR ***v0_brick, ***v1_brick, ***v2_brick;
  FFT_SCALAR ***v3_brick, ***v4_brick, ***v5_brick;
  double *greensfn, *greensfn2;
  double **vg, **vg2;
  double *fkx, *fky, *fkz;
  FFT_SCALAR *density_fft;
  FFT_SCALAR *work1, *work2;

  FFT_SCALAR **rho1d, **rho_coeff, **drho1d,
      **drho_coeff;    // coefficients for the table of spreading function
  double *sf_precoeff1, *sf_precoeff2, *sf_precoeff3;
  double *sf_precoeff4, *sf_precoeff5, *sf_precoeff6;
  double sf_coeff[6];    // coefficients for calculating ad self-forces

  // FFTs and grid communication

  class FFT3d *fft1, *fft2;
  class Remap *remap;
  class Grid3d *gc;

  FFT_SCALAR *gc_buf1, *gc_buf2;
  int ngc_buf1, ngc_buf2, npergrid;

  // group-group interactions

  int group_allocate_flag;
  FFT_SCALAR ***density_A_brick, ***density_B_brick;
  FFT_SCALAR *density_A_fft, *density_B_fft;

  int **part2grid;    // storage for particle -> grid mapping
  int nmax;

  double *boxlo;
  // TIP4P settings
  int typeH, typeO;    // atom types of TIP4P water H and O atoms
  double qdist;        // distance from O site to negative charge
  double alpha;        // geometric factor

  virtual void set_grid_global();
  virtual void set_grid_local();

  virtual void allocate();
  virtual void allocate_peratom();
  virtual void deallocate();
  virtual void deallocate_peratom();
  int factorable(int);
  virtual void compute_gf_ik();
  virtual void compute_gf_ad();
  void compute_sf_precoeff();

  virtual void particle_map();
  virtual void make_rho();
  virtual void brick2fft();

  virtual void poisson();
  virtual void poisson_ik();
  virtual void poisson_ad();

  virtual void fieldforce();
  virtual void fieldforce_ik();
  virtual void fieldforce_ad();

  virtual void poisson_peratom();
  virtual void fieldforce_peratom();
  void procs2grid2d(int, int, int, int *, int *);
  void compute_rho1d(const FFT_SCALAR &, const FFT_SCALAR &, const FFT_SCALAR &);
  void compute_drho1d(const FFT_SCALAR &, const FFT_SCALAR &, const FFT_SCALAR &);
  void compute_rho_coeff();
  virtual void slabcorr();

  // triclinic

  int triclinic;    // domain settings, orthog or triclinic
  void setup_triclinic();
  void compute_gf_ik_triclinic();
  void poisson_ik_triclinic();
  void poisson_groups_triclinic();

  // group-group interactions

  virtual void allocate_groups();
  virtual void deallocate_groups();
  virtual void make_rho_groups(int, int, int);
  virtual void poisson_groups(int);
  virtual void slabcorr_groups(int, int, int);

  double poly_horner(const double x, const double *coeff, const int n) const
  {
    // coeff[0] + coeff[1] x + ... + coeff[n-1] x^(n-1)
    double p = coeff[n - 1];
    for (int i = n - 2; i >= 0; --i) p = p * x + coeff[i];
    return p;
  }

  void poly_and_deriv_horner(const double x, const double *coeff, const int n, double &p,
                             double &dp) const
  {
    // p(x) and dp/dx, Horner form
    p = coeff[n - 1];
    dp = 0.0;
    for (int i = n - 2; i >= 0; --i) {
      dp = dp * x + p;
      p = p * x + coeff[i];
    }
  }

  [[nodiscard]] double spreading_weight2_from_t(const double t) const
  {
    // t = (order * h / 2) * |q| / spreading_select_c
    // returns ( (order/2 * poly(2t-1))^2 ), or 0 if t>1
    if (t > 1.0) return 0.0;
    const double x = 2.0 * t - 1.0;
    // apply Horner rule for polynomial evaluation
    const double appx = poly_horner(x, fourier_spread_poly_coeff, fourier_spreading_order);
    const double w = 0.5 * order * appx;
    return w * w;
  }

  [[nodiscard]] double spreading_weight2_from_abs_index(const int abs_index,
                                                        const double scale) const
  {
    // integer-form helper: t = scale * abs_index
    return spreading_weight2_from_t(scale * (double) abs_index);
  }

  [[nodiscard]] double gf_denom_psw(const double &kx, const double &ky, const double &kz,
                                    const double &hx, const double &hy, const double &hz) const
  {
    int Nmax = (differentiation_flag == 0) ? 2 : 0;

    const double stepx = 2.0 * MY_PI / hx;
    const double stepy = 2.0 * MY_PI / hy;
    const double stepz = 2.0 * MY_PI / hz;

    // sum_{nx,ny,nz} wx*wy*wz = (sum wx)*(sum wy)*(sum wz)
    double sumx = 0.0, sumy = 0.0, sumz = 0.0;

    for (int nx = -Nmax; nx <= Nmax; ++nx) {
      const double qx = kx + stepx * (double) nx;
      const double t = (0.5 * order * hx * fabs(qx)) / spreading_select_c;
      sumx += spreading_weight2_from_t(t);
    }

    for (int ny = -Nmax; ny <= Nmax; ++ny) {
      const double qy = ky + stepy * (double) ny;
      const double t = (0.5 * order * hy * fabs(qy)) / spreading_select_c;
      sumy += spreading_weight2_from_t(t);
    }

    for (int nz = -Nmax; nz <= Nmax; ++nz) {
      const double qz = kz + stepz * (double) nz;
      const double t = (0.5 * order * hz * fabs(qz)) / spreading_select_c;
      sumz += spreading_weight2_from_t(t);
    }

    const double denom = sumx * sumy * sumz;
    return denom * denom;
  };
};

}    // namespace LAMMPS_NS

#endif
#endif
