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
KSpaceStyle(pppm,PPPM);
// clang-format on
#else

#ifndef LMP_PPPM_H
#define LMP_PPPM_H

#include "kspace.h"
#include "lmpfftsettings.h" // IWYU pragma: export

namespace LAMMPS_NS {

class PPPM : public KSpace {
 public:
  PPPM(class LAMMPS *);
  ~PPPM() override;
  void settings(int, char **) override;
  void init() override;
  void setup() override;
  void reset_grid() override;
  void compute(int, int) override;
  int timing_1d(int, double &) override;
  int timing_3d(int, double &) override;
  double memory_usage() override;

  void compute_group_group(int, int, int) override;

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
  double *greensfn;
  double **vg;
  double *fkx, *fky, *fkz;
  FFT_SCALAR *density_fft;
  FFT_SCALAR *work1, *work2;

  double *gf_b;
  FFT_SCALAR **rho1d, **rho_coeff, **drho1d, **drho_coeff;
  double *sf_precoeff1, *sf_precoeff2, *sf_precoeff3;
  double *sf_precoeff4, *sf_precoeff5, *sf_precoeff6;
  double sf_coeff[6];    // coefficients for calculating ad self-forces
  double **acons;

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
  void adjust_gewald();
  virtual double newton_raphson_f();
  double derivf();
  double final_accuracy();

  virtual void allocate();
  virtual void allocate_peratom();
  virtual void deallocate();
  virtual void deallocate_peratom();
  int factorable(int);
  virtual double compute_df_kspace();
  double estimate_ik_error(double, double, bigint);
  virtual double compute_qopt();
  virtual void compute_gf_denom();
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

  // grid communication

  void pack_forward_grid(int, void *, int, int *) override;
  void unpack_forward_grid(int, void *, int, int *) override;
  void pack_reverse_grid(int, void *, int, int *) override;
  void unpack_reverse_grid(int, void *, int, int *) override;

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

  /* ----------------------------------------------------------------------
   denominator for Hockney-Eastwood Green's function
     of x,y,z = sin(kx*deltax/2), etc

            inf                 n-1
   S(n,k) = Sum  W(k+pi*j)**2 = Sum b(l)*(z*z)**l
           j=-inf               l=0

          = -(z*z)**n /(2n-1)! * (d/dx)**(2n-1) cot(x)  at z = sin(x)
   gf_b = denominator expansion coeffs
------------------------------------------------------------------------- */

  inline double gf_denom(const double &x, const double &y, const double &z) const
  {
    double sx, sy, sz;
    sz = sy = sx = 0.0;
    for (int l = order - 1; l >= 0; l--) {
      sx = gf_b[l] + sx * x;
      sy = gf_b[l] + sy * y;
      sz = gf_b[l] + sz * z;
    }
    double s = sx * sy * sz;
    return s * s;
  };
};

}    // namespace LAMMPS_NS

#endif
#endif
