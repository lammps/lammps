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
KSpaceStyle(pppm/disp,PPPMDisp);
// clang-format on
#else

#ifndef LMP_PPPM_DISP_H
#define LMP_PPPM_DISP_H

#include "kspace.h"
#include "lmpfftsettings.h"

namespace LAMMPS_NS {

static constexpr int EWALD_MAXORDER = 6;
static constexpr int EWALD_FUNCS = 4;

class PPPMDisp : public KSpace {
 public:
  PPPMDisp(class LAMMPS *);
  ~PPPMDisp() override;
  void init() override;
  void setup() override;
  void reset_grid() override;
  void settings(int, char **) override;
  void compute(int, int) override;
  int timing_1d(int, double &) override;
  int timing_3d(int, double &) override;
  double memory_usage() override;

 protected:
  int me, nprocs;
  int nfactors;
  int *factors;
  double csumij;
  double csum;
  double *csumi;    //needed as correction term for per atom calculations!
  double *cii;
  int csumflag;
  double cutoff, cutoff_lj;
  double volume;
  double *B;
  double virial_1[6], virial_6[6];
  double sf_coeff[6], sf_coeff_6[6];
  int peratom_allocate_flag;

  int nsplit;
  int nsplit_alloc;
  int function[EWALD_FUNCS];

  double delxinv, delyinv, delzinv, delvolinv;
  double delxinv_6, delyinv_6, delzinv_6, delvolinv_6;

  double shift, shiftone, shiftatom_lo, shiftatom_hi;
  int nxlo_in, nylo_in, nzlo_in, nxhi_in, nyhi_in, nzhi_in;
  int nxlo_out, nylo_out, nzlo_out, nxhi_out, nyhi_out, nzhi_out;
  int nxlo_fft, nylo_fft, nzlo_fft, nxhi_fft, nyhi_fft, nzhi_fft;
  int nlower, nupper;
  int ngrid, nfft_brick, nfft, nfft_both;

  double shift_6, shiftone_6, shiftatom_lo_6, shiftatom_hi_6;
  int nxlo_in_6, nylo_in_6, nzlo_in_6, nxhi_in_6, nyhi_in_6, nzhi_in_6;
  int nxlo_out_6, nylo_out_6, nzlo_out_6, nxhi_out_6, nyhi_out_6, nzhi_out_6;
  int nxlo_fft_6, nylo_fft_6, nzlo_fft_6, nxhi_fft_6, nyhi_fft_6, nzhi_fft_6;
  int nlower_6, nupper_6;
  int ngrid_6, nfft_brick_6, nfft_6, nfft_both_6;

  // the following variables are needed for every structure factor

  FFT_SCALAR ***density_brick;
  FFT_SCALAR ***vdx_brick, ***vdy_brick, ***vdz_brick;
  FFT_SCALAR *density_fft;
  FFT_SCALAR ***u_brick;
  FFT_SCALAR ***v0_brick, ***v1_brick, ***v2_brick, ***v3_brick, ***v4_brick, ***v5_brick;

  FFT_SCALAR ***density_brick_g;
  FFT_SCALAR ***vdx_brick_g, ***vdy_brick_g, ***vdz_brick_g;
  FFT_SCALAR *density_fft_g;
  FFT_SCALAR ***u_brick_g;
  FFT_SCALAR ***v0_brick_g, ***v1_brick_g, ***v2_brick_g, ***v3_brick_g, ***v4_brick_g,
      ***v5_brick_g;

  FFT_SCALAR ***density_brick_a0;
  FFT_SCALAR ***vdx_brick_a0, ***vdy_brick_a0, ***vdz_brick_a0;
  FFT_SCALAR *density_fft_a0;
  FFT_SCALAR ***u_brick_a0;
  FFT_SCALAR ***v0_brick_a0, ***v1_brick_a0, ***v2_brick_a0, ***v3_brick_a0, ***v4_brick_a0,
      ***v5_brick_a0;

  FFT_SCALAR ***density_brick_a1;
  FFT_SCALAR ***vdx_brick_a1, ***vdy_brick_a1, ***vdz_brick_a1;
  FFT_SCALAR *density_fft_a1;
  FFT_SCALAR ***u_brick_a1;
  FFT_SCALAR ***v0_brick_a1, ***v1_brick_a1, ***v2_brick_a1, ***v3_brick_a1, ***v4_brick_a1,
      ***v5_brick_a1;

  FFT_SCALAR ***density_brick_a2;
  FFT_SCALAR ***vdx_brick_a2, ***vdy_brick_a2, ***vdz_brick_a2;
  FFT_SCALAR *density_fft_a2;
  FFT_SCALAR ***u_brick_a2;
  FFT_SCALAR ***v0_brick_a2, ***v1_brick_a2, ***v2_brick_a2, ***v3_brick_a2, ***v4_brick_a2,
      ***v5_brick_a2;

  FFT_SCALAR ***density_brick_a3;
  FFT_SCALAR ***vdx_brick_a3, ***vdy_brick_a3, ***vdz_brick_a3;
  FFT_SCALAR *density_fft_a3;
  FFT_SCALAR ***u_brick_a3;
  FFT_SCALAR ***v0_brick_a3, ***v1_brick_a3, ***v2_brick_a3, ***v3_brick_a3, ***v4_brick_a3,
      ***v5_brick_a3;

  FFT_SCALAR ***density_brick_a4;
  FFT_SCALAR ***vdx_brick_a4, ***vdy_brick_a4, ***vdz_brick_a4;
  FFT_SCALAR *density_fft_a4;
  FFT_SCALAR ***u_brick_a4;
  FFT_SCALAR ***v0_brick_a4, ***v1_brick_a4, ***v2_brick_a4, ***v3_brick_a4, ***v4_brick_a4,
      ***v5_brick_a4;

  FFT_SCALAR ***density_brick_a5;
  FFT_SCALAR ***vdx_brick_a5, ***vdy_brick_a5, ***vdz_brick_a5;
  FFT_SCALAR *density_fft_a5;
  FFT_SCALAR ***u_brick_a5;
  FFT_SCALAR ***v0_brick_a5, ***v1_brick_a5, ***v2_brick_a5, ***v3_brick_a5, ***v4_brick_a5,
      ***v5_brick_a5;

  FFT_SCALAR ***density_brick_a6;
  FFT_SCALAR ***vdx_brick_a6, ***vdy_brick_a6, ***vdz_brick_a6;
  FFT_SCALAR *density_fft_a6;
  FFT_SCALAR ***u_brick_a6;
  FFT_SCALAR ***v0_brick_a6, ***v1_brick_a6, ***v2_brick_a6, ***v3_brick_a6, ***v4_brick_a6,
      ***v5_brick_a6;

  FFT_SCALAR ****density_brick_none;
  FFT_SCALAR ****vdx_brick_none, ****vdy_brick_none, ****vdz_brick_none;
  FFT_SCALAR **density_fft_none;
  FFT_SCALAR ****u_brick_none;
  FFT_SCALAR ****v0_brick_none, ****v1_brick_none, ****v2_brick_none, ****v3_brick_none,
      ****v4_brick_none, ****v5_brick_none;

  // needed for each interaction type

  double *greensfn;
  double **vg;
  double **vg2;

  double *greensfn_6;
  double **vg_6;
  double **vg2_6;

  double *fkx, *fky, *fkz;
  double *fkx2, *fky2, *fkz2;
  double *fkx_6, *fky_6, *fkz_6;
  double *fkx2_6, *fky2_6, *fkz2_6;
  double *gf_b;
  double *gf_b_6;

  double *sf_precoeff1, *sf_precoeff2, *sf_precoeff3, *sf_precoeff4, *sf_precoeff5, *sf_precoeff6;
  double *sf_precoeff1_6, *sf_precoeff2_6, *sf_precoeff3_6, *sf_precoeff4_6, *sf_precoeff5_6,
      *sf_precoeff6_6;
  FFT_SCALAR **rho1d, **rho_coeff;
  FFT_SCALAR **drho1d, **drho_coeff;
  FFT_SCALAR **rho1d_6, **rho_coeff_6;
  FFT_SCALAR **drho1d_6, **drho_coeff_6;
  FFT_SCALAR *work1, *work2;
  FFT_SCALAR *work1_6, *work2_6;

  class FFT3d *fft1, *fft2;
  class FFT3d *fft1_6, *fft2_6;
  class Remap *remap, *remap_6;
  class Grid3d *gc, *gc6;

  FFT_SCALAR *gc_buf1, *gc_buf2, *gc6_buf1, *gc6_buf2;
  int ngc_buf1, ngc_buf2, npergrid;
  int ngc6_buf1, ngc6_buf2, npergrid6;

  int **part2grid;    // storage for particle -> grid mapping
  int **part2grid_6;
  int nmax;

  int triclinic;    // domain settings, orthog or triclinic
  double *boxlo;

  // TIP4P settings

  int typeH, typeO;    // atom types of TIP4P water H and O atoms
  double qdist;        // distance from O site to negative charge
  double alpha;        // geometric factor

  void init_coeffs();
  int qr_alg(double **, double **, int);
  void hessenberg(double **, double **, int);
  void qr_tri(double **, double **, int);
  void mmult(double **, double **, double **, int);
  int check_convergence(double **, double **, double **, double **, double **, double **, int);

  void set_grid_global();
  void set_grid_global_6();
  void set_grid_local(int, int, int, int, double &, double &, double &, double &,
                      int &, int &, int &, int &, int &, int &, int &, int &);
  void set_init_g6();
  void set_n_pppm_6();

  void adjust_gewald();
  void adjust_gewald_6();

  double f();
  double derivf();
  double f_6();
  double derivf_6();
  double final_accuracy();
  void final_accuracy_6(double &, double &, double &);
  double lj_rspace_error();

  double compute_qopt();
  double compute_qopt_6();
  double compute_qopt_ik();
  double compute_qopt_ad();
  double compute_qopt_6_ik();
  double compute_qopt_6_ad();

  void calc_csum();

  virtual void allocate();
  virtual void allocate_peratom();
  virtual void deallocate();
  virtual void deallocate_peratom();

  int factorable(int);
  double rms(double, double, bigint, double, double **);
  double diffpr(double, double, double, double, double **);
  void compute_gf_denom(double *, int);
  double gf_denom(double, double, double, double *, int);

  void compute_sf_precoeff(int, int, int, int, int, int, int, int, int, int, double *, double *,
                           double *, double *, double *, double *);
  void compute_gf();
  void compute_sf_coeff();
  void compute_gf_6();
  void compute_sf_coeff_6();

  virtual void particle_map(double, double, double, double, int **, int, int, int, int, int, int,
                            int, int);
  virtual void particle_map_c(double, double, double, double, int **, int, int, int, int, int, int,
                              int, int);
  virtual void make_rho_c();
  virtual void make_rho_g();
  virtual void make_rho_a();
  virtual void make_rho_none();

  virtual void brick2fft(int, int, int, int, int, int, FFT_SCALAR ***, FFT_SCALAR *, FFT_SCALAR *,
                         LAMMPS_NS::Remap *);
  virtual void brick2fft_a();
  virtual void brick2fft_none();

  virtual void poisson_ik(FFT_SCALAR *, FFT_SCALAR *, FFT_SCALAR *, LAMMPS_NS::FFT3d *,
                          LAMMPS_NS::FFT3d *, int, int, int, int, int, int, int, int, int, int, int,
                          int, int, int, int, int, double &, double *, double *, double *, double *,
                          double *, double *, double *, FFT_SCALAR ***, FFT_SCALAR ***,
                          FFT_SCALAR ***, double *, double **, double **, FFT_SCALAR ***,
                          FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***,
                          FFT_SCALAR ***, FFT_SCALAR ***);

  virtual void poisson_ad(FFT_SCALAR *, FFT_SCALAR *, FFT_SCALAR *, LAMMPS_NS::FFT3d *,
                          LAMMPS_NS::FFT3d *, int, int, int, int, int, int, int, int, int, int, int,
                          int, int, int, int, int, double &, double *, double *, double **,
                          double **, FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***,
                          FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***);

  virtual void poisson_peratom(FFT_SCALAR *, FFT_SCALAR *, LAMMPS_NS::FFT3d *, double **, double **,
                               int, int, int, int, int, int, int, FFT_SCALAR ***, FFT_SCALAR ***,
                               FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***);
  virtual void poisson_2s_ik(FFT_SCALAR *, FFT_SCALAR *, FFT_SCALAR ***, FFT_SCALAR ***,
                             FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***,
                             FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***,
                             FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***,
                             FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***,
                             FFT_SCALAR ***, FFT_SCALAR ***);
  virtual void poisson_2s_ad(FFT_SCALAR *, FFT_SCALAR *, FFT_SCALAR ***, FFT_SCALAR ***,
                             FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***,
                             FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***,
                             FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***);

  virtual void poisson_2s_peratom(FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***,
                                  FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***,
                                  FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***);

  virtual void poisson_none_ad(int, int, FFT_SCALAR *, FFT_SCALAR *, FFT_SCALAR ***, FFT_SCALAR ***,
                               FFT_SCALAR ****, FFT_SCALAR ****, FFT_SCALAR ****, FFT_SCALAR ****,
                               FFT_SCALAR ****, FFT_SCALAR ****);
  virtual void poisson_none_ik(int, int, FFT_SCALAR *, FFT_SCALAR *, FFT_SCALAR ***, FFT_SCALAR ***,
                               FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***,
                               FFT_SCALAR ****, FFT_SCALAR ****, FFT_SCALAR ****, FFT_SCALAR ****,
                               FFT_SCALAR ****, FFT_SCALAR ****, FFT_SCALAR ****);
  virtual void poisson_none_peratom(int, int, FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***,
                                    FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***,
                                    FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***, FFT_SCALAR ***,
                                    FFT_SCALAR ***);

  virtual void fieldforce_c_ik();
  virtual void fieldforce_c_ad();
  virtual void fieldforce_c_peratom();
  virtual void fieldforce_g_ik();
  virtual void fieldforce_g_ad();
  virtual void fieldforce_g_peratom();
  virtual void fieldforce_a_ik();
  virtual void fieldforce_a_ad();
  virtual void fieldforce_a_peratom();
  virtual void fieldforce_none_ik();
  virtual void fieldforce_none_ad();
  virtual void fieldforce_none_peratom();
  void procs2grid2d(int, int, int, int *, int *);
  void compute_rho1d(const FFT_SCALAR &, const FFT_SCALAR &, const FFT_SCALAR &, int, FFT_SCALAR **,
                     FFT_SCALAR **);
  void compute_drho1d(const FFT_SCALAR &, const FFT_SCALAR &, const FFT_SCALAR &, int,
                      FFT_SCALAR **, FFT_SCALAR **);
  void compute_rho_coeff(FFT_SCALAR **, FFT_SCALAR **, int);
  virtual void slabcorr(int);

  // grid communication

  void pack_forward_grid(int, void *, int, int *) override;
  void unpack_forward_grid(int, void *, int, int *) override;
  void pack_reverse_grid(int, void *, int, int *) override;
  void unpack_reverse_grid(int, void *, int, int *) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
