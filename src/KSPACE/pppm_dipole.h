/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef KSPACE_CLASS

KSpaceStyle(pppm/dipole,PPPMDipole)

#else

#ifndef LMP_PPPM_DIPOLE_H
#define LMP_PPPM_DIPOLE_H

#include "pppm.h"

namespace LAMMPS_NS {

class PPPMDipole : public PPPM {
 public:
  PPPMDipole(class LAMMPS *, int, char **);
  virtual ~PPPMDipole();
  void init();
  void setup();
  void setup_grid();
  void compute(int, int);
  int timing_1d(int, double &);
  int timing_3d(int, double &);
  double memory_usage();

 protected:
  void set_grid_global();
  double newton_raphson_f();

  void allocate();
  void allocate_peratom();
  void deallocate();
  void deallocate_peratom();
  void compute_gf_denom();

  void slabcorr();

  // grid communication

  void pack_forward(int, FFT_SCALAR *, int, int *);
  void unpack_forward(int, FFT_SCALAR *, int, int *);
  void pack_reverse(int, FFT_SCALAR *, int, int *);
  void unpack_reverse(int, FFT_SCALAR *, int, int *);

  // dipole

  FFT_SCALAR ***densityx_brick_dipole,***densityy_brick_dipole,***densityz_brick_dipole;
  FFT_SCALAR ***vdxx_brick_dipole,***vdyy_brick_dipole,***vdzz_brick_dipole;
  FFT_SCALAR ***vdxy_brick_dipole,***vdxz_brick_dipole,***vdyz_brick_dipole;
  FFT_SCALAR ***ux_brick_dipole,***uy_brick_dipole,***uz_brick_dipole;
  FFT_SCALAR ***v0x_brick_dipole,***v1x_brick_dipole,***v2x_brick_dipole;
  FFT_SCALAR ***v3x_brick_dipole,***v4x_brick_dipole,***v5x_brick_dipole;
  FFT_SCALAR ***v0y_brick_dipole,***v1y_brick_dipole,***v2y_brick_dipole;
  FFT_SCALAR ***v3y_brick_dipole,***v4y_brick_dipole,***v5y_brick_dipole;
  FFT_SCALAR ***v0z_brick_dipole,***v1z_brick_dipole,***v2z_brick_dipole;
  FFT_SCALAR ***v3z_brick_dipole,***v4z_brick_dipole,***v5z_brick_dipole;
  FFT_SCALAR *work3,*work4;
  FFT_SCALAR *densityx_fft_dipole,*densityy_fft_dipole,*densityz_fft_dipole;
  class GridComm *cg_dipole;
  class GridComm *cg_peratom_dipole;
  int only_dipole_flag;
  double musum,musqsum,mu2;
  double find_gewald_dipole(double, double, bigint, double, double);
  double newton_raphson_f_dipole(double, double, bigint, double, double);
  double derivf_dipole(double, double, bigint, double, double);
  double compute_df_kspace_dipole();
  double compute_qopt_dipole();
  void compute_gf_dipole();
  void make_rho_dipole();
  void brick2fft_dipole();
  void poisson_ik_dipole();
  void poisson_peratom_dipole();
  void fieldforce_ik_dipole();
  void fieldforce_peratom_dipole();
  double final_accuracy_dipole();
  void musum_musq();

};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
