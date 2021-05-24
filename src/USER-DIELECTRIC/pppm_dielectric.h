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

KSpaceStyle(pppm/dielectric,PPPMDielectric)

#else

#ifndef LMP_PPPM_DIELECTRIC_H
#define LMP_PPPM_DIELECTRIC_H

#include "pppm.h"

namespace LAMMPS_NS {

class PPPMDielectric : public PPPM {
 public:
  PPPMDielectric(class LAMMPS *);
  virtual ~PPPMDielectric();
  virtual void init();
  void setup_grid();
  virtual void compute(int, int);
  virtual double memory_usage();

  double** efield;
  double* phi;
  int potflag;   // 1/0 if per-atom electrostatic potential phi is needed

  void qsum_qsq();

 protected:
  FFT_SCALAR ***densityx_brick_dipole,***densityy_brick_dipole,***densityz_brick_dipole;
  FFT_SCALAR ***vdxx_brick_dipole,***vdyy_brick_dipole,***vdzz_brick_dipole;
  FFT_SCALAR ***vdxy_brick_dipole,***vdxz_brick_dipole,***vdyz_brick_dipole;
  FFT_SCALAR ***u_brick_dipole;
  FFT_SCALAR ***ux_brick_dipole,***uy_brick_dipole,***uz_brick_dipole;
  FFT_SCALAR *densityx_fft_dipole,*densityy_fft_dipole,*densityz_fft_dipole;
  FFT_SCALAR *work3,*work4;

  class GridComm *cg_mu;

  virtual void allocate();
  virtual void deallocate();
  void slabcorr();

  void fieldforce_ik();
  void fieldforce_ad();

  // grid communication

  virtual void pack_forward(int, FFT_SCALAR *, int, int *);
  virtual void unpack_forward(int, FFT_SCALAR *, int, int *);
  virtual void pack_reverse(int, FFT_SCALAR *, int, int *);
  virtual void unpack_reverse(int, FFT_SCALAR *, int, int *);

  // dipole

  int mu_flag;
  double musqsum,musum,mu2;
  void make_rho_dipole();
  void brick2fft_dipole();
  void poisson_ik_dipole();
  void fieldforce_ik_dipole();
  void musum_musq();

  class AtomVecDielectric* avec;

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
