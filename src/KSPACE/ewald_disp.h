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
KSpaceStyle(ewald/disp,EwaldDisp);
KSpaceStyle(ewald/disp/dipole,EwaldDisp);
// clang-format on
#else

#ifndef LMP_EWALD_DISP_H
#define LMP_EWALD_DISP_H

#include "ewald_const.h"
#include "kspace.h"

namespace LAMMPS_NS {

// total number of k-space sum channels: 1 (coulomb) + 1 (geometric 1/r^6)
// + 7 (arithmetic 1/r^6 sub-sums) + 1 (dipole)

static constexpr int EWALD_MAX_NSUMS = 10;

class EwaldDisp : public KSpace {
 public:
  EwaldDisp(class LAMMPS *);
  ~EwaldDisp() override;
  void init() override;
  void setup() override;
  void settings(int, char **) override;
  void compute(int, int) override;
  double memory_usage() override;

 private:
  double unitk[6];    // 2 pi times inverse box shape (domain->h_inv layout:
                      // diagonal in [0..2], yz/xz/xy tilt couplings in [3..5])
  int termflag[EwaldConst::EWALD_NTERMS];    // which k-space terms are enabled
  int nterms;                                // number of enabled terms
  int nsums;                                 // number of enabled sum channels
  int coeff_sums_flag;                       // 1 if coeff_sum[] is up-to-date
  int first_output;

  int kcount, kcount_max;    // number of k-vectors in use / allocated
  int kmax;                  // max k-index in any dimension
  int nevec, nevec_max;      // size of the per-atom e^(ikr) table
  int peratom_allocate_flag;
  int nmax;
  double gsqmx, q2, b2, M2;
  double *kenergy;    // energy prefactors: nterms per k-vector
  double *kvirial;    // virial prefactors: 6*nterms per k-vector
  double energy_self[EwaldConst::EWALD_NTERMS];
  double virial_self[EwaldConst::EWALD_NTERMS];
  double **energy_self_peratom;
  double **virial_self_peratom;
  struct cvector *ekr;    // per-atom powers e^(i*k*r), 2*kmax+1 per atom
  struct hvector *hvec;   // k-vectors, cartesian components
  struct kvector *kvec;   // k-vectors, integer indices shifted by kmax

  double mumurd2e, dielectric, *B, volume;

  // global sums of dispersion coefficients (x) and their squares (x2) per
  // sum channel: 0 = charge (unused), 1 = geometric 1/r^6, 2..8 = the seven
  // arithmetic 1/r^6 sub-sums, 9 = dipole (mu^2)

  struct Sum {
    double x, x2;
  } coeff_sum[EWALD_MAX_NSUMS];

  // structure factors: per-processor partial sums and the reduced totals,
  // nsums channels per k-vector, channel order coulomb, geometric 1/r^6,
  // arithmetic 1/r^6 (7 channels), dipole (enabled channels only)

  struct complex *sfac, *sfac_all;

  double rms(int, double, bigint, double, double, double);
  void reallocate();
  void allocate_peratom();
  void reallocate_atoms();
  void deallocate();
  void deallocate_peratom();
  void coeffs();
  void init_coeffs();
  void init_coeff_sums();
  void init_self();
  void init_self_peratom();
  void eik_dot_r();
  void compute_force();
  void compute_energy();
  void compute_energy_peratom();
  void compute_virial();
  void compute_virial_dipole();
  void compute_virial_peratom();
  void slabcorr();
  double NewtonSolve(double, double, bigint, double, double);
  double f(double, double, bigint, double, double);
  double derivf(double, double, bigint, double, double);
};

}    // namespace LAMMPS_NS

#endif
#endif
