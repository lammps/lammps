/* -*- c++ -*- ----------------------------------------------------------
 *
 *                    *** Smooth Mach Dynamics ***
 *
 * This file is part of the MACHDYN package for LAMMPS.
 * Copyright (2014) Georg C. Ganzenmueller, georg.ganzenmueller@emi.fhg.de
 * Fraunhofer Ernst-Mach Institute for High-Speed Dynamics, EMI,
 * Eckerstrasse 4, D-79104 Freiburg i.Br, Germany.
 *
 * ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
 LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
 https://www.lammps.org/, Sandia National Laboratories
 LAMMPS development team: developers@lammps.org

 Copyright (2003) Sandia Corporation.  Under the terms of Contract
 DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 certain rights in this software.  This software is distributed under
 the GNU General Public License.

 See the README file in the top-level LAMMPS directory.
 ------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(smd/ulsph,PairULSPH);
// clang-format on
#else

#ifndef LMP_ULSPH_H
#define LMP_ULSPH_H

#include "pair.h"
#include <Eigen/Eigen>

namespace LAMMPS_NS {

class PairULSPH : public Pair {
 public:
  PairULSPH(class LAMMPS *);
  ~PairULSPH() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  void init_style() override;
  void init_list(int, class NeighList *) override;
  double memory_usage() override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  void AssembleStressTensor();
  void *extract(const char *, int &) override;
  void PreCompute();
  void PreCompute_DensitySummation();
  double effective_shear_modulus(const Eigen::Matrix3d &d_dev,
                                 const Eigen::Matrix3d &deltaStressDev, const double dt,
                                 const int itype);

  Eigen::Vector3d ComputeHourglassForce(const int i, const int itype, const int j, const int jtype,
                                        const Eigen::Vector3d &dv, const Eigen::Vector3d &xij,
                                        const Eigen::Vector3d &g, const double c_ij,
                                        const double mu_ij, const double rho_ij);

 protected:
  double *c0_type;                    // reference speed of sound defined per particle type
  double *rho0;                       // reference mass density per type
  double *Q1;                         // linear artificial viscosity coeff
  int *eos, *viscosity, *strength;    // eos and strength material models
  double **artificial_pressure;       // true/false: use Monaghan's artificial pressure correction?
  double **artificial_stress;         // artificial stress amplitude

  double *onerad_dynamic, *onerad_frozen;
  double *maxrad_dynamic, *maxrad_frozen;

  void allocate();

  int nmax;    // max number of atoms on this proc
  int *numNeighs;
  Eigen::Matrix3d *K;
  double *shepardWeight, *c0, *rho;
  Eigen::Vector3d *smoothVel;
  Eigen::Matrix3d *stressTensor, *L, *F;

  double dtCFL;

 private:
  // enumerate EOSs. MUST BE IN THE RANGE [1000, 2000)
  enum {
    EOS_LINEAR = 1000,
    EOS_PERFECT_GAS = 1001,
    EOS_TAIT = 1002,
  };

  // enumerate physical viscosity models. MUST BE IN THE RANGE [2000, 3000)
  enum { VISCOSITY_NEWTON = 2000 };

  // enumerate strength models. MUST BE IN THE RANGE [3000, 4000)
  enum { STRENGTH_LINEAR = 3000, STRENGTH_LINEAR_PLASTIC = 3001 };

  // enumerate some quantitities and associate these with integer values such that they can be used for lookup in an array structure
  enum {
    NONE = 0,
    BULK_MODULUS = 1,
    HOURGLASS_CONTROL_AMPLITUDE = 2,
    EOS_TAIT_EXPONENT = 3,
    REFERENCE_SOUNDSPEED = 4,
    REFERENCE_DENSITY = 5,
    EOS_PERFECT_GAS_GAMMA = 6,
    SHEAR_MODULUS = 7,
    YIELD_STRENGTH = 8,
    YOUNGS_MODULUS = 9,
    POISSON_RATIO = 10,
    LAME_LAMBDA = 11,
    HEAT_CAPACITY = 12,
    M_MODULUS = 13,
    HARDENING_PARAMETER = 14,
    VISCOSITY_MU = 15,
    MAX_KEY_VALUE = 16
  };
  double **
      Lookup;    // holds per-type material parameters for the quantities defined in enum statement above.

  bool velocity_gradient_required;
  int updateFlag;    // indicates if any relative particle pair movement is significant compared to smoothing length

  bool density_summation, density_continuity, velocity_gradient, gradient_correction_flag;
  double *effm;
};

}    // namespace LAMMPS_NS

#endif
#endif
