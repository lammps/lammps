/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
-------------------------------------------------------------------------
   BOCS written by: Nicholas J. H. Dunn and Michael R. DeLyser
   from The Pennsylvania State University
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(PRESSURE/BOCS,ComputePressureBocs);
// clang-format on
#else

#ifndef LMP_COMPUTE_PRESSURE_BOCS_H
#define LMP_COMPUTE_PRESSURE_BOCS_H

#include "compute.h"

namespace LAMMPS_NS {
// Enumerate the p_basis_type magic values to improve readability:
enum { BASIS_ANALYTIC, BASIS_LINEAR_SPLINE, BASIS_CUBIC_SPLINE };
// Enumerate the data file column names to improve readability
enum { VOLUME, PRESSURE_CORRECTION };
// Declare names for the number of columns in the splines data structure to improve readability
const int NUM_LINEAR_SPLINE_COLUMNS = 2;    // linear spline columns passed to compute
const int NUM_CUBIC_SPLINE_COLUMNS = 5;     // cubic spline columns passed to compute

// ComputePressure -> ComputePressureBocs MRD NJD
class ComputePressureBocs : public Compute {
 public:
  ComputePressureBocs(class LAMMPS *, int, char **);
  ~ComputePressureBocs() override;
  void init() override;
  double compute_scalar() override;
  void compute_vector() override;
  void reset_extra_compute_fix(const char *) override;

  double compute_cg_scalar();
  double get_cg_p_corr(int, double *, int, double, double);
  double get_cg_fluct(double, double);
  void send_cg_info(int, int, double *, int, double);
  void send_cg_info(int, double **, int);
  double get_cg_p_corr(double **, int, double);
  double find_index(double *, double);

 protected:
  double boltz, nktv2p, inv_volume;
  int nvirial, dimension;
  double **vptr;
  double *kspace_virial;
  Compute *temperature;
  char *id_temp;
  double virial[6];
  int keflag, pairflag, bondflag, angleflag, dihedralflag, improperflag;
  int fixflag, kspaceflag;

  // NJD MRD
  int p_basis_type;
  int p_match_flag;
  double vavg;
  int N_mol;
  int N_basis;
  double *phi_coeff;
  double **splines;
  int spline_length;

  void virial_compute(int, int);
};

}    // namespace LAMMPS_NS

#endif
#endif
