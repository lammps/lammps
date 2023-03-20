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
/* ----------------------------------------------------------------------
   Contributing author: Ilya Valuev (JIHT RAS)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(awpmd/cut,PairAWPMDCut);
// clang-format on
#else

#ifndef LMP_PAIR_AWPMD_CUT_H
#define LMP_PAIR_AWPMD_CUT_H

#include "pair.h"

class AWPMD_split;

namespace LAMMPS_NS {

class PairAWPMDCut : public Pair {
  friend class FixNVEAwpmd;

 public:
  PairAWPMDCut(class LAMMPS *);
  ~PairAWPMDCut() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  void min_pointers(double **, double **);
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;

  void min_xf_pointers(int, double **, double **) override;
  void min_xf_get(int) override;
  void min_x_set(int) override;
  double memory_usage() override;

 private:
  int flexible_pressure_flag;
  double cut_global;
  double **cut;

  int nmax;                          // number of additional variables for minimizer
  double *min_var, *min_varforce;    // additional variables for minimizer

  void allocate();

  void virial_eradius_compute();

  AWPMD_split *wpmd;         // solver object
  double ermscale;           // scale of width mass for motion
  double width_pbc;          // setting for width pbc
  double half_box_length;    // calculated by coeff function
};

}    // namespace LAMMPS_NS

#endif
#endif
