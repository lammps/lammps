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

//
// Contributing author, Richard Meng, Queen's University at Kingston, 10.02.25, contact@richardzjm.com
//

#ifdef PAIR_CLASS
// clang-format off
PairStyle(mtp/extrapolation,PairMTPExtrapolation);
// clang-format on
#else

#ifndef LMP_PAIR_MTP_EXTRAPOLATION_H
#define LMP_PAIR_MTP_EXTRAPOLATION_H

#include "pair_mtp.h"

namespace LAMMPS_NS {

class PairMTPExtrapolation : public PairMTP {
 public:
  PairMTPExtrapolation(class LAMMPS *);
  ~PairMTPExtrapolation() override;
  void compute(int, int) override;                        //Workhorse comuptation
  void settings(int, char **) override;                   // Reads args from "pair_style"
  void *extract(const char *, int &) override;            // Provides access to compute grade flag
  void *extract_peratom(const char *, int &) override;    // Provides access to per-atom data

 protected:
  void read_file(FILE *);                    //Parsing file using LAMMPS utils
  double calculate_extrapolation_grade();    // Grades from candidate vector
  void compile_grades();                     // Collect grades across collective
  void evaluate_grades();                    // Evaluate grades against the thresholds
  void write_config();                       // Write to a MLIP-3 preselected compatible file.

  int coeff_count;    // Sum of radial, species and linear coeff count

  int extrapolation_flag;      // Whether to use extrapolation this iteration (MUST BE INT)
  bool mlip3_style = false;    // Whether to write configs with MLIP-3 compatability

  int configuration_mode;     // Is configuration mode?
  double select_threshold;    // Grade threshold for selection
  double break_threshold;     // Grade threshold for termination
  double max_grade;           // Grade of current iteration

  // Active set
  double **active_set;            // Current active set
  double **inverse_active_set;    // Inverse of the current active set

  //Working buffers
  double ***radial_jacobian;         // Jacobian of radial component wrt to basic moment
  double *radial_moment_ders;        //Ders of non-elemnetary moments wrt to basis moments
  double *energy_ders_wrt_coeffs;    // Candidate information vector

  // Only needed for neigbhourhood mode
  int nbh_count = 0;
  double *nbh_extrapolation_grades = nullptr;    // Extrapolation grades of all neighbourhoods

  // Data for compiling configs in a MLIP-3 compatible format
  FILE *preselected_file;                  // Write to preselected file
  fmt::memory_buffer *write_buffer_ptr;    // Write buffer pointer
};

}    // namespace LAMMPS_NS

#endif
#endif
