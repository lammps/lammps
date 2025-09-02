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
   Contributing author: David Immel (d.immel@fz-juelich.de, FZJ, Germany)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(atom_weight/apip,FixAtomWeightAPIP);
// clang-format on
#else

#ifndef LMP_FIX_ATOM_WEIGHT_APIP_H
#define LMP_FIX_ATOM_WEIGHT_APIP_H

#include "fix.h"

namespace LAMMPS_NS {

/**
 *  Small wrapper for the lammps timers to get time intervals.
 */

class APIPtimer : protected Pointers {
 public:
  APIPtimer(class LAMMPS *);
  void init();
  double get_work();

 private:
  double time[4];             ///< value of lammps timers
  double time_interval[4];    ///< time interval between two calls

  void set_time();
};

/**
 *  Fix to compute an atomic wegiht that can be used to load-balance an adaptive-precision potential.
 */

class FixAtomWeightAPIP : public Fix {
 public:
  FixAtomWeightAPIP(class LAMMPS *, int, char **);
  ~FixAtomWeightAPIP() override;
  int setmask() override;
  void post_constructor() override;
  void init() override;
  void end_of_step() override;
  void pre_exchange() override;
  void setup_pre_exchange() override;
  void setup_pre_force(int) override;
  double compute_vector(int) override;

 private:
  bigint last_calc;    ///< last timestep in which weights were calculated

  // user set variables

  // clang-format off
  bool rescale_work;          ///< rescale total work of all atoms according to our model to the measured total work
  char * time_simple_extract_name;  ///< argument for pair->extract to get the simple per_atom_time
  char * time_complex_extract_name; ///< argument for pair->extract to get the complex per_atom_time
  char * time_group_extract_name;   ///< argument for pair->extract to get the time independent of lambda
  char * time_lambda_extract_name;  ///< argument for pair->extract to get the time for lambda max calculation
  char * time_group_name;     ///< name of group corresponding to time_group_extract_name
  int time_group_i;          ///< id of group corresponding to time_group_extract_name
  int time_group_bit;         ///< groupbit corresponding to time_group_extract_name
  // clang-format on

  // calculated variables

  double time_simple_atom;     ///< time per particle of simple pair style
  double time_complex_atom;    ///< time per particle of precise pair style
  double time_group_atom;      ///< time per particle independent of lambda
  double time_lambda_atom;     ///< time per atom with non-simple lambda_own
  int n_simple;                ///< number of simple particles before load balancing z
  int n_complex;               ///< number of complex particles before load balancing z
  int n_group;                 ///< number of group particles before load balancing z
  int n_lambda;                ///< numebr of own_complex particles before load balancing z

  double avg_time_atom[4];    ///< global average of per atom times per type

  // class pointers

  class FixStoreAtom *fixstore;    ///< per-atom weights are stored in FixStore
  class APIPtimer *ap_timer;       ///< wrapper for lammps timers
  class Fix *fix_lambda;           ///< ptr to fix_lambda to extract information

  void calc_work_per_particle();
};

}    // namespace LAMMPS_NS

#endif
#endif
