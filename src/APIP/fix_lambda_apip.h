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
FixStyle(lambda/apip,FixLambdaAPIP);
// clang-format on
#else

#ifndef LMP_FIX_LAMBDA_APIP_H
#define LMP_FIX_LAMBDA_APIP_H

#include "fix.h"

namespace LAMMPS_NS {

class FixLambdaAPIP : public Fix {
  friend class PairLambdaInputAPIP;
  friend class PairLambdaZoneAPIP;

 public:
  FixLambdaAPIP(class LAMMPS *, int, char **);
  ~FixLambdaAPIP() override;
  int modify_param(int, char **) override;
  void init() override;
  int setmask() override;
  void post_constructor() override;
  void setup_pre_force(int) override;
  void post_integrate() override;
  void end_of_step() override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  void write_restart(FILE *) override;
  void restart(char *) override;
  void *extract(const char *, int &) override;

 private:
  enum { FORWARD_MAX, FORWARD_TA };
  int comm_forward_flag;    // flag that determines which variables are communicated in comm forward

  class PairLambdaInputAPIP *pair_lambda_input;
  class PairLambdaZoneAPIP *pair_lambda_zone;

  double cut_lo;    ///< distance at which the cutoff function of the transition zone decays from 1
  double cut_hi;    ///< distance at which the cutoff function of the transition zone is 0
  double cut_width;           ///< cut_hi - cut_lo
  double cut_hi_sq;           ///< cut_hi_sq * cut_hi_sq
  double threshold_lo;        ///< threshold above which the fast potential starts to be turned off
  double threshold_hi;        ///< threshold above which the fast potential is turned off completely
  double threshold_width;     ///< threshold_hi - threshold_lo
  double min_delta_lambda;    ///< minimum steps size of lambda

  double lambda_non_group;    ///< lambda for atoms that are not in the group of this fix

  int history_length;              ///< number of steps of which lambda_input is sterd at max
  int history_used;                ///< number of steps of which lambda_input is currently stored
  int history_last;                ///< index of last written lambda_input value
  class FixStoreAtom *fixstore;    ///< ptr to stored lambda_input values

  int history2_length;              ///< number of steps of which lambda is stored at max
  int history2_used;                ///< number of steps of which lambda is currently stored
  int history2_last;                ///< index of last written lambda value
  class FixStoreAtom *fixstore2;    ///< ptr to stored lambda values

  double switching_function_poly(double);

  bigint invoked_history_update;     ///< last timestep with stored history
  bigint invoked_history2_update;    ///< last timestep with stored history
  bool dump_history_flag;            ///< dump whole stored history (for debugging)

  int group_bit_simple;                    ///< hard coded simple atoms
  int group_bit_complex;                   ///< hard coded complex atoms
  int group_bit_ignore_lambda_input;       ///< ignore lambda_input of this group
  char *group_name_simple;                 ///< hard coded simple atoms
  char *group_name_complex;                ///< hard coded complex atoms
  char *group_name_ignore_lambda_input;    ///< ignore lambda_input of this group

  int nmax_stats;            ///< number of allocated atoms for peratom_stats
  double **peratom_stats;    ///< returned peratom vector

  void calculate_lambda_input_ta();
  void comm_forward_lambda();
  void write_peratom_stats();
  void update_lambda_input_history();
  void communicate_lambda_input_ta();
  void get_lambda_average();
  void update_lambda_history();
};

}    // namespace LAMMPS_NS

#endif
#endif
