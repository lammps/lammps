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

#ifdef FIX_CLASS
// clang-format off
FixStyle(ave/moments,FixAveMoments);
// clang-format on
#else

#ifndef LMP_FIX_AVE_MOMENTS_H
#define LMP_FIX_AVE_MOMENTS_H

#include "fix.h"

namespace LAMMPS_NS {

class FixAveMoments : public Fix {
 public:
  FixAveMoments(class LAMMPS *, int, char **);
  ~FixAveMoments() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void end_of_step() override;
  double compute_scalar() override;
  double compute_vector(int) override;
  double compute_array(int, int) override;

 private:
  struct value_t {
    int which;       // type of data: COMPUTE, FIX, VARIABLE
    int argindex;    // 1-based index if data is vector, else 0
    int iarg;        // argument index in original argument list
    int varlen;      // 1 if value is from variable-length compute
    std::string id;         // compute/fix/variable ID
    std::string keyword;    // column keyword in output
    union {
      class Compute *c;
      class Fix *f;
      int v;
    } val;
  };
  std::vector<value_t> values;
  std::vector<int> moments;

  int nrepeat, nfreq;
  int nvalues;
  bigint nvalid, nvalid_comp_next;

  int startstep;

  int nhistory, iresult;
  double **result_list;

  int iwindow, window_filled;
  double **window_list;

  int consume_moments(int iarg, int narg, char **arg);
  void options(int, int, char **);
  void setnextvalid();

  void get_values(std::vector<double>& scalars);
  void append_values();
  void update_results();
};
}    // namespace LAMMPS_NS
#endif
#endif
