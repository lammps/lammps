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
FixStyle(ave/time,FixAveTime);
// clang-format on
#else

#ifndef LMP_FIX_AVE_TIME_H
#define LMP_FIX_AVE_TIME_H

#include "fix.h"

#include <map>

namespace LAMMPS_NS {

class FixAveTime : public Fix {
 public:
  FixAveTime(class LAMMPS *, int, char **);
  ~FixAveTime() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void end_of_step() override;
  int modify_param(int, char **) override;
  double compute_scalar() override;
  double compute_vector(int) override;
  double compute_array(int, int) override;

 private:
  struct value_t {
    int which;       // type of data: COMPUTE, FIX, VARIABLE
    int argindex;    // 1-based index if data is vector, else 0
    int varlen;      // 1 if value is from variable-length compute
    int offcol;
    std::string id;         // compute/fix/variable ID
    std::string keyword;    // column keyword in output
    union {
      class Compute *c;
      class Fix *f;
      int v;
    } val;
  };
  std::vector<value_t> values;

  int nvalues, nrepeat, nfreq, irepeat;
  bigint nvalid, nvalid_last;

  FILE *fp;
  int nrows;
  int any_variable_length;
  int all_variable_length;
  int lockforever;
  bool yaml_flag, yaml_header;

  int ave, nwindow, startstep, mode;
  int noff, overwrite;
  int *offlist;
  char *format, *format_user;
  char *title1, *title2, *title3;
  bigint filepos;

  std::map<std::string, int> key2col;

  int norm, iwindow, window_limit;
  double *vector;
  double *vector_total;
  double **vector_list;
  double *column;
  double **array;
  double **array_total;
  double ***array_list;

  int column_length(int);
  void invoke_scalar(bigint);
  void invoke_vector(bigint);
  void options(int, int, char **);
  void allocate_arrays();
  bigint nextvalid();
};
}    // namespace LAMMPS_NS
#endif
#endif
