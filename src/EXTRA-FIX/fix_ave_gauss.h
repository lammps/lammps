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
FixStyle(ave/gauss,FixAveGauss)
// clang-format on
#else

#ifndef LMP_FIX_AVE_GAUSS_H
#define LMP_FIX_AVE_GAUSS_H

#include "fix.h"

namespace LAMMPS_NS {

class FixAveGauss : public Fix {
 public:
  FixAveGauss(class LAMMPS *, int, char **);
  ~FixAveGauss();
  int setmask();
  void init();
  void setup(int);
  void end_of_step();
  double compute_scalar();
  double compute_vector(int);
  double compute_array(int,int);

 private:
  struct value_t {
    int which;       // type of data: VARIABLE
    int argindex;    // 1-based index if vector access, else 0
    int iarg;        // argument index in original argument list
    std::string id;  // variable name
    int v;           // variable index
  };
  std::vector<value_t> values;

  int nfreq;
  bigint nvalid, nvalid_last, nfull_next;

  int nvalues;
  std::vector<std::string> value_ids;
  std::vector<int> value_vars;

  int nwindow;
  int iwindow, window_filled;
  double **window_list;

  int nresult, iresult;
  std::vector<int> delays;
  double **result_list;

  void append_values(bigint);
  void update_results(bigint);
  void options(int, int, char **);
};
}    // namespace LAMMPS_NS
#endif
#endif
