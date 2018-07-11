/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(scafacos,FixScafacos)

#else

#ifndef LMP_FIX_SCAFACOS_H
#define LMP_FIX_SCAFACOS_H

#include "fix.h"
#include "fcs.h"
#include <string>

namespace LAMMPS_NS {

class FixScafacos : public Fix {
 public:
  FixScafacos(class LAMMPS *, int, char **);
  virtual ~FixScafacos();
  int setmask();
  void init();
  void init_list(int, NeighList*);
  void setup(int);
  void min_setup(int);
  void setup_pre_reverse(int, int);
  void initial_integrate(int);
  void pre_reverse(int, int);
  void post_force(int);
  void min_post_force(int);
  void final_integrate();
  void reset_dt();
  double compute_scalar();
  double memory_usage();

 protected:
  std::string method;

  // MPI rank
  int rank;

  // source arrays for positions and charges
  double *x, *q;
  // result arrays for potentials and field
  double *pot, *field;

  // box vectors for each dimension
  fcs_float box_x[3], box_y[3], box_z[3];
  // offset of the box from the origin
  fcs_float offset[3];

  // periodicity of the system
  fcs_int periodicity[3];

  // ScaFaCoS handle
  FCS fcs;

  // ScaFaCoS result variable
  FCSResult result;
  
  // function to check results
  bool check_result(FCSResult, int);

  // function to set up handle with common parameters
  void setup_handle();
  
  // function to check if the box parameters changed, so that a new tuning step is required
  bool box_has_changed();

  // store total number of particles (to check if tune needs to be called again)
  fcs_int total_particles;

  // store number of local particles (to be able to readjust the size of result arrays, when needed)
  int local_array_size;

  // should the near field calculations be computed by LAMMPS?
  fcs_int near_field_flag;

  // type of accuracy chosen (if customized)
  fcs_int tolerance_type;

  // value of tolerance
  fcs_float tolerance_value;

  // is tolerance set?
  bool tolerance_set;

  // check if fmm is chosen (ghost particles, since the implementation needs at least 1 particle on each process!)
  bool fmm_chosen;

  // FMM: fmm particle array size
  int fmm_array_size;
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
