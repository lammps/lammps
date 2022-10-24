/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributing Author: David Nicholson (MIT)
------------------------------------------------------------------------- */

#ifndef LMP_FIX_NH_UEF_H
#define LMP_FIX_NH_UEF_H

#include "fix_nh.h"

namespace LAMMPS_NS {
// forward declaration
namespace UEF_utils {
  class UEFBox;
}

class FixNHUef : public FixNH {
 public:
  FixNHUef(class LAMMPS *, int, char **);
  ~FixNHUef() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void pre_exchange() override;
  int pack_restart_data(double *) override;
  void restart(char *) override;
  void end_of_step() override;
  void initial_integrate(int) override;
  void final_integrate() override;
  void initial_integrate_respa(int, int, int) override;
  void final_integrate_respa(int, int) override;
  void post_run() override;
  void get_rot(double[3][3]);
  void get_ext_flags(bool *);
  void get_box(double[3][3]);

 protected:
  void remap() override;
  int size_restart_global() override;
  void nve_x() override;
  void nve_v() override;
  void rotate_x(double[3][3]);
  void inv_rotate_x(double[3][3]);
  void rotate_v(double[3][3]);
  void inv_rotate_v(double[3][3]);
  void rotate_f(double[3][3]);
  void inv_rotate_f(double[3][3]);
  double strain[2], erate[2];    // strain/strain rate : [e_x, e_y]
                                 // always assume traceless e_z = -e_x-e_y

  int rem;    //this is for the narg kluge

  UEF_utils::UEFBox *uefbox;    // interface for the special simulation box

  double rot[3][3];     // rotation matrix
  bool ext_flags[3];    // flags for external "free surfaces"
  bool nearly_equal(double, double, double);
  //bool rotate_output;      // experimental feature. Too many issues for now
};

}    // namespace LAMMPS_NS

#endif
