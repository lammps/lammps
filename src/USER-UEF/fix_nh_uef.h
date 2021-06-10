/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

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
  virtual ~FixNHUef();
  virtual int setmask();
  virtual void init();
  virtual void setup(int);
  virtual void pre_exchange();
  virtual int pack_restart_data(double *);
  virtual void restart(char *);
  virtual void end_of_step();
  virtual void initial_integrate(int);
  virtual void final_integrate();
  virtual void initial_integrate_respa(int, int, int);
  virtual void final_integrate_respa(int, int);
  virtual void post_run();
  void get_rot(double[3][3]);
  void get_ext_flags(bool *);
  void get_box(double[3][3]);

 protected:
  virtual void remap();
  virtual int size_restart_global();
  virtual void nve_x();
  virtual void nve_v();
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

/* ERROR/WARNING messages:

This is a base class for FixNH so it will inherit most of its error/warning messages along with the following:

E: Illegal fix nvt/npt/uef command

Self-explanatory

E: Keyword erate must be set for fix nvt/npt/uef command

Self-explanatory.

E: Simulation box must be triclinic for fix/nvt/npt/uef

Self-explanatory.

E: Only normal stresses can be controlled with fix/nvt/npt/uef

The keywords xy xz and yz cannot be used for pressure control

E: The ext keyword may only be used with iso pressure control

Self-explanatory

E: All controlled stresses must have the same value in fix/nvt/npt/uef

Stress control is only possible when the stress specified for each dimension is the same

E: Dimensions with controlled stresses must have same strain rate in fix/nvt/npt/uef

Stress-controlled dimensions with the same strain rate must have the same target stress.

E: Can't use another fix which changes box shape with fix/nvt/npt/uef

The fix npt/nvt/uef command must have full control over the box shape. You cannot use a simultaneous fix deform command, for example.

E: Pressure ID for fix/nvt/uef doesn't exist

The compute pressure introduced via fix_modify does not exist

E: Using fix nvt/npt/uef without a compute pressure/uef

Self-explanatory.

E: Using fix nvt/npt/uef without a compute temp/uef

Self-explanatory.

E: Initial box is not close enough to the expected uef box

The initial box does not correspond to the shape required by the value of the strain keyword. If the default strain value of zero was used, the initial box is not cubic.

*/
