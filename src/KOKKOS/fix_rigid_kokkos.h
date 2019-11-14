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

FixStyle(rigid/kk,FixRigidKokkos<LMPDeviceType>)
FixStyle(rigid/kk/device,FixRigidKokkos<LMPDeviceType>)
FixStyle(rigid/kk/host,FixRigidKokkos<LMPHostType>)

#else

#ifndef LMP_FIX_RIGID_KOKKOS_H
#define LMP_FIX_RIGID_KOKKOS_H

#include "fix_rigid.h"
#include "kokkos_type.h"
#include "Kokkos_Random.hpp"

namespace LAMMPS_NS {

template <class DeviceType>
class FixRigidKokkos;

template <class DeviceType>
class FixRigidKokkos : public FixRigid {
 public:

  // These did not exist or at least I could not find them:
  typedef Kokkos::DualView<F_FLOAT*[4], Kokkos::LayoutRight, LMPDeviceType> tdual_quat_array;
  typedef Kokkos::DualView<F_FLOAT*[6], Kokkos::LayoutRight, LMPDeviceType> tdual_sum_array;
  typedef Kokkos::DualView<T_INT*[4], Kokkos::LayoutRight, LMPDeviceType> tdual_int4_array;


  FixRigidKokkos(class LAMMPS *, int, char **);
  virtual ~FixRigidKokkos();

  // virtual int setmask(); // Masks remain same
  virtual void init();
  virtual void setup(int);
  virtual void initial_integrate(int);
  virtual void post_force(int);
  virtual void final_integrate();

  // pre_neighbor gets called explicitly during init. At this time, not all
  // kokkos-able arrays and stuff is set. We have to bypass this somehow.
  // No need for explicit setup_pre_neighbor, it only calls this method
  // which is virtual.
  virtual void pre_neighbor();
  virtual double compute_scalar();


  // void initial_integrate_respa(int, int, int);
  // void final_integrate_respa(int, int);
  // void write_restart_file(char *);
  // virtual double compute_scalar();

  virtual int dof(int);
  //void deform(int);
  //void enforce2d();
  //void reset_dt();
  //void zero_momentum();
  //void zero_rotation();

  void grow_arrays(int);
  void set_xv_kokkos(); // Original set_xv and set_v are also protected.
  void set_v_kokkos();
  void compute_forces_and_torques_kokkos();
  void image_shift_kokkos();
  void apply_langevin_thermostat_kokkos();

  template <int NEIGHFLAG>
  void v_tally(EV_FLOAT &ev, const int &i, double v_arr[6]) const;

  enum SYNC_MODIFY_FLAGS { HOST = 0, DEVICE = 1 };
  template <int space> void sync_all();
  template <int space> void modify_all();

  virtual double extract_ke();
  virtual double extract_erotational();



 private:
  // We need Kokkos style containers for everything in the innner loops:
  DAT::tdual_x_array k_xcm;
  DAT::tdual_v_array k_vcm;
  DAT::tdual_f_array k_fcm;

  DAT::tdual_f_array k_tflag;
  DAT::tdual_f_array k_fflag;

  // Careful. This fix' omega, angmom and torque are defined in fix_rigid.
  // They are not the same as those in atom_vec and atom_vec_kokkos!
  DAT::tdual_v_array k_omega;
  DAT::tdual_v_array k_angmom;
  DAT::tdual_f_array k_torque;
  DAT::tdual_x_array k_inertia;

  // k_quat has to be a special array because it is a quaternion!
  tdual_quat_array k_quat;
  tdual_int4_array k_remapflag;

  DAT::tdual_x_array k_ex_space, k_ey_space, k_ez_space;
  DAT::tdual_x_array k_displace;

  tdual_sum_array k_sum, k_all;
  tdual_sum_array k_langextra;

  DAT::tdual_int_1d k_body, k_eflags;
  DAT::tdual_imageint_1d k_xcmimage;
  DAT::tdual_imageint_1d k_imagebody;
  DAT::tdual_float_1d k_masstotal;
  DAT::tdual_int_1d k_nrigid;


  DAT::tdual_x_array k_orient;
  DAT::tdual_x_array k_dorient;
  DAT::tdual_float_1d k_virial;

  // Needed if we apply langvin forces:
  Kokkos::Random_XorShift64_Pool<DeviceType> rand_pool;
  typedef typename Kokkos::Random_XorShift64_Pool<DeviceType>::generator_type rand_type;


  bool bypass_pre_neighbor;

}; // class FixRigidKokkos


} // namespace LAMMPS_NS




#endif // LMP_FIX_RIGID_KOKKOS_H
#endif // FIX_CLASS

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix rigid custom requires previously defined property/atom

UNDOCUMENTED

E: Fix rigid custom requires integer-valued property/atom

UNDOCUMENTED

E: Variable name for fix rigid custom does not exist

UNDOCUMENTED

E: Fix rigid custom variable is no atom-style variable

UNDOCUMENTED

E: Unsupported fix rigid custom property

UNDOCUMENTED

E: Fix rigid molecule requires atom attribute molecule

Self-explanatory.

E: Too many molecules for fix rigid

The limit is 2^31 = ~2 billion molecules.

E: Could not find fix rigid group ID

A group ID used in the fix rigid command does not exist.

E: One or more atoms belong to multiple rigid bodies

Two or more rigid bodies defined by the fix rigid command cannot
contain the same atom.

E: No rigid bodies defined

The fix specification did not end up defining any rigid bodies.

E: Fix rigid z force cannot be on for 2d simulation

Self-explanatory.

E: Fix rigid xy torque cannot be on for 2d simulation

Self-explanatory.

E: Fix rigid langevin period must be > 0.0

Self-explanatory.

E: Fix rigid npt/nph dilate group ID does not exist

Self-explanatory.

E: One or zero atoms in rigid body

Any rigid body defined by the fix rigid command must contain 2 or more
atoms.

W: More than one fix rigid

It is not efficient to use fix rigid more than once.

E: Rigid fix must come before NPT/NPH fix

NPT/NPH fix must be defined in input script after all rigid fixes,
else the rigid fix contribution to the pressure virial is
incorrect.

W: Cannot count rigid body degrees-of-freedom before bodies are initialized

This means the temperature associated with the rigid bodies may be
incorrect on this timestep.

W: Computing temperature of portions of rigid bodies

The group defined by the temperature compute does not encompass all
the atoms in one or more rigid bodies, so the change in
degrees-of-freedom for the atoms in those partial rigid bodies will
not be accounted for.

E: Fix rigid atom has non-zero image flag in a non-periodic dimension

Image flags for non-periodic dimensions should not be set.

E: Insufficient Jacobi rotations for rigid body

Eigensolve for rigid body was not sufficiently accurate.

E: Fix rigid: Bad principal moments

The principal moments of inertia computed for a rigid body
are not within the required tolerances.

E: Cannot open fix rigid inpfile %s

The specified file cannot be opened.  Check that the path and name are
correct.

E: Unexpected end of fix rigid file

A read operation from the file failed.

E: Fix rigid file has no lines

Self-explanatory.

E: Incorrect rigid body format in fix rigid file

The number of fields per line is not what expected.

E: Invalid rigid body ID in fix rigid file

The ID does not match the number of an existing ID of rigid bodies
that are defined by the fix rigid command.

E: Cannot open fix rigid restart file %s

The specified file cannot be opened.  Check that the path and name are
correct.

*/
