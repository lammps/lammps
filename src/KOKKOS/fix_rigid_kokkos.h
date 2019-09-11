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

namespace LAMMPS_NS {

template <class DeviceType>
class FixRigidKokkos;

template <class DeviceType>
class FixRigidKokkos : public FixRigid {
 public:
  FixRigidKokkos(class LAMMPS *, int, char **);
  virtual ~FixRigidKokkos();
  void cleanup_copy();

  // virtual int setmask();
  // virtual void init();
  // virtual void setup(int);
  virtual void initial_integrate(int);
  //void post_force(int);
  //virtual void final_integrate();
  // void initial_integrate_respa(int, int, int);
  // void final_integrate_respa(int, int);
  // void write_restart_file(char *);
  // virtual double compute_scalar();



  KOKKOS_INLINE_FUNCTION
  void initial_integrate_item(int) const;
  KOKKOS_INLINE_FUNCTION
  void initial_integrate_rmass_item(int) const;


 private:

  typename ArrayTypes<DeviceType>::t_x_array x;
  typename ArrayTypes<DeviceType>::t_v_array v;
  typename ArrayTypes<DeviceType>::t_f_array f;
  typename ArrayTypes<DeviceType>::t_float_1d rmass;
  typename ArrayTypes<DeviceType>::t_float_1d mass;
  typename ArrayTypes<DeviceType>::t_int_1d type;
  typename ArrayTypes<DeviceType>::t_int_1d mask;
	
	
}; // class FixRigidKokkos

	
template <class DeviceType, int RMass>
struct FixRigidKokkosInitialIntegrateFunctor  {
  typedef DeviceType  device_type ;
  FixRigidKokkos<DeviceType> c;

  FixRigidKokkosInitialIntegrateFunctor(FixRigidKokkos<DeviceType>* c_ptr):
    c(*c_ptr) {c.cleanup_copy();};

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    if (RMass) c.initial_integrate_rmass_item(i);
    else c.initial_integrate_item(i);
  }
};

	/*
template <class DeviceType, int RMass>
struct FixRigidKokkosFinalIntegrateFunctor  {
  typedef DeviceType  device_type ;
  FixNVEKokkos<DeviceType> c;

  FixNVEKokkosFinalIntegrateFunctor(FixNVEKokkos<DeviceType>* c_ptr):
  c(*c_ptr) {c.cleanup_copy();};
  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    if (RMass) c.final_integrate_rmass_item(i);
    else c.final_integrate_item(i);
  }
};
	*/


	
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
