/* ----------------------------------------------------------------------
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

FixStyle(store/state,FixStoreState)

#else

#ifndef LMP_FIX_STORE_STATE_H
#define LMP_FIX_STORE_STATE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixStoreState : public Fix {
 public:
  FixStoreState(class LAMMPS *, int, char **);
  ~FixStoreState();
  int setmask();
  void init();
  void setup(int);
  void end_of_step();

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  int pack_restart(int, double *);
  void unpack_restart(int, int);
  int size_restart(int);
  int maxsize_restart();

 private:
  int nvalues;
  int *which,*argindex,*value2index;
  char **ids;
  double **values;         // archived atom properties
  double *vbuf;            // 1d ptr to values

  int comflag;
  double cm[3];            // center of mass

  int kflag,cfv_flag,firstflag;
  int cfv_any;             // 1 if any compute/fix/variable specified

  typedef void (FixStoreState::*FnPtrPack)(int);
  FnPtrPack *pack_choice;              // ptrs to pack functions

  void pack_id(int);
  void pack_molecule(int);
  void pack_type(int);
  void pack_mass(int);

  void pack_x(int);
  void pack_y(int);
  void pack_z(int);
  void pack_xs(int);
  void pack_ys(int);
  void pack_zs(int);
  void pack_xs_triclinic(int);
  void pack_ys_triclinic(int);
  void pack_zs_triclinic(int);
  void pack_xu(int);
  void pack_yu(int);
  void pack_zu(int);
  void pack_xu_triclinic(int);
  void pack_yu_triclinic(int);
  void pack_zu_triclinic(int);
  void pack_ix(int);
  void pack_iy(int);
  void pack_iz(int);

  void pack_vx(int);
  void pack_vy(int);
  void pack_vz(int);
  void pack_fx(int);
  void pack_fy(int);
  void pack_fz(int);
  void pack_q(int);
  void pack_mux(int);
  void pack_muy(int);
  void pack_muz(int);
  void pack_radius(int);
  void pack_omegax(int);
  void pack_omegay(int);
  void pack_omegaz(int);
  void pack_angmomx(int);
  void pack_angmomy(int);
  void pack_angmomz(int);
  void pack_tqx(int);
  void pack_tqy(int);
  void pack_tqz(int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix store/state for atom property that isn't allocated

Self-explanatory.

E: Compute ID for fix store/state does not exist

Self-explanatory.

E: Fix store/state compute does not calculate per-atom values

Computes that calculate global or local quantities cannot be used
with fix store/state.

E: Fix store/state compute does not calculate a per-atom vector

The compute calculates a per-atom vector.

E: Fix store/state compute does not calculate a per-atom array

The compute calculates a per-atom vector.

E: Fix store/state compute array is accessed out-of-range

Self-explanatory.

E: Fix ID for fix store/state does not exist

Self-explanatory

E: Fix store/state fix does not calculate per-atom values

Fixes that calculate global or local quantities cannot be used with
fix store/state.

E: Fix store/state fix does not calculate a per-atom vector

The fix calculates a per-atom array.

E: Fix store/state fix does not calculate a per-atom array

The fix calculates a per-atom vector.

E: Fix store/state fix array is accessed out-of-range

Self-explanatory.

E: Fix for fix store/state not computed at compatible time

Fixes generate their values on specific timesteps.  Fix store/state
is requesting a value on a non-allowed timestep.

E: Variable name for fix store/state does not exist

Self-explanatory.

E: Fix store/state variable is not atom-style variable

Only atom-style variables calculate per-atom quantities.

*/
