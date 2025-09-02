/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(hmc,FixHMC)
// clang-format on
#else

#ifndef LMP_FIX_HMC_H
#define LMP_FIX_HMC_H

#include "fix.h"

namespace LAMMPS_NS {

// forward declarations
class Compute;
class FixRigidSmall;
class RanPark;

class FixHMC : public Fix {
 public:
  FixHMC(class LAMMPS *, int, char **);
  ~FixHMC() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void end_of_step() override;
  double compute_scalar() override;
  double compute_vector(int) override;
  double memory_usage() override;

 private:
  void setup_arrays_and_pointers();
  void save_current_state();
  void restore_saved_state();
  void random_velocities();
  void rigid_body_random_velocities();
  void grow_store(int, int);

  int nstore, maxstore, bufextra;
  double *buf_store;

  int resample_on_accept_flag, mom_flag, flag_rigid;
  int first_init_complete, first_setup_complete;

  char *id_rigid;
  FixRigidSmall *fix_rigid;

  int nattempts, naccepts;
  double KT, mbeta;
  double PE, KE;
  double DeltaPE, DeltaKE;

  RanPark *random;
  RanPark *random_equal;

  int ne;
  int neg;
  double *eglobal;
  double **eglobalptr;

  int nv;
  double **vglobal;
  double ***vglobalptr;

  Compute *pe;
  Compute *ke;
  Compute *peatom;
  Compute *press;
  Compute *pressatom;

  int peatom_flag;
  int press_flag;
  int pressatom_flag;

  int nvalues;
};

}    // namespace LAMMPS_NS

#endif
#endif
