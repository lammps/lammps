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
FixStyle(store/state,FixStoreState);
// clang-format on
#else

#ifndef LMP_FIX_STORE_STATE_H
#define LMP_FIX_STORE_STATE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixStoreState : public Fix {
 public:
  FixStoreState(class LAMMPS *, int, char **);
  ~FixStoreState() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void end_of_step() override;

  double memory_usage() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  int pack_restart(int, double *) override;
  void unpack_restart(int, int) override;
  int size_restart(int) override;
  int maxsize_restart() override;

 private:
  struct value_t {
    int which;
    int argindex;
    std::string id;
    union {
      class Compute *c;
      class Fix *f;
      int v;
      int d;
      int i;
    } val;
    void (FixStoreState::* pack_choice)(int);    // ptr to pack function
  };
  std::vector<value_t> values;

  double **avalues;    // archived atom properties
  double *vbuf;       // 1d ptr to values

  int comflag;
  double cm[3];    // center of mass

  int kflag, cfv_flag, firstflag;
  int cfv_any;    // 1 if any compute/fix/variable specified

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
  void pack_xsu(int);
  void pack_ysu(int);
  void pack_zsu(int);
  void pack_xsu_triclinic(int);
  void pack_ysu_triclinic(int);
  void pack_zsu_triclinic(int);

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
  void pack_mu(int);
  void pack_radius(int);
  void pack_diameter(int);
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
}    // namespace LAMMPS_NS
#endif
#endif
