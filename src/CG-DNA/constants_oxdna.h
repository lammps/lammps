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

#ifndef CONSTANTS_OXDNA_H
#define CONSTANTS_OXDNA_H

#include "pointers.h"

namespace LAMMPS_NS {

class ConstantsOxdna : protected Pointers {
 public:
  ConstantsOxdna(class LAMMPS *lmp);
  ~ConstantsOxdna() override = default;

  // oxDNA1 getters
  static double get_dx_cbk_oxdna1() { return dx_cbk_oxdna1; }
  static double get_dx_cstk_oxdna1() { return dx_cstk_oxdna1; }
  static double get_dx_cbs_oxdna1() { return dx_cbs_oxdna1; }

  // oxDNA2 getters
  static double get_dx_cbk_oxdna2() { return dx_cbk_oxdna2; }
  static double get_dy_cbk_oxdna2() { return dy_cbk_oxdna2; }

  // oxDNA3 getters
  static double get_dx_cstk_oxdna3() { return dx_cstk_oxdna3; }
  static double get_dx_cbs_pur_oxdna3() { return dx_cbs_pur_oxdna3; }
  static double get_dx_cbs_pyr_oxdna3() { return dx_cbs_pyr_oxdna3; }

  // oxRNA2 getters
  static double get_dx_cbk_oxrna2() { return dx_cbk_oxrna2; }
  static double get_dz_cbk_oxrna2() { return dz_cbk_oxrna2; }
  static double get_dx_cstk_3p_oxrna2() { return dx_cstk_3p_oxrna2; }
  static double get_dy_cstk_3p_oxrna2() { return dy_cstk_3p_oxrna2; }
  static double get_dx_cstk_5p_oxrna2() { return dx_cstk_5p_oxrna2; }
  static double get_dy_cstk_5p_oxrna2() { return dy_cstk_5p_oxrna2; }

  // electrostatic getters
  static double get_lambda_dh_one_prefactor() { return lambda_dh_one_prefactor; }
  static double get_qeff_dh_pf_one_prefactor() { return qeff_dh_pf_one_prefactor; }

 private:
  std::string units;
  bool real_flag;
  void set_real_units();

  // oxDNA1 parameters
  static double dx_cbk_oxdna1, dx_cstk_oxdna1, dx_cbs_oxdna1;

  // oxDNA2 parameters
  static double dx_cbk_oxdna2, dy_cbk_oxdna2;

  // oxDNA3 parameters
  static double dx_cbk_oxdna3, dy_cbk_oxdna3;
  static double dx_cstk_oxdna3;
  static double dx_cbs_pur_oxdna3, dx_cbs_pyr_oxdna3;

  // oxRNA2 parameters
  static double dx_cbk_oxrna2, dz_cbk_oxrna2;
  static double dx_cstk_3p_oxrna2, dy_cstk_3p_oxrna2;
  static double dx_cstk_5p_oxrna2, dy_cstk_5p_oxrna2;

  // electrostatic parameters
  static double lambda_dh_one_prefactor, qeff_dh_pf_one_prefactor;
};

}    // namespace LAMMPS_NS
#endif
