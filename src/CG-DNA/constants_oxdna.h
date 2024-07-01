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
  virtual ~ConstantsOxdna(){};

  // oxDNA 1 getters
  static double get_d_cs() { return d_cs; }
  static double get_d_cst() { return d_cst; }
  static double get_d_chb() { return d_chb; }
  static double get_d_cb() { return d_cb; }

  // oxDNA 2 getters
  static double get_d_cs_x() { return d_cs_x; }
  static double get_d_cs_y() { return d_cs_y; }
  static double get_lambda_dh_one_prefactor() { return lambda_dh_one_prefactor; }
  static double get_qeff_dh_pf_one_prefactor() { return qeff_dh_pf_one_prefactor; }

  // oxRNA 2 getters
  static double get_d_cs_z() { return d_cs_z; }
  static double get_d_cst_x_3p() { return d_cst_x_3p; }
  static double get_d_cst_y_3p() { return d_cst_y_3p; }
  static double get_d_cst_x_5p() { return d_cst_x_5p; }
  static double get_d_cst_y_5p() { return d_cst_y_5p; }

 private:
  std::string units;
  bool real_flag;
  void set_real_units();

  // oxDNA 1 parameters
  static double d_cs, d_cst, d_chb, d_cb;

  // oxDNA 2 parameters
  static double d_cs_x, d_cs_y;
  static double lambda_dh_one_prefactor, qeff_dh_pf_one_prefactor;

  // oxRNA 2 parameters
  static double d_cs_z;
  static double d_cst_x_3p, d_cst_y_3p;
  static double d_cst_x_5p, d_cst_y_5p;
};

}    // namespace LAMMPS_NS

#endif
