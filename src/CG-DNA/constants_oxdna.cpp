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
/* ----------------------------------------------------------------------
   Contributing authors: Oliver Henrich (University of Strathclyde, Glasgow)
                         Kierran Falloon (University of Strathclyde, Glasgow)
------------------------------------------------------------------------- */

#include "constants_oxdna.h"

#include "update.h"

namespace LAMMPS_NS {

ConstantsOxdna::ConstantsOxdna(class LAMMPS *lmp) : Pointers(lmp)
{
  // set oxDNA units
  units = update->unit_style;
  real_flag = utils::strmatch(units, "^real");
  if (real_flag) set_real_units();
}

// default to lj units
// oxDNA1 parameters
double ConstantsOxdna::dx_cbk_oxdna1 = -0.4;
double ConstantsOxdna::dx_cstk_oxdna1 = +0.34;
double ConstantsOxdna::dx_cbs_oxdna1 = +0.4;

// oxDNA2 parameters
double ConstantsOxdna::dx_cbk_oxdna2 = -0.34;
double ConstantsOxdna::dy_cbk_oxdna2 = +0.3408;

// oxDNA3 parameters
double ConstantsOxdna::dx_cstk_oxdna3 = +0.37;
double ConstantsOxdna::dx_cbs_pur_oxdna3 = +0.43;
double ConstantsOxdna::dx_cbs_pyr_oxdna3 = +0.37;

// oxRNA2 parameters
double ConstantsOxdna::dx_cbk_oxrna2 = -0.4;
double ConstantsOxdna::dz_cbk_oxrna2 = +0.2;
double ConstantsOxdna::dx_cstk_3p_oxrna2 = +0.4;
double ConstantsOxdna::dy_cstk_3p_oxrna2 = +0.1;
double ConstantsOxdna::dx_cstk_5p_oxrna2 = +0.124906078525;
double ConstantsOxdna::dy_cstk_5p_oxrna2 = -0.00866274917473;

// electrostatic parameters
double ConstantsOxdna::lambda_dh_one_prefactor = +0.3616455075438555;      // = C1
double ConstantsOxdna::qeff_dh_pf_one_prefactor = +0.08173808693529228;    // = C2

void ConstantsOxdna::set_real_units()
{
  // oxDNA1 parameters in real units
  dx_cbk_oxdna1 = -3.4072;
  dx_cstk_oxdna1 = +2.89612;
  dx_cbs_oxdna1 = +3.4072;

  // oxDNA2 parameters in real units
  dx_cbk_oxdna2 = -2.89612;
  dy_cbk_oxdna2 = +2.9029344;

  // oxDNA3 parameters in real units
  dx_cstk_oxdna3 = +3.15166;
  dx_cbs_pur_oxdna3 = +3.66274;
  dx_cbs_pyr_oxdna3 = +3.15166;

  // oxRNA2 parameters in real units
  dx_cbk_oxrna2 = -3.4072;
  dz_cbk_oxrna2 = +1.7036;
  dx_cstk_3p_oxrna2 = +3.4072;
  dy_cstk_3p_oxrna2 = +0.8518;
  dx_cstk_5p_oxrna2 = +1.063949977;
  dy_cstk_5p_oxrna2 = -0.07378929747;

  // electrostatic parameters in real units
  lambda_dh_one_prefactor = +0.05624154892;        // = C1 * 8.518 * sqrt(k_B/4.142e-20)
  qeff_dh_pf_one_prefactor = +4.15079634587587;    // = C2 * 5.961689060210325 * 8.518
};

}    // namespace LAMMPS_NS
