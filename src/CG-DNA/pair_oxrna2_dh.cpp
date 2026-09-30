/* ----------------------------------------------------------------------
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
   Contributing author: Oliver Henrich (University of Strathclyde, Glasgow)
------------------------------------------------------------------------- */

#include "pair_oxrna2_dh.h"

#include "constants_oxdna.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   compute vector COM-sugar-phosphate backbone interaction site in oxRNA2
------------------------------------------------------------------------- */
void PairOxrna2Dh::compute_backbone_site(double e1[3], double /*e2*/[3],
  double e3[3], double r[3]) const
{
  double dx_cbk_oxdna2 = ConstantsOxdna::get_dx_cbk_oxdna1();
  double dz_cbk_oxrna2 = ConstantsOxdna::get_dz_cbk_oxrna2();

  r[0] = dx_cbk_oxdna2 * e1[0] + dz_cbk_oxrna2 * e3[0];
  r[1] = dx_cbk_oxdna2 * e1[1] + dz_cbk_oxrna2 * e3[1];
  r[2] = dx_cbk_oxdna2 * e1[2] + dz_cbk_oxrna2 * e3[2];
}
