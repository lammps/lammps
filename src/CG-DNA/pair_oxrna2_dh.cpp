/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

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

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
    compute vector COM-sugar-phosphate backbone interaction site in oxRNA2
------------------------------------------------------------------------- */
void PairOxrna2Dh::compute_interaction_sites(double e1[3], double /*e2*/[3], double e3[3],
                                             double r[3])
{
  double d_cs_x = -0.4, d_cs_z = +0.2;

  r[0] = d_cs_x * e1[0] + d_cs_z * e3[0];
  r[1] = d_cs_x * e1[1] + d_cs_z * e3[1];
  r[2] = d_cs_x * e1[2] + d_cs_z * e3[2];
}
