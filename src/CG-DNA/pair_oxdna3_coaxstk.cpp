// clang-format off
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

#include "pair_oxdna3_coaxstk.h"
#include "nucleotide_oxdna.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairOxdna3Coaxstk::PairOxdna3Coaxstk(LAMMPS *lmp) : PairOxdna2Coaxstk(lmp)
{

  // sequence-specific coaxial stacking strength
  // A:0 C:1 G:2 T:3, 3'- [i] X [j] -5'

  eta_cxst[0][0] = 1.1217958408368172;
  eta_cxst[1][0] = 1.0712851690057155;
  eta_cxst[2][0] = 1.1161603311902566;
  eta_cxst[3][0] = 1.0052361315065244;

  eta_cxst[0][1] = 1.1217958408368172;
  eta_cxst[1][1] = 0.7892685731520542;
  eta_cxst[2][1] = 1.1022201982984874;
  eta_cxst[3][1] = 0.8658975520778347;

  eta_cxst[0][2] = 1.1217958408368172;
  eta_cxst[1][2] = 0.9896542231533637;
  eta_cxst[2][2] = 1.1088392608169480;
  eta_cxst[3][2] = 1.1217958408368172;

  eta_cxst[0][3] = 0.9300223683636719;
  eta_cxst[1][3] = 0.7694592613578328;
  eta_cxst[2][3] = 1.0007533199170144;
  eta_cxst[3][3] = 0.8593983791552220;

  single_enable = 0;
  writedata = 0;
  trim_flag = 0;
}
