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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(airebo/bc,PairAIREBObc);
// clang-format on
#else

#ifndef LMP_PAIR_AIREBO_BC_H
#define LMP_PAIR_AIREBO_BC_H

#include "pair_airebo.h"

namespace LAMMPS_NS {

/* ----------------------------------------------------------------------
   AIREBO-BC: bond-centric modification of the AIREBO bond-order potential
   (Hur & Stuart, J. Chem. Phys. 137, 054102 (2012)).

   The only physical change relative to AIREBO is that the P_CC coordination
   correction is made bond-centric: P_CC becomes a function of the bond-averaged
   coordination numbers Nbar^t = 1/2 (N_ij^t + N_ji^t), stored on a half-integer
   spline grid (Table III of the paper).  All of that logic lives in the
   PairAIREBO base class and is selected by bcflag, so this variant is simply
   the base class with bcflag = 1 (mirroring how airebo/morse sets morseflag).
------------------------------------------------------------------------- */

class PairAIREBObc : public PairAIREBO {
 public:
  PairAIREBObc(class LAMMPS *);
};

}    // namespace LAMMPS_NS

#endif
#endif
