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

#ifdef PairInclude
#include "pair_airebo.h"
#include "pair_eam.h"
#include "pair_eam_alloy.h"
#include "pair_eam_fs.h"
#include "pair_sw.h"
#include "pair_tersoff.h"
#include "pair_tersoff_zbl.h"
#endif

#ifdef PairClass
PairStyle(airebo,PairAIREBO)
PairStyle(eam,PairEAM)
PairStyle(eam/alloy,PairEAMAlloy)
PairStyle(eam/fs,PairEAMFS)
PairStyle(sw,PairSW)
PairStyle(tersoff,PairTersoff)
PairStyle(tersoff/zbl,PairTersoffZBL)
#endif
