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

#ifdef KSpaceInclude
#include "ewald.h"
#include "pppm.h"
#include "pppm_tip4p.h"
#endif

#ifdef KSpaceClass
KSpaceStyle(ewald,Ewald)
KSpaceStyle(pppm,PPPM)
KSpaceStyle(pppm/tip4p,PPPMTIP4P)
#endif

#ifdef PairInclude
#include "pair_born_coul_long.h"
#include "pair_buck_coul_long.h"
#include "pair_coul_long.h"
#include "pair_lj_cut_coul_long.h"
#include "pair_lj_cut_coul_long_tip4p.h"
#include "pair_lj_charmm_coul_long.h"
#endif

#ifdef PairClass
PairStyle(born/coul/long,PairBornCoulLong)
PairStyle(buck/coul/long,PairBuckCoulLong)
PairStyle(coul/long,PairCoulLong)
PairStyle(lj/cut/coul/long,PairLJCutCoulLong)
PairStyle(lj/cut/coul/long/tip4p,PairLJCutCoulLongTIP4P)
PairStyle(lj/charmm/coul/long,PairLJCharmmCoulLong)
#endif
