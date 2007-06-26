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
#include "pair_eam_opt.h"
#include "pair_eam_alloy_opt.h"
#include "pair_eam_fs_opt.h"
#include "pair_lj_charmm_coul_long_opt.h"
#include "pair_lj_cut_opt.h"
#include "pair_morse_opt.h"
#endif

#ifdef PairClass
PairStyle(eam/opt,PairEAMOpt)
PairStyle(eam/alloy/opt,PairEAMAlloyOpt)
PairStyle(eam/fs/opt,PairEAMFSOpt)
PairStyle(lj/cut/opt,PairLJCutOpt)
PairStyle(lj/charmm/coul/long/opt,PairLJCharmmCoulLongOpt)
PairStyle(morse/opt,PairMorseOpt)
#endif
