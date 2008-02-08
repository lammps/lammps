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

// add new include files in appropriate Include ifdef
// add new style keywords and class names in appropriate Class ifdef
// see style.h for examples

#ifdef AngleInclude
#include "angle_cg_cmm.h"
#endif

#ifdef AngleClass
AngleStyle(cg/cmm,AngleCGCMM)
#endif

#ifdef PairInclude
#include "pair_cg_cmm.h"
#include "pair_cg_cmm_coul_cut.h"
#include "pair_cg_cmm_coul_long.h"
#endif

#ifdef PairClass
PairStyle(cg/cmm,PairCGCMM)
PairStyle(cg/cmm/coul/cut,PairCGCMMCoulCut)
PairStyle(cg/cmm/coul/long,PairCGCMMCoulLong)
#endif
