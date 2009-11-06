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

#ifdef AtomInclude
#include "atom_vec_colloid.h"
#endif

#ifdef AtomClass
AtomStyle(colloid,AtomVecColloid)
#endif

#ifdef FixInclude
#include "fix_wall_colloid.h"
#endif

#ifdef FixClass
FixStyle(wall/colloid,FixWallColloid)
#endif

#ifdef PairInclude
#include "pair_colloid.h"
#include "pair_lubricate.h"
#include "pair_yukawa_colloid.h"
#endif

#ifdef PairClass
PairStyle(colloid,PairColloid)
PairStyle(lubricate,PairLubricate)
PairStyle(yukawa/colloid,PairYukawaColloid)
#endif
