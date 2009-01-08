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
#include "atom_vec_granular.h"
#endif

#ifdef AtomClass
AtomStyle(granular,AtomVecGranular)
# endif

#ifdef FixInclude
#include "fix_freeze.h"
#include "fix_pour.h"
#include "fix_shear_history.h"
#include "fix_wall_gran.h"
#endif

#ifdef FixClass
FixStyle(freeze,FixFreeze)
FixStyle(pour,FixPour)
FixStyle(SHEAR_HISTORY,FixShearHistory)
FixStyle(wall/gran,FixWallGran)
#endif

#ifdef PairInclude
#include "pair_gran_hertz_history.h"
#include "pair_gran_hooke.h"
#include "pair_gran_hooke_history.h"
#endif

#ifdef PairClass
PairStyle(gran/hertz/history,PairGranHertzHistory)
PairStyle(gran/hooke,PairGranHooke)
PairStyle(gran/hooke/history,PairGranHookeHistory)
#endif
