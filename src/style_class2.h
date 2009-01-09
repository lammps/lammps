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

#ifdef AngleInclude
#include "angle_class2.h"
#endif

#ifdef AngleClass
AngleStyle(class2,AngleClass2)
#endif

#ifdef BondInclude
#include "bond_class2.h"
#endif

#ifdef BondClass
BondStyle(class2,BondClass2)
#endif

#ifdef DihedralInclude
#include "dihedral_class2.h"
#endif

#ifdef DihedralClass
DihedralStyle(class2,DihedralClass2)
#endif

#ifdef ImproperInclude
#include "improper_class2.h"
#endif

#ifdef ImproperClass
ImproperStyle(class2,ImproperClass2)
#endif

#ifdef PairInclude
#include "pair_lj_class2.h"
#include "pair_lj_class2_coul_cut.h"
#include "pair_lj_class2_coul_long.h"
#endif

#ifdef PairClass
PairStyle(lj/class2,PairLJClass2)
PairStyle(lj/class2/coul/cut,PairLJClass2CoulCut)
PairStyle(lj/class2/coul/long,PairLJClass2CoulLong)
#endif
