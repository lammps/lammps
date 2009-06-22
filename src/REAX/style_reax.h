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

#ifdef FixInclude
#include "fix_reax_bonds.h"
#endif

#ifdef PairInclude
#include "pair_reax.h"
#endif

#ifdef FixClass
FixStyle(reax/bonds,FixReaxBonds)
#endif

#ifdef PairClass
PairStyle(reax,PairREAX)
#endif
