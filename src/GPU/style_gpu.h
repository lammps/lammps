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
#endif

#ifdef AtomClass
# endif

#ifdef ComputeInclude
#endif

#ifdef ComputeClass
#endif

#ifdef FixInclude
#endif

#ifdef FixClass
#endif

#ifdef PairInclude
#include "pair_lj_cut_gpu.h"
#include "pair_gayberne_gpu.h"
#endif

#ifdef PairClass
PairStyle(lj/cut/gpu,PairLJCutGPU)
PairStyle(gayberne/gpu,PairGayBerneGPU)
#endif

