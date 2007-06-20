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
#include "atom_vec_ellipsoid.h"
#endif

#ifdef AtomClass
AtomStyle(ellipsoid,AtomVecEllipsoid)
# endif

#ifdef ComputeInclude
#include "compute_temp_asphere.h"
#endif

#ifdef ComputeClass
ComputeStyle(temp/asphere,ComputeTempAsphere)
#endif

#ifdef FixInclude
#include "fix_nve_asphere.h"
#include "fix_nvt_asphere.h"
#include "fix_npt_asphere.h"
#endif

#ifdef FixClass
FixStyle(nve/asphere,FixNVEASphere)
FixStyle(nvt/asphere,FixNVTASphere)
FixStyle(npt/asphere,FixNPTASphere)
#endif

#ifdef PairInclude
#include "pair_gayberne.h"
#endif

#ifdef PairClass
PairStyle(gayberne,PairGayBerne)
#endif
