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

#ifdef CommandInclude
#include "prd.h"
#endif

#ifdef CommandClass
CommandStyle(prd,PRD)
#endif

#ifdef ComputeInclude
#include "compute_event_displace.h"
#endif

#ifdef ComputeClass
ComputeStyle(event/displace,ComputeEventDisplace)
#endif

#ifdef FixInclude
#include "fix_event.h"
#endif

#ifdef FixClass
FixStyle(EVENT,FixEvent)
#endif
