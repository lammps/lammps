/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: W. Michael Brown (Intel)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(nvt/sllod/intel,FixNVTSllodIntel)

#else

#ifndef LMP_FIX_NVTSLLOD_INTEL_H
#define LMP_FIX_NVTSLLOD_INTEL_H

#include "fix_nh_intel.h"

namespace LAMMPS_NS {

class FixNVTSllodIntel : public FixNHIntel {
 public:
  FixNVTSllodIntel(class LAMMPS *, int, char **);
  ~FixNVTSllodIntel() {}
  void init();

 private:
  int nondeformbias;

  void nh_v_temp();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Temperature control must be used with fix nvt/sllod

Self-explanatory.

E: Pressure control can not be used with fix nvt/sllod

Self-explanatory.

E: Temperature for fix nvt/sllod does not have a bias

The specified compute must compute temperature with a bias.

E: Using fix nvt/sllod with inconsistent fix deform remap option

Fix nvt/sllod requires that deforming atoms have a velocity profile
provided by "remap v" as a fix deform option.

E: Using fix nvt/sllod with no fix deform defined

Self-explanatory.

*/

