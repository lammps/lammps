/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(nvt/sllod,FixNVTSllod)

#else

#ifndef LMP_FIX_NVT_SLLOD_H
#define LMP_FIX_NVT_SLLOD_H

#include "fix_nh.h"

namespace LAMMPS_NS {

class FixNVTSllod : public FixNH {
 public:
  FixNVTSllod(class LAMMPS *, int, char **);
  ~FixNVTSllod() {}
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
