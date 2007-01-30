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

#ifndef REGION_INTERSECT_H
#define REGION_INTERSECT_H

#include "region.h"

namespace LAMMPS_NS {

class RegIntersect : public Region {
 public:
  RegIntersect(class LAMMPS *, int, char **);
  ~RegIntersect();
  int match(double, double, double);

 private:
  int nregion;
  int *list;
};

}

#endif
