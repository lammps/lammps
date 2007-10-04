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

#ifndef ANGLE_HYBRID_H
#define ANGLE_HYBRID_H

#include "stdio.h"
#include "angle.h"

namespace LAMMPS_NS {

class AngleHybrid : public Angle {
 public:
  AngleHybrid(class LAMMPS *);
  ~AngleHybrid();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, int, char **);
  double equilibrium_angle(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  double single(int, int, int, int);
  double memory_usage();

 private:
  int nstyles;                  // # of different angle styles
  Angle **styles;               // class list for each Angle style
  char **keywords;              // keyword for each Angle style
  int *map;                     // which style each angle type points to

  int *nanglelist;              // # of angles in sub-style anglelists
  int *maxangle;                // max # of angles sub-style lists can store
  int ***anglelist;             // anglelist for each sub-style
  
  void allocate();
};

}

#endif
