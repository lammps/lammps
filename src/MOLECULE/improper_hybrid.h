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

#ifndef IMPROPER_HYBRID_H
#define IMPROPER_HYBRID_H

#include "stdio.h"
#include "improper.h"

namespace LAMMPS_NS {

class ImproperHybrid : public Improper {
 public:
  ImproperHybrid(class LAMMPS *);
  ~ImproperHybrid();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, int, char **);
  void write_restart(FILE *);
  void read_restart(FILE *);
  double memory_usage();

 private:
  int nstyles;                  // # of different improper styles
  Improper **styles;            // class list for each Improper style
  char **keywords;              // keyword for each improper style
  int *map;                     // which style each improper type points to

  int *nimproperlist;           // # of impropers in sub-style improperlists
  int *maximproper;             // max # of impropers sub-style lists can store
  int ***improperlist;          // improperlist for each sub-style
  
  void allocate();
};

}

#endif
