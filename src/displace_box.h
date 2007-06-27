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

#ifndef DISPLACE_BOX_H
#define DISPLACE_BOX_H

#include "pointers.h"

namespace LAMMPS_NS {

class DisplaceBox : protected Pointers {
 public:
  DisplaceBox(class LAMMPS *);
  void command(int, char **);

 private:
  int remapflag,scaleflag;

  struct Set {
    int style,substyle;
    double flo,fhi,ftilt;
    double dlo,dhi,dtilt;
    double scale;
    double lo_start,hi_start;
    double lo_stop,hi_stop;
    double tilt_start,tilt_stop;
    double vol_start;
    int fixed,dynamic1,dynamic2;
  };
  Set *set;

  void options(int, char **);
};

}

#endif
