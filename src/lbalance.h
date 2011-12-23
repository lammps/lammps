/* -*- c++ -*- ----------------------------------------------------------
LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
Transfer Simulations

www.liggghts.com | www.cfdem.com
Christoph Kloss, christoph.kloss@cfdem.com


LIGGGHTS is based on LAMMPS
LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
http://lammps.sandia.gov, Sandia National Laboratories
Steve Plimpton, sjplimp@sandia.gov


Copyright (2003) Sandia Corporation. Under the terms of Contract
DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
certain rights in this software. This software is distributed under
the GNU General Public License.


See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

//NP pure virtual class

#ifndef LMP_LBALANCE_H
#define LMP_LBALANCE_H

#include "neighbor.h"
#include "string.h"

namespace LAMMPS_NS {

class Lbalance : protected Pointers {

 public:
  Lbalance(class LAMMPS *lmp, int narg, char **arg): Pointers(lmp) {
    int n = strlen(arg[0]) + 1;
    style = new char[n];
    strcpy(style,arg[0]);
    iarg = 1;
  }

  virtual ~Lbalance() {
    delete[] style;
  }

  virtual void loadbalance_local_boxes()=0;

  double cutneighmax() {return neighbor->cutneighmax;}

 protected:
  int iarg;
  char *style;
};
}

#endif
