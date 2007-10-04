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

#ifndef KSPACE_H
#define KSPACE_H

#include "pointers.h"

namespace LAMMPS_NS {

class KSpace : protected Pointers {
 public:
  double energy;
  double virial[6];

  double g_ewald;
  double slab_volfactor;
  int gridflag,gewaldflag;
  int nx_pppm,ny_pppm,nz_pppm;
  int order;
  int slabflag;
 
  KSpace(class LAMMPS *, int, char **);
  virtual ~KSpace() {}
  void modify_params(int, char **);
  virtual void init() = 0;
  virtual void setup() = 0;
  virtual void compute(int, int) = 0;
  virtual void timing(int, double &, double &) {}
  virtual double memory_usage() {return 0.0;}
};

}

#endif
