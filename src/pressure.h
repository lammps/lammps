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

#ifndef PRESSURE_H
#define PRESSURE_H

#include "lammps.h"

class Temperature;

class Pressure : public LAMMPS {
 public:
  double p_total;
  double p_tensor[6],virial[6];

  Pressure() {}
  ~Pressure() {}
  void init();
  void compute(Temperature *);

 private:
  double boltz,nktv2p;
  double *pair_virial,*bond_virial,*angle_virial;
  double *dihedral_virial,*improper_virial,*kspace_virial;
  double *shake_virial,*rigid_virial,*poems_virial;
  int tensorflag;
  int pairflag,bondflag,angleflag,dihedralflag,improperflag,kspaceflag;
  int shakeflag,bodyflag,rigidflag,poemsflag;
};

#endif
