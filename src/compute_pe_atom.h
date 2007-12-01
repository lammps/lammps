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

#ifndef COMPUTE_PE_ATOM_H
#define COMPUTE_PE_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputePEAtom : public Compute {
 public:
  ComputePEAtom(class LAMMPS *, int, char **);
  ~ComputePEAtom();
  void init() {}
  void compute_peratom();
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double memory_usage();

 private:
  int pairflag,bondflag,angleflag,dihedralflag,improperflag;
  int nmax;
  double *energy;
};

}

#endif
