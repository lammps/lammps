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

#ifndef COMPUTE_TEMP_COM_H
#define COMPUTE_TEMP_COM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeTempCoM : public Compute {
 public:
  ComputeTempCoM(class LAMMPS *, int, char **);
  ~ComputeTempCoM();
  void init();
  double compute_scalar();
  void compute_vector();
  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);

 private:
  int fix_dof;
  double tfactor;
  double dof1;
  double tfactor1;

  void recount();
};

}

#endif
