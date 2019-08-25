/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(saed,ComputeSAED)

#else

#ifndef LMP_COMPUTE_SAED_H
#define LMP_COMPUTE_SAED_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeSAED : public Compute {
 public:
  ComputeSAED(class LAMMPS *, int, char **);
  ~ComputeSAED();
  void    init();
  void    compute_vector();
  double  memory_usage();
//testing
  double  saed_var[10];

 private:
  int     me;
  int     *ztype;            // Atomic number of the different atom types
  double  c[3];              // Parameters controlling resolution of reciprocal space explored
  double  dR_Ewald;          // Thickness of Ewald sphere slice
  double  prd_inv[3];        // Inverse spacing of unit cell
  bool    echo;              // echo compute_array progress
  bool    manual;            // Turn on manual recpiprocal map
  int     nRows;             // Number of relp explored

  double  Zone[3];           // Zone axis to view SAED
  double  R_Ewald;           // Radius of Ewald sphere (distance units)
  double  lambda;            // Radiation wavelenght (distance units)
  double  dK[3];             // spacing of reciprocal points in each dimension
  int     Knmax[3];          // maximum integer value for K points in each dimension
  double  Kmax;              // Maximum reciprocal distance to explore

  int ntypes;
  int nlocalgroup;
  int *store_tmp;

};

}

#endif
#endif
