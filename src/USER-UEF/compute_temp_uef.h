/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributing author: David Nicholson (MIT)
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(temp/uef,ComputeTempUef)

#else

#ifndef LMP_COMPUTE_TEMP_UEF_H
#define LMP_COMPUTE_TEMP_UEF_H

#include "compute_temp.h"

namespace LAMMPS_NS {

class ComputeTempUef : public ComputeTemp {
 public:
  ComputeTempUef(class LAMMPS *, int, char **);
  virtual ~ComputeTempUef(){}
  virtual void init();
  virtual void compute_vector();
  void yes_rot();
  void no_rot();


 protected:
  bool rot_flag;
  void virial_rot(double*,const double[3][3]);
  int ifix_uef;
};


}

#endif
#endif

/* ERROR/WARNING messages:

This class inherits most of the warnings from ComputePressure. The
only addition is:

E: Can't use compute temp/uef without defining a fix nvt/npt/uef

Self-explanatory.

*/
