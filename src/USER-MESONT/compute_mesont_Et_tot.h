/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributing author: Maxim Shugaev (UVA), mvs9t@virginia.edu
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(mesont/Et_tot,ComputeMESONT_Et_tot)

#else

#ifndef LMP_COMPUTE_MESONT_ET_H
#define LMP_COMPUTE_MESONT_ET_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeMESONT_Et_tot : public Compute {
 public:
  ComputeMESONT_Et_tot(class LAMMPS *, int, char **);
  ~ComputeMESONT_Et_tot() {}
  void init() {}
  double compute_scalar();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal compute mesont/Et_tot command

Incorrect argument list in the compute init.

E: Compute mesont/Et_tot must use group all

UNSPECIFIED.

E: mesont/Et_tot is allowed only with mesont/tpm pair style

Use mesont/tpm pair style.

*/
