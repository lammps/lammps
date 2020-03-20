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

ComputeStyle(mesont/Es_tot,ComputeMESONT_Es_tot)

#else

#ifndef LMP_COMPUTE_MESONT_ES_H
#define LMP_COMPUTE_MESONT_ES_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeMESONT_Es_tot : public Compute {
 public:
  ComputeMESONT_Es_tot(class LAMMPS *, int, char **);
  ~ComputeMESONT_Es_tot() {}
  void init() {}
  double compute_scalar();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal compute mesont/Es_tot command

Incorrect argument list in the compute init.

E: Compute mesont/Es_tot must use group all

UNSPECIFIED.

E: mesont/Es_tot is allowed only with mesont/tpm pair style

Use mesont/tpm pair style.

*/
