/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Akhlak Mahmood

   Contact:
     Department of Materials Science and Engineering,
     North Carolina State University,
     Raleigh, NC, USA

     amahmoo3@ncsu.edu; mahmoodakhlak@gmail.com
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(mspin, ComputeMspin)

#else

#ifndef LMP_COMPUTE_MSPIN_H
#define LMP_COMPUTE_MSPIN_H

#include "compute.h"

namespace LAMMPS_NS {
class ComputeMspin : public Compute {
 public:
  ComputeMspin(class LAMMPS *, int, char **);
  ~ComputeMspin();
  void init();
  void compute_vector();

 private:
  int irfix;
  char *rfix;

  int ibody, jbody;
};
}    // namespace LAMMPS_NS

#endif
#endif
