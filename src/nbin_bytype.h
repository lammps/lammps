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

#ifdef NBIN_CLASS

NBinStyle(bytype,
          NBinBytype,
          NB_BYTYPE)

#else

#ifndef LMP_NBIN_BYTYPE_H
#define LMP_NBIN_BYTYPE_H

#include "nbin.h"

namespace LAMMPS_NS {

class NBinBytype : public NBin {
 public:

  NBinBytype(class LAMMPS *);
  ~NBinBytype();
  void bin_atoms_setup(int);
  void setup_bins(int);
  void bin_atoms();  

  int coord2bin(double *x, int itype);
  bigint memory_usage();

 private:
  int maxtypes;
  int * maxbins_type;

  void setup_types();
  int itype_min();
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
