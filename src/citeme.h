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

#ifndef LMP_CITEME_H
#define LMP_CITEME_H

#include "pointers.h"
#include <set>

namespace LAMMPS_NS {

class CiteMe : protected Pointers {
 public:
  CiteMe(class LAMMPS *);
  virtual ~CiteMe();
  void add(const char *);       // print out and register publication

 private:
  FILE *fp;                    // opaque pointer to log.cite file object
  typedef std::set<const char *> citeset;
  citeset *cs;                 // registered set of publications
};

}

#endif

/* ERROR/WARNING messages:

E: Cannot open log.cite file

This file is created when you use some LAMMPS features, to indicate
what paper you should cite on behalf of those who implemented
the feature.  Check that you have write privileges into the directory
you are running in.

*/
