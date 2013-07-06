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

#include "lmptype.h"
#include "pointers.h"

namespace LAMMPS_NS {

class CiteMe : protected Pointers {
 public:
  CiteMe(class LAMMPS *);
  virtual ~CiteMe();

  void on();                    // turn citing on
  void off();                   // turn citing off

  void add(int);         // add a paper to the list of citations

  // constants for references
  enum {
    FIRST_ENTRY=-1,
    PLIMPTON_1995=0,
    PLIMPTON_1997,
    AUHL_2003,
    JANSSENS_2006,
    INTVELD_2008,
    PARKS_2008,
    MUKHERJEE_2008,
    BROWN_2009,
    THOMPSON_2009,
    PETERSEN_2010,
    BROWN_2011,
    BROWN_2012,
    JARAMILLO_BOTERO_2011,
    KONG_2011,
    AKTULGA_2012,
    PLIMPTON_2012,
    SIRK_2013,
    FIORIN_2013,
    LAST_ENTRY
  };

 private:
  void *_pubs;
  bool _active;
};

}

#endif

/* ERROR/WARNING messages:


*/
