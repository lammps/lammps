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

#ifndef LMP_LABEL_MAP_H
#define LMP_LABEL_MAP_H

#include "pointers.h"  // IWYU pragma: export

namespace LAMMPS_NS {

class LabelMap : protected Pointers {
 public:
   int natomtypes,nbondtypes,nangletypes;
   int ndihedraltypes,nimpropertypes;
   char **typelabel,**btypelabel,**atypelabel;
   char **dtypelabel,**itypelabel;

   LabelMap(LAMMPS *lmp);
   ~LabelMap();

   void allocate_type_labels();
   int find_type(char *, char **, int);

 protected:
};

}

#endif

/* ERROR/WARNING messages:


*/
