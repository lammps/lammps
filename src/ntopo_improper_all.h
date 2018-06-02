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

#ifdef NTOPO_CLASS

NTopoStyle(NTOPO_IMPROPER_ALL,NTopoImproperAll)

#else

#ifndef LMP_TOPO_IMPROPER_ALL_H
#define LMP_TOPO_IMPROPER_ALL_H

#include "ntopo.h"

namespace LAMMPS_NS {

class NTopoImproperAll : public NTopo {
 public:
  NTopoImproperAll(class LAMMPS *);
  ~NTopoImproperAll() {}
  void build();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Improper atoms %d %d %d %d missing on proc %d at step %ld

UNDOCUMENTED

W: Improper atoms missing at step %ld

UNDOCUMENTED

*/
