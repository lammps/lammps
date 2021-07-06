/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef NTOPO_CLASS
// clang-format off
NTopoStyle(NTOPO_DIHEDRAL_TEMPLATE,NTopoDihedralTemplate);
// clang-format on
#else

#ifndef LMP_TOPO_DIHEDRAL_TEMPLATE_H
#define LMP_TOPO_DIHEDRAL_TEMPLATE_H

#include "ntopo.h"

namespace LAMMPS_NS {

class NTopoDihedralTemplate : public NTopo {
 public:
  NTopoDihedralTemplate(class LAMMPS *);
  ~NTopoDihedralTemplate() {}
  void build();
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Dihedral atoms %d %d %d %d missing on proc %d at step %ld

UNDOCUMENTED

W: Dihedral atoms missing at step %ld

UNDOCUMENTED

*/
