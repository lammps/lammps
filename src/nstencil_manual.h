/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_NSTENCIL_MANUAL_H
#define LMP_NSTENCIL_MANUAL_H

#include "nstencil.h"

namespace LAMMPS_NS {

class NStencilManual : public NStencil {
 public:
  NStencilManual(class LAMMPS *);
  void create() override;

  int half;

  void assign_neighbor_info(double);
 protected:
  int tri, dim_3d;
};

}    // namespace LAMMPS_NS

#endif
