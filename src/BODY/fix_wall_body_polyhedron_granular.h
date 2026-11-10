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

#ifdef FIX_CLASS
// clang-format off
FixStyle(wall/body/polyhedron/granular,FixWallBodyPolyhedronGranular);
// clang-format on
#else

#ifndef LMP_FIX_WALL_BODY_POLYHEDRON_GRANULAR_H
#define LMP_FIX_WALL_BODY_POLYHEDRON_GRANULAR_H

#include "fix_wall_body_polyhedron.h"

namespace LAMMPS_NS {

namespace Granular_NS {
  class GranularModel;
}

class FixWallBodyPolyhedronGranular : public FixWallBodyPolyhedron {
 public:
  FixWallBodyPolyhedronGranular(class LAMMPS *, int, char **);
  ~FixWallBodyPolyhedronGranular() override;
  void init() override;

 protected:
  Granular_NS::GranularModel *model;    // contact model of the vertices with the wall

  void contact_geometry(int i, double R);
  double vertex_force(int i, double R, double vnnr, const double *n, const double *vn,
                      const double *vt, double *fw) override;
  void tangential_spring(int i, const double *n, const double *vt, double fne, double *xi,
                         double *ft) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
