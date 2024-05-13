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

#ifdef REGION_CLASS
// clang-format off
RegionStyle(prism,RegPrism);
// clang-format on
#else

#ifndef LMP_REGION_PRISM_H
#define LMP_REGION_PRISM_H

#include "region.h"

namespace LAMMPS_NS {

class RegPrism : public Region {
  friend class CreateBox;

 public:
  RegPrism(class LAMMPS *, int, char **);
  ~RegPrism() override;
  int inside(double, double, double) override;
  int surface_interior(double *, double) override;
  int surface_exterior(double *, double) override;

 private:
  double xlo, xhi, ylo, yhi, zlo, zhi;
  double xy, xz, yz;
  double h[3][3], hinv[3][3];
  int dimension;
  double a[3], b[3], c[3];    // edge vectors of region
  double clo[3], chi[3];      // opposite corners of prism
  double face[6][3];          // unit normals of 6 prism faces
  double corners[8][3];       // 8 corner pts of prism
  int tri[12][3];             // 3 corner pts of 12 triangles (2 per face)

  void find_nearest(double *, double &, double &, double &);
  int inside_tri(double *, double *, double *, double *, double *);
  double closest(double *, double *, double *, double);
};

}    // namespace LAMMPS_NS

#endif
#endif
