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

#ifndef LMP_SURF_EXTRA_H
#define LMP_SURF_EXTRA_H

namespace SurfExtra {

int overlap_sphere_line(double *, double, double *, double *, double *, double *, double &);
int overlap_sphere_tri(double *, double, double *, double *, double *, double *, double *, double *,
                       double &);
int nearest_point_line(double *, double *, double *, double *);
}    // namespace SurfExtra

#endif
