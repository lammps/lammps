/* -*- c++ -*- ----------------------------------------------------------
 *
 *                    *** Smooth Mach Dynamics ***
 *
 * This file is part of the MACHDYN package for LAMMPS.
 * Copyright (2014) Georg C. Ganzenmueller, georg.ganzenmueller@emi.fhg.de
 * Fraunhofer Ernst-Mach Institute for High-Speed Dynamics, EMI,
 * Eckerstrasse 4, D-79104 Freiburg i.Br, Germany.
 *
 * ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(smd/tri_surface,PairTriSurf);
// clang-format on
#else

#ifndef LMP_SMD_TRI_SURFACE_H
#define LMP_SMD_TRI_SURFACE_H

#include "pair.h"
#include <Eigen/Eigen>

namespace LAMMPS_NS {

class PairTriSurf : public Pair {
 public:
  PairTriSurf(class LAMMPS *);
  ~PairTriSurf() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  void init_style() override;
  void init_list(int, class NeighList *) override;
  double memory_usage() override;
  void PointTriangleDistance(const Eigen::Vector3d &P, const Eigen::Vector3d &TRI1,
                             const Eigen::Vector3d &TRI2, const Eigen::Vector3d &TRI3,
                             Eigen::Vector3d &CP, double &dist);
  double clamp(const double a, const double min, const double max);
  void *extract(const char *, int &) override;

 protected:
  double **bulkmodulus;
  double **kn;

  double *onerad_dynamic, *onerad_frozen;
  double *maxrad_dynamic, *maxrad_frozen;

  double scale;
  double stable_time_increment;    // stable time step size

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
