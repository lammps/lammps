/* -*- c++ -*- ----------------------------------------------------------
 *
 *                    *** Smooth Mach Dynamics ***
 *
 * This file is part of the USER-SMD package for LAMMPS.
 * Copyright (2014) Georg C. Ganzenmueller, georg.ganzenmueller@emi.fhg.de
 * Fraunhofer Ernst-Mach Institute for High-Speed Dynamics, EMI,
 * Eckerstrasse 4, D-79104 Freiburg i.Br, Germany.
 *
 * ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(smd/tri_surface,PairTriSurf)

#else

#ifndef LMP_SMD_TRI_SURFACE_H
#define LMP_SMD_TRI_SURFACE_H

#include "pair.h"
#include <Eigen/Eigen>


namespace LAMMPS_NS {

class PairTriSurf : public Pair {
 public:
  PairTriSurf(class LAMMPS *);
  virtual ~PairTriSurf();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void init_style();
  void init_list(int, class NeighList *);
  virtual double memory_usage();
  void PointTriangleDistance(const Eigen::Vector3d P, const Eigen::Vector3d TRI1, const Eigen::Vector3d TRI2, const Eigen::Vector3d TRI3,
        Eigen::Vector3d &CP, double &dist);
  double clamp(const double a, const double min, const double max);
  void *extract(const char *, int &);

 protected:
  double **bulkmodulus;
  double **kn;

  double *onerad_dynamic,*onerad_frozen;
  double *maxrad_dynamic,*maxrad_frozen;

  double scale;
  double stable_time_increment; // stable time step size

  void allocate();
};

}

#endif
#endif
