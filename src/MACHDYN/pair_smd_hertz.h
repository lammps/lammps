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
PairStyle(smd/hertz,PairHertz);
// clang-format on
#else

#ifndef LMP_SMD_HERTZ_H
#define LMP_SMD_HERTZ_H

#include "pair.h"

namespace LAMMPS_NS {

class PairHertz : public Pair {
 public:
  PairHertz(class LAMMPS *);
  ~PairHertz() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;
  void init_style() override;
  void init_list(int, class NeighList *) override;
  double memory_usage() override;
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
