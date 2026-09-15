/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef TEST_ANALYTIC_MODELS_H
#define TEST_ANALYTIC_MODELS_H

#include "test_config.h"

namespace LAMMPS_NS {
class LAMMPS;
}

// Evaluate the closed-form ("analytic") model selected in the YAML config
// after run segment `segment` has completed.  Does nothing unless
// cfg.analytic_enable is set and `segment` matches cfg.analytic_segment
// (a negative analytic_segment means "the last segment").  Parameters are
// taken from the cfg.variables block; the relevant per-atom state is read
// from the live LAMMPS instance and compared with cfg.analytic_tol.
void check_analytic_model(const TestConfig &cfg, LAMMPS_NS::LAMMPS *lmp, int segment);

#endif
