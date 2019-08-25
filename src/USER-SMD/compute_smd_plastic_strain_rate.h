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

#ifdef COMPUTE_CLASS

ComputeStyle(smd/plastic/strain/rate,ComputeSMDPlasticStrainRate)

#else

#ifndef LMP_COMPUTE_SMD_PLASTIC_STRAIN_RATE_H
#define LMP_COMPUTE_SMD_PLASTIC_STRAIN_RATE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeSMDPlasticStrainRate : public Compute {
 public:
  ComputeSMDPlasticStrainRate(class LAMMPS *, int, char **);
  ~ComputeSMDPlasticStrainRate();
  void init();
  void compute_peratom();
  double memory_usage();

 private:
  int nmax;
  double *plastic_strain_rate_vector;
};

}

#endif
#endif
