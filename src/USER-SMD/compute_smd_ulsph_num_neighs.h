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

ComputeStyle(smd/ulsph/num/neighs,ComputeSMDULSPHNumNeighs)

#else

#ifndef LMP_COMPUTE_SMD_ULSPH_NUM_NEIGHS_H
#define LMP_COMPUTE_SMD_ULSPH_NUM_NEIGHS_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeSMDULSPHNumNeighs : public Compute {
 public:
  ComputeSMDULSPHNumNeighs(class LAMMPS *, int, char **);
  ~ComputeSMDULSPHNumNeighs();
  void init();
  void compute_peratom();
  double memory_usage();

 private:
  int nmax;
  double *numNeighsOutput;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

W: More than one compute ke/atom

It is not efficient to use compute ke/atom more than once.

*/
