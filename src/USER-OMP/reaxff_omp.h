/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Hasan Metin Aktulga, Purdue University
   (now at Lawrence Berkeley National Laboratory, hmaktulga@lbl.gov)

   Heavily modified and adapted for LAMMPS by the LAMMPS developers.
------------------------------------------------------------------------- */

#ifndef LMP_REAXFF_OMP_H
#define LMP_REAXFF_OMP_H

#include "reaxff_types.h"

namespace ReaxFF 
{
  // uncomment to enabled collecting ReaxFF OpenMP timing data
  // #define OMP_TIMING 1

#ifdef OMP_TIMING
  // pkcoff timing fields
  enum { COMPUTEINDEX=0,
    COMPUTEWLINDEX,
    COMPUTEBFINDEX,
    COMPUTEQEQINDEX,
    COMPUTENBFINDEX,
    COMPUTEIFINDEX,
    COMPUTETFINDEX,
    COMPUTEBOINDEX,
    COMPUTEBONDSINDEX,
    COMPUTEATOMENERGYINDEX,
    COMPUTEVALENCEANGLESBOINDEX,
    COMPUTETORSIONANGLESBOINDEX,
    COMPUTEHBONDSINDEX,
    COMPUTECG1INDEX,
    COMPUTECG2INDEX,
    COMPUTECGCOMPUTEINDEX,
    COMPUTECALCQINDEX,
    COMPUTEINITMVINDEX,
    COMPUTEMVCOMPINDEX,
    LASTTIMINGINDEX
  };

  extern double ompTimingData[LASTTIMINGINDEX];
  extern int ompTimingCount[LASTTIMINGINDEX];
  extern int ompTimingCGCount[LASTTIMINGINDEX];
#endif

  // exported Functions

  // forces OpenMP

  extern void Compute_ForcesOMP(reax_system *, control_params *,
                                simulation_data *, storage *,
                                reax_list **, output_controls *);

  // init md OpenMP

  extern void InitializeOMP(reax_system *, control_params *,
                            simulation_data *, storage *, reax_list **,
                            output_controls *, MPI_Comm);
}

#endif
