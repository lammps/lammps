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

#ifndef LMP_ACCELERATOR_OMP_H
#define LMP_ACCELERATOR_OMP_H

// true interface to USER-OMP, used in Neighbor class header file
// used when USER-OMP is installed

#ifdef LMP_USER_OMP

  void half_nsq_no_newton_omp(class NeighList *);
  void half_nsq_newton_omp(class NeighList *);

  void half_bin_no_newton_omp(class NeighList *);
  void half_bin_newton_omp(class NeighList *);
  void half_bin_newton_tri_omp(class NeighList *);

  void half_multi_no_newton_omp(class NeighList *);
  void half_multi_newton_omp(class NeighList *);
  void half_multi_newton_tri_omp(class NeighList *);

  void full_nsq_omp(class NeighList *);
  void full_nsq_ghost_omp(class NeighList *);
  void full_bin_omp(class NeighList *);
  void full_bin_ghost_omp(class NeighList *);
  void full_multi_omp(class NeighList *);

  void half_from_full_no_newton_omp(class NeighList *);
  void half_from_full_newton_omp(class NeighList *);

  void granular_nsq_no_newton_omp(class NeighList *);
  void granular_nsq_newton_omp(class NeighList *);
  void granular_bin_no_newton_omp(class NeighList *);
  void granular_bin_newton_omp(class NeighList *);
  void granular_bin_newton_tri_omp(class NeighList *);

  void respa_nsq_no_newton_omp(class NeighList *);
  void respa_nsq_newton_omp(class NeighList *);
  void respa_bin_no_newton_omp(class NeighList *);
  void respa_bin_newton_omp(class NeighList *);
  void respa_bin_newton_tri_omp(class NeighList *);

#else

// dummy interface to USER-OMP
// needed for compiling Neighbor class when USER-OMP is not installed

  void half_nsq_no_newton_omp(class NeighList *) {}
  void half_nsq_newton_omp(class NeighList *) {}

  void half_bin_no_newton_omp(class NeighList *) {}
  void half_bin_newton_omp(class NeighList *) {}
  void half_bin_newton_tri_omp(class NeighList *) {}

  void half_multi_no_newton_omp(class NeighList *) {}
  void half_multi_newton_omp(class NeighList *) {}
  void half_multi_newton_tri_omp(class NeighList *) {}

  void full_nsq_omp(class NeighList *) {}
  void full_nsq_ghost_omp(class NeighList *) {}
  void full_bin_omp(class NeighList *) {}
  void full_bin_ghost_omp(class NeighList *) {}
  void full_multi_omp(class NeighList *) {}

  void half_from_full_no_newton_omp(class NeighList *) {}
  void half_from_full_newton_omp(class NeighList *) {}

  void granular_nsq_no_newton_omp(class NeighList *) {}
  void granular_nsq_newton_omp(class NeighList *) {}
  void granular_bin_no_newton_omp(class NeighList *) {}
  void granular_bin_newton_omp(class NeighList *) {}
  void granular_bin_newton_tri_omp(class NeighList *) {}

  void respa_nsq_no_newton_omp(class NeighList *) {}
  void respa_nsq_newton_omp(class NeighList *) {}
  void respa_bin_no_newton_omp(class NeighList *) {}
  void respa_bin_newton_omp(class NeighList *) {}
  void respa_bin_newton_tri_omp(class NeighList *) {}

#endif
#endif
