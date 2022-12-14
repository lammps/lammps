/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

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

#if defined(_OPENMP)
#include <omp.h>
#endif

namespace ReaxFF {
// exported Functions

// bond orders OpenMP

extern void Add_dBond_to_ForcesOMP(reax_system *, int, int, storage *, reax_list **);
extern void Add_dBond_to_Forces_NPTOMP(reax_system *, int, int, storage *, reax_list **);
extern int BOp_OMP(storage *, reax_list *, double, int, int, far_neighbor_data *,
                   single_body_parameters *, single_body_parameters *, two_body_parameters *, int,
                   double, double, double, double, double, double, double);

extern void BOOMP(reax_system *, storage *, reax_list **);

// bonds OpenMP

extern void BondsOMP(reax_system *, simulation_data *, storage *, reax_list **);

// forces OpenMP

extern void Compute_ForcesOMP(reax_system *, control_params *, simulation_data *, storage *,
                              reax_list **);

// hydrogen bonds

extern void Hydrogen_BondsOMP(reax_system *, control_params *, simulation_data *, storage *,
                              reax_list **);

// init md OpenMP

extern void InitializeOMP(reax_system *, control_params *, simulation_data *, storage *,
                          reax_list **, MPI_Comm);

// multi body

extern void Atom_EnergyOMP(reax_system *, simulation_data *, storage *, reax_list **);

// nonbonded

extern void vdW_Coulomb_Energy_OMP(reax_system *, control_params *, simulation_data *, storage *,
                                   reax_list **);
extern void Tabulated_vdW_Coulomb_Energy_OMP(reax_system *, control_params *, simulation_data *,
                                             storage *, reax_list **);
extern void LR_vdW_CoulombOMP(reax_system *, storage *, control_params *, int, int, double,
                              LR_data *);

// torsion angles

extern void Torsion_AnglesOMP(reax_system *, control_params *, simulation_data *, storage *,
                              reax_list **);

// valence angles

extern void Calculate_ThetaOMP(rvec, double, rvec, double, double *, double *);
extern void Calculate_dCos_ThetaOMP(rvec, double, rvec, double, rvec *, rvec *, rvec *);
extern void Valence_AnglesOMP(reax_system *, control_params *, simulation_data *, storage *,
                              reax_list **);

// OpenMP helpers

inline int get_tid()
{
#if defined(_OPENMP)
  return omp_get_thread_num();
#else
  return 0;
#endif
}
}    // namespace ReaxFF

#endif
