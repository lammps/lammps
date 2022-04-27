/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
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

#ifndef LMP_REAXFF_API_H
#define LMP_REAXFF_API_H

#include "reaxff_types.h"    // IWYU pragma: export

#include <cmath>

namespace ReaxFF {
// per instance data
struct API {
  control_params *control;
  reax_system *system;
  simulation_data *data;
  storage *workspace;
  reax_list *lists;
};

// exported Functions

// allocate

extern void Allocate_Workspace(control_params *, storage *, int);
extern void DeAllocate_System(reax_system *);
extern void DeAllocate_Workspace(control_params *, storage *);
extern void PreAllocate_Space(reax_system *, storage *);
extern void ReAllocate(reax_system *, control_params *, simulation_data *, storage *, reax_list **);

// bond orders

extern void BO(reax_system *, storage *, reax_list **);
extern int BOp(storage *, reax_list *, double, int, int, far_neighbor_data *,
               single_body_parameters *, single_body_parameters *, two_body_parameters *);
extern void Add_dBond_to_Forces(reax_system *, int, int, storage *, reax_list **);

// bonds

extern void Bonds(reax_system *, simulation_data *, storage *, reax_list **);

// control

extern void Read_Control_File(const char *, control_params *);

// ffield

extern void Read_Force_Field(const char *, reax_interaction *, control_params *, MPI_Comm);

// forces

extern void Compute_Forces(reax_system *, control_params *, simulation_data *, storage *,
                           reax_list **);
extern void Estimate_Storages(reax_system *, control_params *, reax_list **, int *, int *, int *,
                              int *);

// hydrogen bonds

extern void Hydrogen_Bonds(reax_system *, simulation_data *, storage *, reax_list **);

// init md

extern void Init_System(reax_system *, control_params *);
extern void Init_Simulation_Data(simulation_data *);
extern void Init_Workspace(reax_system *, control_params *, storage *);
extern void Initialize(reax_system *, control_params *, simulation_data *, storage *, reax_list **,
                       MPI_Comm);

// lists

extern void Make_List(int, int, int, reax_list *);
extern void Delete_List(reax_list *);

inline int Start_Index(int i, reax_list *l)
{
  return l->index[i];
}
inline int End_Index(int i, reax_list *l)
{
  return l->end_index[i];
}
inline void Set_Start_Index(int i, int val, reax_list *l)
{
  l->index[i] = val;
}
inline void Set_End_Index(int i, int val, reax_list *l)
{
  l->end_index[i] = val;
}
inline int Num_Entries(int i, reax_list *l)
{
  return l->end_index[i] - l->index[i];
}

// lookup

extern void Init_Lookup_Tables(reax_system *, control_params *, storage *, MPI_Comm);
extern void Deallocate_Lookup_Tables(reax_system *);
extern void Natural_Cubic_Spline(LAMMPS_NS::Error *, const double *, const double *,
                                 cubic_spline_coef *, unsigned int);
extern void Complete_Cubic_Spline(LAMMPS_NS::Error *, const double *, const double *, double v0,
                                  double vlast, cubic_spline_coef *coef, unsigned int n);

// multi body

extern void Atom_Energy(reax_system *, control_params *, simulation_data *, storage *,
                        reax_list **);

// nonbonded

extern void Compute_Polarization_Energy(reax_system *, simulation_data *, storage *);
extern void vdW_Coulomb_Energy(reax_system *, control_params *, simulation_data *, storage *,
                               reax_list **);
extern void Tabulated_vdW_Coulomb_Energy(reax_system *, control_params *, simulation_data *,
                                         storage *, reax_list **);
extern void LR_vdW_Coulomb(reax_system *, storage *, control_params *, int, int, double, LR_data *);

// reset tools

extern void Reset(reax_system *, control_params *, simulation_data *, storage *, reax_list **);
extern void Reset_Simulation_Data(simulation_data *);
extern void Reset_Workspace(reax_system *, storage *);

// toolbox

extern void *scalloc(LAMMPS_NS::Error *, rc_bigint, rc_bigint, const std::string &);
extern void *smalloc(LAMMPS_NS::Error *, rc_bigint, const std::string &);
extern void sfree(LAMMPS_NS::Error *, void *, const std::string &);

// torsion angles

extern double Calculate_Omega(rvec, double, rvec, double, rvec, double, rvec, double,
                              three_body_interaction_data *, three_body_interaction_data *, rvec,
                              rvec, rvec, rvec);
extern void Torsion_Angles(reax_system *, control_params *, simulation_data *, storage *,
                           reax_list **);

// valence angles

extern void Calculate_Theta(rvec, double, rvec, double, double *, double *);
extern void Calculate_dCos_Theta(rvec, double, rvec, double, rvec *, rvec *, rvec *);
extern void Valence_Angles(reax_system *, control_params *, simulation_data *, storage *,
                           reax_list **);

// vector

inline void rvec_Add(rvec ret, rvec v)
{
  ret[0] += v[0];
  ret[1] += v[1];
  ret[2] += v[2];
}

inline void rvec_Copy(rvec dest, rvec src)
{
  dest[0] = src[0];
  dest[1] = src[1];
  dest[2] = src[2];
}

inline void rvec_Cross(rvec ret, rvec v1, rvec v2)
{
  ret[0] = v1[1] * v2[2] - v1[2] * v2[1];
  ret[1] = v1[2] * v2[0] - v1[0] * v2[2];
  ret[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

inline double rvec_Dot(rvec v1, rvec v2)
{
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

inline double rvec_Norm(rvec v)
{
  return sqrt(SQR(v[0]) + SQR(v[1]) + SQR(v[2]));
}

inline void rvec_Scale(rvec ret, double c, rvec v)
{
  ret[0] = c * v[0];
  ret[1] = c * v[1];
  ret[2] = c * v[2];
}

inline void rvec_ScaledAdd(rvec ret, double c, rvec v)
{
  ret[0] += c * v[0], ret[1] += c * v[1], ret[2] += c * v[2];
}

inline void rvec_ScaledSum(rvec ret, double c1, rvec v1, double c2, rvec v2)
{
  ret[0] = c1 * v1[0] + c2 * v2[0];
  ret[1] = c1 * v1[1] + c2 * v2[1];
  ret[2] = c1 * v1[2] + c2 * v2[2];
}

inline void ivec_MakeZero(ivec v)
{
  v[0] = v[1] = v[2] = 0;
}

inline void ivec_Copy(ivec dest, ivec src)
{
  dest[0] = src[0], dest[1] = src[1], dest[2] = src[2];
}

inline void ivec_Scale(ivec dest, double C, ivec src)
{
  dest[0] = (int) (C * src[0]);
  dest[1] = (int) (C * src[1]);
  dest[2] = (int) (C * src[2]);
}

inline void ivec_Sum(ivec dest, ivec v1, ivec v2)
{
  dest[0] = v1[0] + v2[0];
  dest[1] = v1[1] + v2[1];
  dest[2] = v1[2] + v2[2];
}
}    // namespace ReaxFF

#endif
