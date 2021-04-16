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

#ifndef LMP_REAXFF_API_H
#define LMP_REAXFF_API_H

#include "reaxff_types.h"

#include <cmath>

namespace ReaxFF 
{
  // per instance data
  struct API {
    control_params *control;
    reax_system *system;
    output_controls *out_control;
    simulation_data *data;
    storage *workspace;
    reax_list *lists;
  };

  // exported Functions

  // allocate

  extern void DeAllocate_System(reax_system *);
  extern void DeAllocate_Workspace(control_params *, storage *);
  extern int  PreAllocate_Space(reax_system *, control_params *, storage *);
  extern void ReAllocate(reax_system *, control_params *, simulation_data *,
                         storage *, reax_list **);

  // control

  extern void Read_Control_File(const char *, control_params *, output_controls *);

  // ffield

  extern void Read_Force_Field(const char *, reax_interaction *, control_params *);

  // forces

  extern void Compute_Forces(reax_system *, control_params *, simulation_data *,
                             storage *, reax_list **, output_controls *);

  // init md

  extern void Initialize(reax_system *, control_params *, simulation_data *,
                         storage *, reax_list **, output_controls *, MPI_Comm);

  // io tools

  extern void Close_Output_Files(reax_system *, output_controls *);
  extern void Output_Results(reax_system *, control_params *, simulation_data *,
                             reax_list **, output_controls *, MPI_Comm);

  // lists

  inline int Start_Index(int i, reax_list *l) { return l->index[i]; }
  inline int End_Index(int i, reax_list *l) { return l->end_index[i]; }
  inline void Set_Start_Index(int i, int val, reax_list *l) {
    l->index[i] = val;
  }
  inline void Set_End_Index(int i, int val, reax_list *l) {
    l->end_index[i] = val;
  }
  extern void Delete_List(reax_list*);
  extern int  Make_List(int, int, int, reax_list *);

  // lookup

  extern void Deallocate_Lookup_Tables(reax_system *);
  extern void Natural_Cubic_Spline(LAMMPS_NS::Error*, const double *,
                                   const double *, cubic_spline_coef *,
                                   unsigned int);
  extern void Complete_Cubic_Spline(LAMMPS_NS::Error*, const double *,
                                    const double *, double v0, double vlast,
                                    cubic_spline_coef *coef, unsigned int n);

  // reset tools

  extern void Reset(reax_system *, control_params *, simulation_data *,
                    storage *, reax_list **);

  // toolbox

  extern void *scalloc(LAMMPS_NS::Error *, rc_bigint, rc_bigint, const char *);
  extern void *smalloc(LAMMPS_NS::Error *, rc_bigint, const char *);
  extern void sfree(LAMMPS_NS::Error *, void *, const char *);

  // vector

  inline void rvec_Add(rvec ret, rvec v) {
    ret[0] += v[0]; ret[1] += v[1]; ret[2] += v[2];
  }

  inline void rvec_Copy(rvec dest, rvec src) {
    dest[0] = src[0]; dest[1] = src[1]; dest[2] = src[2];
  }

  inline void rvec_Cross(rvec ret, rvec v1, rvec v2)
  {
    ret[0] = v1[1] * v2[2] - v1[2] * v2[1];
    ret[1] = v1[2] * v2[0] - v1[0] * v2[2];
    ret[2] = v1[0] * v2[1] - v1[1] * v2[0];
  }

  inline double rvec_Dot(rvec v1, rvec v2) {
    return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
  }

  inline double rvec_Norm( rvec v ) {
    return sqrt( SQR(v[0]) + SQR(v[1]) + SQR(v[2]) );
  }

  inline void rvec_Scale(rvec ret, double c, rvec v) {
    ret[0] = c * v[0]; ret[1] = c * v[1]; ret[2] = c * v[2];
  }

  inline void rvec_ScaledAdd(rvec ret, double c, rvec v) {
    ret[0] += c * v[0], ret[1] += c * v[1], ret[2] += c * v[2];
  }

  inline void rvec_ScaledSum(rvec ret, double c1, rvec v1 ,double c2, rvec v2)
  {
    ret[0] = c1 * v1[0] + c2 * v2[0];
    ret[1] = c1 * v1[1] + c2 * v2[1];
    ret[2] = c1 * v1[2] + c2 * v2[2];
  }

  inline void ivec_MakeZero(ivec v) { v[0] = v[1] = v[2] = 0; }
}

#endif
