/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(saed/vtk,FixSAEDVTK);
// clang-format on
#else

#ifndef LMP_FIX_SAED_VTK_H
#define LMP_FIX_SAED_VTK_H

#include "fix.h"

namespace LAMMPS_NS {

class FixSAEDVTK : public Fix {
 public:
  FixSAEDVTK(class LAMMPS *, int, char **);
  ~FixSAEDVTK() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void end_of_step() override;
  double compute_vector(int) override;
  void reset_timestep(bigint);

 private:
  int nrepeat, nfreq, irepeat;
  bigint nvalid;
  char *ids;
  FILE *fp;
  int nrows;

  int ave, nwindow, nsum, startstep;

  int norm, iwindow, window_limit;
  double *vector;
  double *vector_total;
  double **vector_list;

  void invoke_scalar(bigint);
  void invoke_vector(bigint);
  void options(int, char **);

  bigint nextvalid();

  class ComputeSAED *compute_saed;

  double Zone[3];    // Zone axis to view SAED
  double R_Ewald;    // Radius of Ewald sphere (distance units)
  double lambda;     // Radiation wavelenght (distance units)
  double dK[3];      // spacing of reciprocal points in each dimension
  int Knmax[3];      // maximum integer value for K points in each dimension
  int Knmin[3];      // minimum integer value for K points in each dimension

  int KnSlice[6];       // min 0-2 max 2-5 hkl index using zone
  double Kmax;          // Maximum reciprocal distance to explore
  double c[3];          // Parameters controlling resolution of reciprocal space explored
  double dR_Ewald;      // Thickness of Ewald sphere slice
  double prd_inv[3];    // Inverse spacing of unit cell

  char *filename;    // user-specified file
  int nOutput;
  int Dim[3];
  bool manual;    // Turn on manual recpiprocal map
};

}    // namespace LAMMPS_NS

#endif
#endif
