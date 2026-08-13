/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(xrd,ComputeXRD);
// clang-format on
#else

#ifndef LMP_COMPUTE_XRD_H
#define LMP_COMPUTE_XRD_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeXRD : public Compute {
 public:
  ComputeXRD(class LAMMPS *, int, char **);
  ~ComputeXRD() override;
  void init() override;
  void compute_array() override;
  double memory_usage() override;

 protected:
  void set_spacing();
  int update_reciprocal();
  void refresh_angles();

  int me;
  int *ztype;           // Atomic number of the different atom types
  double Min2Theta;     // Minimum 2theta value (input in 2theta rad)
  double Max2Theta;     // Maximum 2theta value (input in 2theta rad)
  double Kmax;          // Maximum reciprocal distance to explore
  double c[3];          // Resolution parameters for reciprocal space explored
  int Knmax[3];         // maximum integer value for K points in each dimension
  double dK[3];         // Parameters controlling resolution of reciprocal space explored
  double prd_inv[3];    // Inverse spacing of unit cell
  double dK_orig[3];    // Spacing that fixed the set of nodes explored
  double prd_last[3];   // Box dimensions dK was last computed from
  int warned_range;     // 1 once the out of range warning has been given
  int LP;               // Switch to turn on Lorentz-Polarization factor 1=on
  bool echo;            // echo compute_array progress
  bool manual;          // Turn on manual recpiprocal map

  int ntypes;
  int nlocalgroup;
  double lambda;    // Radiation wavelenght (distance units)
  int radflag;
  int *store_tmp;

  // parsed here so that the derived compute xrd/fft style, which builds on this
  // constructor to guarantee an identical set of reciprocal lattice nodes, can
  // accept its own keywords.  Unused by compute xrd itself.

  int nufft_order;            // width of the Kaiser-Bessel spreading stencil
  double nufft_oversample;    // oversampling factor of the FFT mesh
};

}    // namespace LAMMPS_NS

#endif
#endif
