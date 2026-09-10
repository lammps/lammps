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
  double h_last[6];      // Box matrix the reciprocal lattice was last built from
  int warned_range;     // 1 once the out of range warning has been given
  int triclinic;        // 1 if the simulation cell is triclinic

  // step in reciprocal space per unit of each node index, scaled by the c
  // parameters.  the rows are the reciprocal lattice vectors of the cell, so
  // for an orthogonal box only the diagonal is nonzero and equals dK.

  double rlv[3][3];

  // reciprocal space position of a node, K = i*rlv[0] + j*rlv[1] + k*rlv[2].
  // rlv is upper triangular, exactly as the LAMMPS box matrix is.

  inline void kvector(const double g[3][3], int i, int j, int k, double *K) const {
    K[0] = i*g[0][0];
    K[1] = i*g[0][1] + j*g[1][1];
    K[2] = i*g[0][2] + j*g[1][2] + k*g[2][2];
  }
  int LP;               // Switch to turn on Lorentz-Polarization factor 1=on
  bool echo;            // echo compute_array progress
  bool manual;          // Turn on manual recpiprocal map

  int ntypes;
  int nlocalgroup;
  double lambda;    // Radiation wavelenght (distance units)
  int radflag;
  int *store_tmp;

  // a keyword this style does not recognize may belong to a style derived from
  // it, which parses the command again once this constructor has finished, so
  // the leftovers are recorded rather than rejected here.  whichever class ends
  // up parsing last rejects the ones nothing claimed.

  int nunclaimed;      // number of argument indices below
  int *unclaimed;      // indices in argv of the words no style has claimed

  void reject_unclaimed(char **);
};

}    // namespace LAMMPS_NS

#endif
#endif
