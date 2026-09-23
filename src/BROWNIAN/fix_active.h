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

#ifdef FIX_CLASS

FixStyle(active, FixActive);

#else

#ifndef LMP_FIX_ACTIVE_H
#define LMP_FIX_ACTIVE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixActive : public Fix {
 public:
  FixActive(class LAMMPS *, int, char **);
  ~FixActive() override;
  int setmask() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void setup(int) override;
  void post_force(int) override;
  double memory_usage() override;

  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;

  // Forward-communicate per-atom orientation to ghost atoms so that the
  // Vicsek and nematic alignment kernels see the orientations of neighbours
  // owned by other MPI ranks.  Without this, alignment partners across
  // subdomain boundaries are silently dropped.
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

 private:

  enum ActiveModel {
    ABP,
    AOUP,
    CHIRAL,
    RTP,
    VICSEK,
    ROD,
    SPINNER,
    IABP,
    MIPS,
    NEMATIC
  };
  ActiveModel model;

  // common parameters
  double v0;
  double Dr;
  double Dt;
  int seed;
  int dimension;

  // model specific
  double omega;          // chiral angular velocity
  double tumble_rate;    // RTP Poisson rate
  double Rcut;           // alignment / density cutoff (Vicsek, MIPS, Nematic)
  double alignment;      // Vicsek alignment strength in [0,1]
  double aspect;         // rod mobility anisotropy mu_par / mu_perp
  double omega_spin;     // spinner angular velocity
  double eta_odd;        // odd (Hall) viscosity
  double gamma;          // friction coefficient
  double tau;            // AOUP persistence time
  double Da;             // AOUP active diffusion amplitude
  double activity;       // active nematic stress amplitude
  double K;              // nematic alignment stiffness
  double kT;             // thermal energy for IABP fluctuation-dissipation
  double rho_max;        // MIPS density at which v(rho) -> 0

  // per-atom storage
  int nmax_alloc;
  double **theta;        // 2D orientation angle [nmax][1]
  double **quat;         // 3D quaternion [nmax][4]
  double **va;           // AOUP active velocity [nmax][3]

  class RanMars *random;
  class NeighList *list;

  // model implementations
  void compute_abp();
  void compute_aoup();
  void compute_chiral();
  void compute_rtp();
  void compute_vicsek();
  void compute_rod();
  void compute_spinner();
  void compute_iabp();
  void compute_mips();
  void compute_nematic();
};

}

#endif
#endif
