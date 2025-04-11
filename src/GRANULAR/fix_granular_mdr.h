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
// clang-format off
FixStyle(GRANULAR/MDR,FixGranularMDR);
// clang-format on
#else

#ifndef LMP_FIX_GRANULAR_MDR_H
#define LMP_FIX_GRANULAR_MDR_H

#include "fix.h"

namespace LAMMPS_NS {
namespace Granular_MDR_NS {

  enum HistoryIndex {
    DELTA_0 = 0,       // apparent overlap
    DELTA_1,
    DELTAO_0,          // displacement
    DELTAO_1,
    DELTA_MDR_0,       // MDR apparent overlap
    DELTA_MDR_1,
    DELTA_BULK_0,      // bulk displacement
    DELTA_BULK_1,
    DELTAMAX_MDR_0,    // maximum MDR apparent overlap
    DELTAMAX_MDR_1,
    YFLAG_0,           // yield flag
    YFLAG_1,
    DELTAY_0,          // yield displacement
    DELTAY_1,
    CA_0,              // contact area intercept
    CA_1,
    AADH_0,            // adhesive contact radius
    AADH_1,
    AC_0,              // contact area
    AC_1,
    EPS_BAR_0,         // volume-averaged   infinitesimal sor
    EPS_BAR_1,
    PENALTY,           // contact penalty
    DELTA_MAX,
    DELTAP_0,
    DELTAP_1
  };

}    // namespace Granular_MDR_NS

class FixGranularMDR : public Fix {
 public:
  FixGranularMDR(class LAMMPS *, int, char **);
  ~FixGranularMDR() override;
  int setmask() override;
  void post_constructor() override;
  void setup_pre_force(int) override;
  void pre_force(int) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  void set_arrays(int) override;

 private:
  int comm_stage;
  char *id_fix;
  double psi_b_coeff;
  class PairGranular *pair;
  class FixNeighHistory *fix_history;
  std::vector<Fix *> fix_wall_list;

  void mean_surf_disp();
  void calculate_contact_penalty();
  void update_fix_gran_wall();

  int index_Ro;                 // initial radius
  int index_Vgeo;               // geometric particle volume of apparent particle afterremoving spherical cap volume
  int index_Velas;              // particle volume from linear elasticity
  int index_Vcaps;              // spherical cap volume from intersection of apparentradius particle and contact planes
  int index_eps_bar;            // volume-averaged infinitesimal strain tensor
  int index_dRnumerator;        // summation of numerator terms in calculation of dR
  int index_dRdenominator;      // summation of denominator terms in calculation of dR
  int index_Acon0;              // total area involved in contacts: Acon^{n}
  int index_Acon1;              // total area involved in contacts: Acon^{n+1}
  int index_Atot;               // total particle area
  int index_Atot_sum;           // running sum of contact area minus cap area
  int index_ddelta_bar;         // change in mean surface displacement
  int index_psi;                // ratio of free surface area to total surface area
  int index_sigmaxx;            // xx-component of the stress tensor, not necessary forforce calculation
  int index_sigmayy;            // yy-component of the stress tensor, not necessary forforce calculation
  int index_sigmazz;            // zz-component of the stress tensor, not necessary forforce calculation
  int index_history_setup_flag; // flag to check if history variables have beeninitialized
  int index_contacts;           // total contacts on particle
  int index_adhesive_length;    // total length of adhesive contact on a particle
};

}    // namespace LAMMPS_NS

#endif
#endif
