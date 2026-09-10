/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Ludwig Ahrens-Iwers (TUHH), Shern Tee (UQ), Robert Meissner (TUHH)
------------------------------------------------------------------------- */

#ifdef KSPACE_CLASS
// clang-format off
KSpaceStyle(pppm/electrode/tip4p, PPPMElectrodeTIP4P);
// clang-format on
#else

#ifndef LMP_PPPM_ELECTRODE_TIP4P_H
#define LMP_PPPM_ELECTRODE_TIP4P_H

#include "pppm_electrode.h"
#include "pppm.h"

namespace LAMMPS_NS {

class PPPMElectrodeTIP4P : public PPPMElectrode {
 public:
  PPPMElectrodeTIP4P(class LAMMPS *);
 protected:



  void compute_boundary_corr(double, int, int, double &, double *) override;
  void compute_vector_boundary_corr(double *, int, int, bool) override;
   void init_tip4p() override;
 private:
  // TIP4P support: atom types and geometric parameters for the fictitious M-site

  void particle_map() override;
  void make_rho() override;
  void fieldforce_ik() override;
  void fieldforce_ad() override;
  void fieldforce_peratom() override;


  // compute the fictitious TIP4P charge site position
  void find_M(int i, int &iH1, int &iH2, double *xM);
  void make_rho_in_brick(int, FFT_SCALAR ***, bool) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
