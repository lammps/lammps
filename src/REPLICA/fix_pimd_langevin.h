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
FixStyle(pimd/langevin,FixPIMDLangevin);
// clang-format on
#else

#ifndef FIX_PIMD_LANGEVIN_H
#define FIX_PIMD_LANGEVIN_H

#include "fix_pimd_nve.h"

namespace LAMMPS_NS {

class FixPIMDLangevin : public FixPIMDNVE {
 public:
  FixPIMDLangevin(class LAMMPS *, int, char **);
  ~FixPIMDLangevin() override;

  enum { PHYSICAL, NORMAL };
  enum { BAOAB, OBABO };
  enum { ISO, ANISO, TRICLINIC };
  enum { PILE_L };
  enum { MTTK, BZP };
  enum { NVE, NVT, NPH, NPT };
  enum { SINGLE_PROC, MULTI_PROC };

  void end_of_step() override;

 protected:
  bool parse_keyword(int, char **, int &) override;
  void setup_subclass_state() override;
  int seed;

  // System setting variables
  int thermostat;          // PILE_L
  int barostat;            // BZP or MTTK
  int ensemble;            // nve or nvt or nph or npt

  double fixedpoint[3];    // location of dilation fixed-point

  /* Langevin integration */

  double gamma, c1, c2, tau;
  double *tau_k, *c1_k, *c2_k;
  double pilescale;
  double Lan_temp;

  class RanMars *random;

  int tstat_flag;    // 1 if a thermostat is used
  void langevin_init();
  void qc_step() override;    // integrate for dt/2 for the centroid mode (x <- x + v * dt/2)
  void o_step() override;    // thermostat velocities and optional barostat momentum
  void b_step() override;

  /* Bussi-Zykova-Parrinello barostat */

  int pstat_flag;    // pstat_flag = 1 if barostat is used
  int pstyle;        // pstyle = ISO or ANISO (will support TRICLINIC in the future)
  double W, tau_p, Pext, p_hydro, totenthalpy, Vcoeff;
  int pdim;
  int p_flag[6];
  double p_target[6];
  double vw[6];               // barostat velocity
  double ke_tensor[6];        // kinetic energy tensor
  double c_vir_tensor[6];     // centroid-virial tensor
  double stress_tensor[6];    // path integral centroid-virial stress tensor

  void baro_init();
  void press_v_step();
  void press_o_step();

  /* centroid-virial estimator computation */
  double vol0 = 0.0;

  /* Langevin-specific estimators */
  void compute_stress_tensor();
  void compute_cvir() override;
  void compute_totenthalpy();
  double compute_subclass_vector(int) const override;
  int base_restart_size() const override;
  int pack_base_restart(double *) const override;
  int unpack_base_restart(const double *) override;
};
}    // namespace LAMMPS_NS
#endif
#endif
