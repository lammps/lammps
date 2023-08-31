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

#include "fix.h"

namespace LAMMPS_NS {

class FixPIMDLangevin : public Fix {
 public:
  FixPIMDLangevin(class LAMMPS *, int, char **);
  ~FixPIMDLangevin() override;

  int setmask() override;

  void init() override;
  void setup(int) override;
  void post_force(int) override;
  void initial_integrate(int) override;
  void final_integrate() override;
  void end_of_step() override;

  double compute_vector(int) override;

 protected:
  // System setting variables
  int method;                              // PIMD or NMPIMD or CMD
  int fmmode;                              // physical or normal
  int np;                                  // number of beads
  double inverse_np;                       // 1.0/np
  double temp;                             // temperature
  double hbar;                             // Planck's constant
  double lj_epsilon, lj_sigma, lj_mass;    // LJ unit energy, length, and mass scales
  double other_planck;
  double other_mvv2e;
  double kt;               // k_B * temp
  double beta, beta_np;    // beta = 1./kBT beta_np = 1./kBT/np
  int thermostat;          // NHC or PILE_L
  int barostat;            // BZP
  int integrator;          // obabo or baoab
  int ensemble;            // nve or nvt or nph or npt
  int mapflag;             // should be 1 if number of beads > 1
  int removecomflag;
  double masstotal;

  double fixedpoint[3];    // location of dilation fixed-point

  // ring-polymer model

  double omega_np, fbond, spring_energy, sp;

  // fictitious mass

  double fmass, *mass;

  // inter-partition communication

  MPI_Comm rootworld;
  int me, nprocs, ireplica, nreplica, nprocs_universe;
  int ntotal, maxlocal;

  int cmode;
  int sizeplan;
  int *plansend, *planrecv;

  tagint *tagsend, *tagrecv;
  double **bufsend, **bufrecv, **bufbeads;
  double **bufsorted, **bufsortedall;
  double **outsorted, **buftransall;

  tagint *tagsendall, *tagrecvall;
  double **bufsendall, **bufrecvall;

  int *counts, *displacements;

  void comm_init();
  void inter_replica_comm(double **ptr);

  /* normal-mode operations */

  double *lam, **M_x2xp, **M_xp2x, **M_f2fp, **M_fp2f;
  int *modeindex;

  void reallocate();
  void nmpimd_init();
  void nmpimd_transform(double **, double **, double *);

  /* Langevin integration */

  double dtv, dtf, dtv2, dtv3;
  double gamma, c1, c2, tau;
  double *tau_k, *c1_k, *c2_k;
  double pilescale;
  double Lan_temp;
  double *_omega_k, *Lan_s, *Lan_c;    // sin(omega_k*dt*0.5), cos(omega_k*dt*0.5)

  class RanMars *random;

  int tstat_flag;    // tstat_flat = 1 if thermostat if used
  void langevin_init();
  void b_step();    // integrate for dt/2 according to B part (v <- v + f * dt/2)
  void
  a_step();    // integrate for dt/2 according to A part (non-centroid mode, harmonic force between replicas)
  void qc_step();    // integrate for dt/2 for the centroid mode (x <- x + v * dt/2)
  void o_step();     // integrate for dt according to O part (O-U process, for thermostating)

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
  double **xc, *xcall;
  int maxxc;
  int maxunwrap;
  double **x_unwrap;
  void reallocate_x_unwrap();
  void reallocate_xc();
  void collect_xc();
  void remove_com_motion();
  double vir, vir_, centroid_vir;
  double t_prim, t_vir, t_cv, p_prim, p_vir, p_cv, p_md;

  /* Computes */
  double pote, tote, totke;
  double ke_bead, se_bead, pe_bead;
  double total_spring_energy;
  char *id_pe;
  char *id_press;
  class Compute *c_pe;
  class Compute *c_press;

  void compute_totke();            // 1: kinetic energy
  void compute_spring_energy();    // 2: spring elastic energy
  void compute_pote();             // 3: potential energy
  void compute_tote();             // 4: total energy: 1+2+3 for all the beads
  void compute_stress_tensor();
  void compute_t_prim();
  void compute_t_vir();
  void compute_t_cv();
  void compute_p_prim();
  void compute_p_cv();    // centroid-virial pressure estimator
  void compute_vir();
  void compute_cvir();
  void compute_totenthalpy();

  void write_restart(FILE *fp) override;
  int size_restart_global();
  int pack_restart_data(double *list);
  void restart(char *buf) override;
};
}    // namespace LAMMPS_NS
#endif
#endif
