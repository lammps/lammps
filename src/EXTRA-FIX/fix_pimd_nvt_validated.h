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
FixStyle(pimd/nvt/validated,FixPIMDNVTValidated);
// clang-format on
#else

#ifndef LMP_FIX_PIMD_NVT_VALIDATED_H
#define LMP_FIX_PIMD_NVT_VALIDATED_H

#include "fix.h"

#include <functional>

namespace LAMMPS_NS {

class Compute;
class FixPIMDNVTValidated : public Fix {
 public:
  FixPIMDNVTValidated(class LAMMPS *, int, char **);
  ~FixPIMDNVTValidated() override;

  int setmask() override;

  void init() override;
  void setup(int) override;
  void post_force(int) override;
  void initial_integrate(int) override;
  void final_integrate() override;
  void end_of_step() override;

  double compute_vector(int) override;

 protected:
  using KeywordParser = std::function<bool(int, char **, int &)>;

  FixPIMDNVTValidated(class LAMMPS *, int, char **, bool);
  void init_defaults();
  void finish_nuclear_constructor_setup();
  void parse_arguments(int, char **, const KeywordParser &);
  bool parse_common_keyword(int, char **, int &);

  // System setting variables
  int method;                              // TPPIMD or TPRPMD
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
  int thermostat;          // NHC
  int integrator;          // obabo or baoab
  int ensemble;            // NVT
  int mapflag;             // should be 1 if number of beads > 1
  int removecomflag;
  double masstotal;

  double fixedpoint[3];    // location of dilation fixed-point

  // NHC

  double *eta, *eta_dot;    // chain thermostat for particles
  double *eta_dotdot;
  double *eta_mass;

  int mtchain;                 // length of chain
  int nc_tchain;
  double factor_eta;
  double drag, tdrag_factor;     // drag factor on particle thermostat
  double t_freq;
  double t_period;
  double tdof;
  double ke_target;

  // ring-polymer model

  double omega_np, fbond, spring_energy, sp;

  // fictitious mass

  double fmass, *mass;

  // inter-partition communication

  MPI_Comm rootworld;
  int me, nprocs, ireplica, nreplica, nprocs_universe;
  int ntotal, maxlocal;

  int x_last, x_next;

  int cmode;
  int sizeplan;
  int *plansend, *planrecv;

  tagint *tagsend, *tagrecv;
  double **bufsend, **bufrecv, **bufbeads;
  double **bufsorted, **bufsortedall;

  tagint *tagsendall, *tagrecvall;
  double **bufsendall, **bufrecvall;

  int *counts, *displacements;

  void comm_init();
  void inter_replica_comm(double **ptr);

  /* normal-mode operations */

  double *lam, **M_x2xp, **M_xp2x;
  int *modeindex;

  void reallocate();
  void nmpimd_init();
  void nmpimd_transform(double **, double **, double *);

  /* Ring-polymer integration helpers */

  double dtv, dtf, dtv2, dtv3, dthalf, dt4, dt8;
  double tau;
  double *tau_k;
  double pilescale;
  double *_omega_k, *Lan_s, *Lan_c;    // sin(omega_k*dt*0.5), cos(omega_k*dt*0.5)

  int tstat_flag;    // tstat_flat = 1 if thermostat if used
  void nhc_init();
  void b_step();    // integrate for dt/2 according to B part (v <- v + f * dt/2)
  void
  a_step();    // integrate for dt/2 according to A part (non-centroid mode, harmonic force between replicas)
  void qc_step();    // integrate for dt/2 for the centroid mode (x <- x + v * dt/2)
  void o_step();     // integrate for dt according to O part (O-U process, for thermostating)
  void nhc_temp_integrate();
  double compute_nuclear_kinetic_energy() const;
  double chain_target_energy() const;
  virtual bool thermostat_chain_active() const;
  void update_chain0_acceleration(double);
  void propagate_chain_tail_halfstep(double);
  double propagate_chain0_halfstep(double, bool);
  void update_scaled_nuclear_kinetic(double &, double &) const;
  void advance_chain_positions(double);
  void complete_chain0_halfstep(double, double);
  void update_outer_chain_accelerations(double);

  virtual void setup_subclass_state();
  virtual void after_force_transform_hook();
  virtual void thermostat_step();
  virtual void force_half_step();
  virtual void centroid_position_half_step();
  virtual void nh_v_temp();
  bool nuclear_thermostat_off() const;
  double chain0_target_energy() const;
  int nuclear_restart_size() const;
  int pack_nuclear_restart(double *) const;
  int unpack_nuclear_restart(const double *);
  virtual int subclass_restart_size() const;
  virtual int pack_subclass_restart(double *, int) const;
  virtual int unpack_subclass_restart(const double *, int);
  int nuclear_vector_size() const;
  double compute_nuclear_vector(int) const;
  virtual int subclass_vector_size() const;
  virtual double compute_subclass_vector(int) const;

  /* centroid-virial estimator computation */
  double **xc, *xcall;
  int maxxc;
  int maxunwrap;
  double **x_unwrap;
  void reallocate_x_unwrap();
  void reallocate_xc();
  void collect_xc();
  void remove_com_motion();
  double vir, vir_, centroid_vir;
  double t_prim, t_vir, t_cv, p_prim, p_cv, p_md;

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
  void compute_t_prim();
  void compute_t_vir();
  void compute_p_prim();
  void compute_p_cv();    // centroid-virial pressure estimator
  void compute_vir();
  void compute_xf_vir();
  void compute_cvir();
  void write_restart(FILE *fp) override;
  int size_restart_global();
  int pack_restart_data(double *list);
  void restart(char *buf) override;
};
}    // namespace LAMMPS_NS
#endif
#endif
