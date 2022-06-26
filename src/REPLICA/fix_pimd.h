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
FixStyle(pimd,FixPIMD);
// clang-format on
#else

#ifndef FIX_PIMD_H
#define FIX_PIMD_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPIMD : public Fix {
 public:
  FixPIMD(class LAMMPS *, int, char **);
  ~FixPIMD() override;

  int setmask() override;

  void init() override;
  void setup(int) override;
  void post_force(int) override;
  void initial_integrate(int) override;
  void final_integrate() override;
  void end_of_step() override;

  double memory_usage() override;
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  int pack_restart(int, double *) override;
  void unpack_restart(int, int) override;
  int maxsize_restart() override;
  int size_restart(int) override;
  double compute_vector(int) override;

  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

  /* System setting variables */
  int method; // PIMD or NMPIMD or CMD
  int fmmode; // physical or normal
  int np; // number of beads
  double inverse_np; // 1.0/np
  double temp; // temperature
  double hbar; // Planck's constant
  double kBT; // k_B * temp
  double beta, beta_np; // beta = 1./kBT beta_np = 1./kBT/np
  int thermostat; // NHC or PILE_L
  int barostat; // BZP
  int integrator; // obabo or baoab
  int ensemble; // nve or nvt or nph or npt
  int mapflag; // should be 1 if number of beads > 1
  int removecomflag;
  double masstotal;

  void remove_com_motion(); 

  double fixedpoint[3];    // location of dilation fixed-point
  /* ring-polymer model */

  double omega_np, fbond, spring_energy, sp;
  int x_last, x_next;

  void spring_force();

  /* fictitious mass */

  double fmass, *mass;

  /* inter-partition communication */

  double *sorted;

  int max_nsend;
  tagint *tag_send;
  double *buf_send;

  int max_nlocal;
  double *buf_recv, **buf_beads;

  int size_plan;
  int *plan_send, *plan_recv;
  double **comm_ptr;

  void comm_init();
  void comm_exec(double **);

  /* normal-mode operations */

  double *lam, **M_x2xp, **M_xp2x, **M_f2fp, **M_fp2f;
  int *mode_index;

  void nmpimd_init();
  void nmpimd_fill(double **);
  void nmpimd_transform(double **, double **, double *);

  /* Nose-hoover chain integration */

  int nhc_offset_one_1, nhc_offset_one_2;
  int nhc_size_one_1, nhc_size_one_2;
  int nhc_nchain;
  bool nhc_ready;
  double nhc_temp, t_sys;

  double **nhc_eta;        /* coordinates of NH chains for ring-polymer beads */
  double **nhc_eta_dot;    /* velocities of NH chains                         */
  double **nhc_eta_dotdot; /* acceleration of NH chains                       */
  double **nhc_eta_mass;   /* mass of NH chains                               */

  void nhc_init();
  void nhc_update_v();
  void nhc_update_x();

  /* Langevin integration */

  double dtv, dtf, dtv2, dtv3;
  double gamma, c1, c2, tau;
  double *tau_k, *c1_k, *c2_k;
  double pilescale=1.0;
  double Lan_temp;
  double r1, r2, r3;
  double _omega_np, *_omega_k, *Lan_s, *Lan_c; // sin(omega_k*dt*0.5), cos(omega_k*dt*0.5)

  class RanMars *random;
  int seed=975481;
  FILE *frand;

  int tstat_flag; // tstat_flat = 1 if thermostat if used
  void Langevin_init();
  void b_step(); // integrate for dt/2 according to B part (v <- v + f * dt/2)
  void a_step(); // integrate for dt/2 according to A part (non-centroid mode, harmonic force between replicas)
  void qc_step(); // integrate for dt/2 for the centroid mode (x <- x + v * dt/2)
  void o_step(); // integrate for dt according to O part (O-U process, for thermostating)
  
  /* Bussi-Zykova-Parrinello barostat */

  double f_omega, mtk_term1;
  int pstat_flag; // pstat_flag = 1 if barostat is used
  int pstyle; // pstyle = ISO or ANISO (will support TRICLINIC in the future)
  double W, tau_p, Pext, totenthalpy = 0.0, Vcoeff;
  double vw[6]; // barostat velocity
  double ke_tensor[6]; // kinetic energy tensor
  double c_vir_tensor[6]; // centroid-virial tensor
  double stress_tensor[6]; // path integral centroid-virial stress tensor

  void baro_init();
  void press_v_step();
  void press_o_step();

  /* centroid-virial estimator computation */
  double inv_volume = 0.0, vol_ = 0.0, vol0 = 0.0;
  double volume = 0.0;
  double *xc, *fc;
  int n_unwrap;
  double *x_unwrap;
  void update_x_unwrap();
  void compute_xc();
  // void compute_fc();
  double xf, vir, xcfc, centroid_vir, t_vir, t_cv, p_vir, p_cv, p_cv_, p_md;
  double vir_, xcf, vir2;

  /* Computes */
  double kine, pote, tote, totke;
  double ke_bead, se_bead, pe_bead, pot_energy_partition;
  double total_spring_energy;
  double t_prim, p_prim;
  char *id_pe;
  char *id_press;
  class Compute *c_pe;
  class Compute *c_press;

  void compute_totke(); // 1: kinetic energy 
  void compute_spring_energy(); // 2: spring elastic energy
  void compute_pote(); // 3: potential energy
  void compute_tote(); // 4: total energy: 1+2+3 for all the beads
  void compute_t_prim();  // 5: primitive kinetic energy estimator
  void compute_p_prim();  // primitive pressure estimator
  void compute_stress_tensor();
  void compute_t_vir();  // centroid-virial kinetic energy estimator
  void compute_p_cv();  // centroid-virial pressure estimator
  void compute_p_vir();  // centroid-virial pressure estimator
  void compute_vir();
  void compute_vir_();
  void compute_totenthalpy();
};

}    // namespace LAMMPS_NS

#endif
#endif
