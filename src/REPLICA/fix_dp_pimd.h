/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(dp_pimd,FixDPPimd)

#else

#ifndef FIX_DP_PIMD_H
#define FIX_DP_PIMD_H

#include "fix.h"
#include <cstdio>
#include <fstream>

namespace LAMMPS_NS {

class FixDPPimd : public Fix {
 public:
  FixDPPimd(class LAMMPS *, int, char **);
  virtual ~FixDPPimd();

  int setmask();

  void init();
  //void setup_pre_force(int);
  //void setup_pre_exchange();
  void setup(int);
  void post_force(int);
  void initial_integrate(int);
  void final_integrate();
  void end_of_step();

  double compute_vector(int);
//   double compute_scalar();

  // std::fstream flog;
  // FILE * flog;
  int method, fmmode;
  int np;
  double inverse_np;
  double temp, hbar, kBT, beta, beta_np;
  int thermostat;
  int barostat;
  int integrator;
  int ensemble;
  int mapflag;
  int removecomflag;
  void remove_com_motion(); double masstotal;

  /* ring-polymer model */

  double omega_np, fbond, spring_energy, sp;
  int x_last, x_next;

  // void spring_force();
  void compute_spring_energy();

  /* fictitious mass */

  double fmass, *mass;

  /* inter-partition communication */
  double *sorted;

  // int max_nsend;
  // tagint *tag_search, *tag_send, *tag_recv;
  // double *buf_send, *buf_recv, **buf_beads;
  double *buf_beads;

  // int max_nlocal;

  // int size_plan;
  // int *plan_send, *plan_recv;
  double **comm_ptr;

  void comm_init();
  void comm_exec(double **);

  // double **coords, **forces;
  // int nsearch, nfound, nsend, nrecv;
  // tagint* tags_send;
  // double *coords_send, *coords_recv;
  // double *forces_send, *forces_recv;
  
  // void comm_coords();
  // void comm_forces();

  /* centroid-virial estimator computation */
  double inv_volume = 0.0, vol_ = 0.0, vol0 = 0.0;
  double volume = 0.0;
  double *xc, *fc;
  int n_unwrap;
  double *x_unwrap;
  void update_x_unwrap();
  void compute_xc();
  // void compute_fc();
  void compute_vir();
  void compute_vir_();
  double xf, vir, xcfc, centroid_vir, t_vir, t_cv, p_vir, p_cv, p_cv_, p_md;
  double vir_, xcf, vir2;
  
  void compute_t_vir();  // centroid-virial kinetic energy estimator
  void compute_p_cv();  // centroid-virial pressure estimator
  void compute_p_vir();  // centroid-virial pressure estimator

  /* primitive kinetic energy estimator computation */
  double total_spring_energy;
  double t_prim, p_prim;

  void compute_t_prim();  // primitive kinetic energy estimator
  void compute_p_prim();  // primitive pressure estimator

  /* normal-mode operations */

  double *lam, **M_x2xp, **M_xp2x;
  int *mode_index;

  void nmpimd_init();
  void nmpimd_fill(double**);
  void nmpimd_transform(double*, double**, double*);

  /* Langevin thermostat BAOAB integration */

  double dtv, dtf, dtv2, dtv3;
  bool baoab_ready;
  double gamma, c1, c2, tau;
  double *tau_k, *c1_k, *c2_k;
  double pilescale=1.0;
  double Lan_temp;

  class RanMars *random;
  int seed=975481;
  FILE *frand;

  void Langevin_init();
  // void baoab_update_v();
  // void baoab_update_x();
  // void random_v();
  void b_step(); // integrate the dynamics for dt/2 time according to B part (external force)
  void a_step(); // integrate the dynamics for dt/2 time according to A part (harmonic force between replicas)
  void o_step(); // integrate the dynamics for dt time according to O part (O-U process, for thermostating)
  void svr_step(MPI_Comm which); // stochastic velocity rescaling thermostat
 
  double r1, r2, r3;
  double **eta;
  double ke_centroid, alpha2, sgn_, sgn, alpha;

  double _omega_np, *_omega_k, *Lan_s, *Lan_c; // sin(omega_k*dt*0.5), cos(omega_k*dt*0.5)
  
  /* Bussi-Zykova-Parrinello barostat */

  double *omega_dot;
  double f_omega, mtk_term1;
  void press_v_step();
  void compute_stress_tensor();
  // void v_press_step();
  // void x_press_step();
  void press_remap();
  //void press_x_step();
  void qc_step();
  void press_o_step();
  int pextflag, pstyle;
  double W, tau_p, Pext, totenthalpy = 0.0, Vcoeff;
  double vw[6];
  double ke_tensor[6];
  double c_vir_tensor[6];
  double stress_tensor[6];
  void compute_totenthalpy();

  /* harmonic oscillator model system */
  int harmonicflag;
  double omega;

  /* potential energy and total energy of the extended system */
  double kine, pote, tote, totke;
  double ke_bead, se_bead, pe_bead, pot_energy_partition;
  
  void compute_totke();
  void compute_pote();
  void compute_tote();

  /* add a new compute pe to compute pote */
  char *id_temp;
  char *id_pe;
  char *id_press;
  char *id_press2;
  class Compute *c_temp;
  class Compute *c_pe;
  class Compute *c_press;
  class Compute *c_press2;
  double virial[6];
  double tdof, t_current;

  /* thermodynamic integration */
  int tiflag;
  int timethod;
  double lambda, dfdl;
  double **x_scaled;
  
  // void compute_xscaled();
  // void compute_dfdl();

};


}

#endif
#endif
