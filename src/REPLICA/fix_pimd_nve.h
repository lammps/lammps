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
FixStyle(pimd/nve,FixPIMDNVE);
// clang-format on
#else

#ifndef LMP_FIX_PIMD_NVE_H
#define LMP_FIX_PIMD_NVE_H

#include "fix.h"

#include <functional>

namespace LAMMPS_NS {

class Compute;

class FixPIMDNVE : public Fix {
 public:
  enum { PIMD, NMPIMD, CMD };

  FixPIMDNVE(class LAMMPS *, int, char **);
  ~FixPIMDNVE() override;

  int setmask() override;

  void init() override;
  void setup(int) override;
  void post_force(int) override;
  void initial_integrate(int) override;
  void final_integrate() override;
  void end_of_step() override;

  double compute_vector(int) override;
  void write_restart(FILE *) override;
  void restart(char *) override;

 protected:
  using KeywordParser = std::function<bool(int, char **, int &)>;

  FixPIMDNVE(class LAMMPS *, int, char **, bool);
  void init_defaults();
  void finish_constructor_setup();
  void parse_arguments(int, char **, const KeywordParser &);
  bool parse_common_keyword(int, char **, int &);

  int method;
  int integrator;
  int fmmode;
  int np;
  double inverse_np;
  double temp;
  double hbar;
  double lj_epsilon, lj_sigma, lj_mass;
  double other_planck;
  double other_mvv2e;
  double kt;
  double beta, beta_np;
  int mapflag;
  int removecomflag;
  double masstotal;
  double omega_np, fbond, spring_energy, sp;
  double fmass, *mass;

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

  double *lam, **M_x2xp, **M_xp2x;
  int *modeindex;

  double dtv, dtf, dtv2, dtv3;
  double *_omega_k, *Lan_s, *Lan_c;

  double **xc, *xcall;
  int maxxc;
  int maxunwrap;
  double **x_unwrap;

  double vir, vir_, centroid_vir;
  double t_prim, t_vir, t_cv, p_prim, p_cv, p_md;
  double pote, tote, totke;
  double ke_bead, se_bead, pe_bead;
  double total_spring_energy;
  char *id_pe;
  char *id_press;
  class Compute *c_pe;
  class Compute *c_press;

  void comm_init();
  virtual void inter_replica_comm(double **);
  void reallocate();
  void reallocate_x_unwrap();
  void reallocate_xc();

  void nmpimd_init();
  void nmpimd_transform(double **, double **, double *);

  void collect_xc();
  void b_step();
  void apply_force_velocity_kick();
  void q_step();
  virtual void qc_step();
  virtual void a_step();
  virtual void spring_force();
  void remove_com_motion();
  void unmap_coordinates(double **, imageint *);
  void remap_coordinates(double **, imageint *);
  virtual double **normal_mode_transform_buffer();
  void forward_normal_mode_transform(double **);
  void backward_normal_mode_transform(double **);
  void finalize_setup_normal_mode_coordinates();
  void begin_normal_mode_coordinate_propagation();
  void propagate_normal_mode_coordinate_halfstep();
  void finalize_normal_mode_coordinate_propagation();
  void prepare_common_virial_state();
  void prepare_normal_mode_forces();
  void schedule_common_computes();

  double estimator_atom_count(bool) const;
  double local_kinetic_energy_sum(bool) const;
  double local_normal_mode_spring_energy_sum(bool) const;
  double local_xf_virial_sum(bool) const;
  double local_centroid_virial_sum(bool) const;
  void reduce_bead_and_total(double, double &, double &) const;
  double reduce_partition_scalar(double) const;

  void compute_xf_vir();
  void compute_cvir();
  void compute_vir();
  void compute_totke();
  virtual void compute_spring_energy();
  void compute_pote();
  void compute_tote();
  virtual void compute_t_prim();
  void compute_t_vir();
  virtual void compute_p_prim();
  void compute_p_cv();

  virtual void setup_subclass_state();
  virtual void after_force_transform_hook();
  virtual int subclass_restart_size() const;
  virtual int pack_subclass_restart(double *, int) const;
  virtual int unpack_subclass_restart(const double *, int);
  virtual int subclass_vector_size() const;
  virtual double compute_subclass_vector(int) const;

  virtual int nuclear_vector_size() const;
  virtual double compute_nuclear_vector(int) const;
  virtual int base_restart_size() const;
  virtual int pack_base_restart(double *) const;
  virtual int unpack_base_restart(const double *);
  int size_restart_global();
  int pack_restart_data(double *);
};

}    // namespace LAMMPS_NS

#endif
#endif
