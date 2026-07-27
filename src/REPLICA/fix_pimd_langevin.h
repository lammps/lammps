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
  double lj_epsilon, lj_sigma, lj_mass;    // LJ unit energy, length, and mass scales
  double other_planck;
  double other_mvv2e;
  int thermostat;          // NHC or PILE_L
  int barostat;            // BZP
  int integrator;          // obabo or baoab
  int ensemble;            // nve or nvt or nph or npt
  int removecomflag;

  double fixedpoint[3];    // location of dilation fixed-point

  // ring-polymer model

  double spring_energy;

  // fictitious mass

  void comm_init();
  void comm_init_multirank();
  bool use_base_single_rank_comm() const;
  bool use_langevin_multirank_comm() const;
  virtual void prepare_coordinates();
  double **normal_mode_transform_buffer() override;
  void inter_replica_comm(double **ptr) override;
  void inter_replica_comm_multirank(double **ptr);
  void ring_collect(const std::vector<tagint> &miss_tag,
                                            double **ptr,
                                            std::vector<tagint> &rep_tag,
                                            std::vector<double> &rep_val);
  void spring_force() override;

  /* normal-mode operations */

  double **M_f2fp, **M_fp2f;

  void reallocate();
  void reallocate_multirank();

  int maxsend;
  int multirank_sizeplan;
  int *multirank_plansend, *multirank_planrecv;
  int *multirank_modeindex;
  tagint *multirank_tagsend;
  double *multirank_bufsend, *multirank_bufrecv;
  double **multirank_bufbeads;

  /* Langevin integration */

  double gamma, c1, c2, tau;
  double *tau_k, *c1_k, *c2_k;
  double pilescale;
  double Lan_temp;

  class RanMars *random;

  int tstat_flag;    // tstat_flat = 1 if thermostat if used
  void langevin_init();
  void b_step();    // integrate for dt/2 according to B part (v <- v + f * dt/2)
  void
  a_step();    // integrate for dt/2 according to A part (non-centroid mode, harmonic force between replicas)
  void qc_step();    // integrate for dt/2 for the centroid mode (x <- x + v * dt/2)
  void o_step();     // integrate for dt according to O part (O-U process, for thermostating)
  void q_step();     // integrate for dt/2 for all the beads (x <- x + v * dt/2)

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
  void remove_com_motion();
  double vir, vir_, centroid_vir;
  double p_vir;

  /* Computes */
  char *id_pe;
  char *id_press;
  class Compute *c_pe;
  class Compute *c_press;

  void compute_totke();                    // 1: kinetic energy
  virtual void compute_spring_energy();    // 2: spring elastic energy
  void compute_pote();                     // 3: potential energy
  void compute_tote();                     // 4: total energy: 1+2+3 for all the beads
  void compute_stress_tensor();
  virtual void compute_t_prim();
  void compute_t_vir();
  void compute_t_cv();
  void compute_p_prim();
  void compute_p_cv();    // centroid-virial pressure estimator
  void compute_vir();
  void compute_xf_vir();
  void compute_cvir();
  void compute_totenthalpy();
  void schedule_common_computes();
  int subclass_vector_size() const override;
  double compute_subclass_vector(int) const override;
  int base_restart_size() const override;
  int pack_base_restart(double *) const override;
  int unpack_base_restart(const double *) override;
};
}    // namespace LAMMPS_NS
#endif
#endif
