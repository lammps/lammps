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
FixStyle(pimd,FixPIMDNVT);
FixStyle(pimd/nvt,FixPIMDNVT);
// clang-format on
#else

#ifndef FIX_PIMD_NVT_H
#define FIX_PIMD_NVT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPIMDNVT : public Fix {
 public:
  FixPIMDNVT(class LAMMPS *, int, char **);
  ~FixPIMDNVT() override;

  int setmask() override;

  void init() override;
  void setup(int) override;
  void post_force(int) override;
  void initial_integrate(int) override;
  void final_integrate() override;

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

  int method;
  int np;
  double inverse_np;

  /* ring-polymer model */

  double omega_np, fbond, spring_energy, sp, virial;
  int x_last, x_next;

  void spring_force();

  /* fictitious mass */

  double fmass, *mass;

  /* inter-partition communication */

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
  double nhc_temp, dtv, dtf, t_sys;

  double **nhc_eta;        /* coordinates of NH chains for ring-polymer beads */
  double **nhc_eta_dot;    /* velocities of NH chains                         */
  double **nhc_eta_dotdot; /* acceleration of NH chains                       */
  double **nhc_eta_mass;   /* mass of NH chains                               */

  void nhc_init();
  void nhc_update_v();
  void nhc_update_x();
};

}    // namespace LAMMPS_NS

#endif
#endif
