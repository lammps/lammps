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

#ifdef KSPACE_CLASS
// clang-format off
KSpaceStyle(msm,MSM);
// clang-format on
#else

#ifndef LMP_MSM_H
#define LMP_MSM_H

#include "kspace.h"

namespace LAMMPS_NS {

class MSM : public KSpace {
 public:
  MSM(class LAMMPS *);
  ~MSM() override;
  void init() override;
  void setup() override;
  void settings(int, char **) override;
  void compute(int, int) override;
  double memory_usage() override;

 protected:
  int me, nprocs;
  double precision;
  int nfactors;
  int *factors;
  double qqrd2e;
  double cutoff;
  double volume;
  double *delxinv, *delyinv, *delzinv;
  double h_x, h_y, h_z;
  double C_p;

  int *nx_msm, *ny_msm, *nz_msm;
  int *nxlo_in, *nylo_in, *nzlo_in;
  int *nxhi_in, *nyhi_in, *nzhi_in;
  int *nxlo_out, *nylo_out, *nzlo_out;
  int *nxhi_out, *nyhi_out, *nzhi_out;
  int *ngrid, *active_flag;
  int *alpha, *betax, *betay, *betaz;
  int nxlo_out_all, nylo_out_all, nzlo_out_all;
  int nxhi_out_all, nyhi_out_all, nzhi_out_all;
  int nxlo_direct, nxhi_direct, nylo_direct;
  int nyhi_direct, nzlo_direct, nzhi_direct;
  int nmax_direct;
  int nlower, nupper;
  int peratom_allocate_flag;
  int levels;

  MPI_Comm *world_levels;

  double ****qgrid;
  double ****egrid;
  double ****v0grid, ****v1grid, ****v2grid;
  double ****v3grid, ****v4grid, ****v5grid;
  double **g_direct;
  double **v0_direct, **v1_direct, **v2_direct;
  double **v3_direct, **v4_direct, **v5_direct;
  double *g_direct_top;
  double *v0_direct_top, *v1_direct_top, *v2_direct_top;
  double *v3_direct_top, *v4_direct_top, *v5_direct_top;

  double **phi1d, **dphi1d;

  int procgrid[3];            // procs assigned in each dim of 3d grid
  int myloc[3];               // which proc I am in each dim
  int ***procneigh_levels;    // my 6 neighboring procs, 0/1 = left/right

  class Grid3d *gcall;        // GridComm class for finest level grid
  class Grid3d **gc;          // GridComm classes for each hierarchical level

  double *gcall_buf1, *gcall_buf2;
  double **gc_buf1, **gc_buf2;
  int ngcall_buf1, ngcall_buf2, npergrid;
  int *ngc_buf1, *ngc_buf2;

  int current_level;

  int **part2grid;    // storage for particle -> grid mapping
  int nmax;

  int triclinic;
  double *boxlo;

  void set_grid_global();
  void set_proc_grid(int);
  void set_grid_local();
  void reset_grid() override;
  double estimate_1d_error(double, double);
  double estimate_3d_error();
  double estimate_total_error();
  void allocate();
  void allocate_peratom();
  void deallocate();
  void deallocate_peratom();
  void allocate_levels();
  void deallocate_levels();
  int factorable(int, int &, int &);
  virtual void particle_map();
  virtual void make_rho();
  virtual void direct(int);
  void direct_peratom(int);
  void direct_top(int);
  void direct_peratom_top(int);
  void restriction(int);
  void prolongation(int);
  void grid_swap_forward(int, double ***&);
  void grid_swap_reverse(int, double ***&);
  virtual void fieldforce();
  virtual void fieldforce_peratom();
  void compute_phis(const double &, const double &, const double &);
  void compute_phis_and_dphis(const double &, const double &, const double &);
  inline double compute_phi(const double &);
  inline double compute_dphi(const double &);
  void get_g_direct();
  void get_virial_direct();
  void get_g_direct_top(int);
  void get_virial_direct_top(int);

  // grid communication

  void pack_forward_grid(int, void *, int, int *) override;
  void unpack_forward_grid(int, void *, int, int *) override;
  void pack_reverse_grid(int, void *, int, int *) override;
  void unpack_reverse_grid(int, void *, int, int *) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
