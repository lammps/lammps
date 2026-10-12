/* -*- c++ -*- --------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* Public C++ interface between LAMMPS src/GPU and lib/gpu.
   All entry points live in namespace LAMMPS_GPU so that the compiler
   validates each call-site declaration against the single definition in
   the corresponding lib/gpu/lal_*_ext.cpp or lal_device.cpp.

   Prerequisites: mpi.h (for MPI_Comm) and either lmptype.h (LAMMPS side)
   or lal_precision.h (lib/gpu side) must be included before this header
   for the tagint typedef; <string> must be visible for std::string.        */

#ifndef LMP_LAMMPS_GPU_H
#define LMP_LAMMPS_GPU_H

#include <cstdio>
#include <mpi.h>
#include <string>

#if !defined(LMP_LMPTYPE_H) && !defined(LAL_PRECISION_H)
#error Must include either lal_precision.h or lmptype.h before lammps_gpu.h
#endif

namespace LAMMPS_GPU {

/* When included from the LAMMPS src side, tagint lives in LAMMPS_NS.
   Bring it into this namespace so the function declarations below can use it
   unqualified.  From the lib/gpu side (lal_precision.h), tagint is already at
   global scope and is found by ordinary lookup, so this block is skipped. */
#if defined(LMP_LMPTYPE_H) && !defined(LAL_PRECISION_H)
using LAMMPS_NS::tagint;
#endif

// Device-level functions (lal_device.cpp)
extern bool lmp_has_compatible_gpu_device();
extern bool lmp_gpu_requires_host_neighbor();
extern std::string lmp_gpu_device_info();
extern void lmp_gpu_defer_device_clear(int flag);
extern int lmp_init_device(MPI_Comm world, MPI_Comm replica, const int ngpu, const int first_gpu_id,
                           const int gpu_mode, const int t_per_atom, const double user_cell_size,
                           char *opencl_config, const int ocl_platform, char *device_type_flags,
                           const int block_pair);
extern void lmp_clear_device();
extern double lmp_gpu_forces(double **f, double **tor, double *eatom, double **vatom,
                             double *virial, double &ecoul, int &error_flag);
extern double lmp_gpu_update_bin_size(const double subx, const double suby, const double subz,
                                      const int nlocal, const double cut);
extern bool lmp_gpu_config(const std::string &category, const std::string &setting);

// Per-style entry points

// pair amoeba
extern int amoeba_gpu_init(const int ntypes, const int max_amtype, const int max_amclass,
                           const double *host_pdamp, const double *host_thole,
                           const double *host_dirdamp, const int *host_amtype2class,
                           const double *host_special_hal, const double *host_special_repel,
                           const double *host_special_disp, const double *host_special_mpole,
                           const double *host_special_polar_wscale,
                           const double *host_special_polar_piscale,
                           const double *host_special_polar_pscale, const double *host_csix,
                           const double *host_adisp, const int nlocal, const int nall,
                           const int max_nbors, const int maxspecial, const int maxspecial15,
                           const double cell_size, int &gpu_mode, FILE *screen,
                           const double polar_dscale, const double polar_uscale);
extern void amoeba_gpu_clear();
extern int **amoeba_gpu_precompute(const int ago, const int inum_full, const int nall,
                                   double **host_x, int *host_type, int *host_amtype,
                                   int *host_amgroup, double **host_rpole, double **host_uind,
                                   double **host_uinp, double *host_pval, double *sublo,
                                   double *subhi, tagint *tag, int **nspecial, tagint **special,
                                   int *nspecial15, tagint **special15, const bool eflag_in,
                                   const bool vflag_in, const bool eatom, const bool vatom,
                                   int **ilist, int **jnum, bool &success, double *host_q,
                                   double *boxlo, double *prd);
extern void amoeba_gpu_compute_multipole_real(
    const int ago, const int inum, const int nall, double **host_x, int *host_type,
    int *host_amtype, int *host_amgroup, double **host_rpole, double *sublo, double *subhi,
    tagint *tag, int **nspecial, tagint **special, int *nspecial15, tagint **special15,
    const bool eflag, const bool vflag, const bool eatom, const bool vatom, int **ilist, int **jnum,
    bool &success, const double aewald, const double felec, const double off2, double *host_q,
    double *boxlo, double *prd, void **tq_ptr);
extern void amoeba_gpu_compute_udirect2b(int *host_amtype, int *host_amgroup, double **host_rpole,
                                         double **host_uind, double **host_uinp,
                                         const double aewald, const double off2, void **fieldp_ptr);
extern void amoeba_gpu_compute_umutual2b(int *host_amtype, int *host_amgroup, double **host_rpole,
                                         double **host_uind, double **host_uinp,
                                         const double aewald, const double off2, void **fieldp_ptr);
extern void amoeba_gpu_update_fieldp(void **fieldp_ptr);
extern void amoeba_gpu_precompute_kspace(const int inum_full, const int bsorder,
                                         double ***host_thetai1, double ***host_thetai2,
                                         double ***host_thetai3, int **igrid, const int nzlo_out,
                                         const int nzhi_out, const int nylo_out, const int nyhi_out,
                                         const int nxlo_out, const int nxhi_out);
extern void amoeba_gpu_fphi_uind(double ****host_grid_brick, void **host_fdip_phi1,
                                 void **host_fdip_phi2, void **host_fdip_sum_phi);
extern void amoeba_gpu_fphi_mpole(double ***host_grid_brick, void **host_fdip_sum_phi,
                                  const double felec);
extern void amoeba_gpu_compute_polar_real(int *host_amtype, int *host_amgroup, double **host_rpole,
                                          double **host_uind, double **host_uinp, const bool eflag,
                                          const bool vflag, const bool eatom, const bool vatom,
                                          const double aewald, const double felec,
                                          const double off2, void **tq_ptr);
extern double amoeba_gpu_bytes();

// pair beck
extern int beck_gpu_init(const int ntypes, double **cutsq, double **host_aa, double **alpha,
                         double **beta, double **AA, double **BB, double *special_lj,
                         const int nlocal, const int nall, const int max_nbors,
                         const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen);
extern void beck_gpu_clear();
extern int **beck_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                                int *host_type, double *sublo, double *subhi, tagint *tag,
                                int **nspecial, tagint **special, const bool eflag,
                                const bool vflag, const bool eatom, const bool vatom, int **ilist,
                                int **jnum, bool &success, double *prd, int *periodicity);
extern void beck_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                             int *host_type, int *ilist, int *numj, int **firstneigh,
                             const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                             bool &success);
extern double beck_gpu_bytes();

// pair born/coul/long/cs
extern int bornclcs_gpu_init(const int ntypes, double **cutsq, double **host_rhoinv,
                             double **host_born1, double **host_born2, double **host_born3,
                             double **host_a, double **host_c, double **host_d, double **sigma,
                             double **offset, double *special_lj, const int inum, const int nall,
                             const int max_nbors, const int maxspecial, const double cell_size,
                             int &gpu_mode, FILE *screen, double **host_cut_ljsq,
                             double host_cut_coulsq, double *host_special_coul, const double qqrd2e,
                             const double g_ewald);
extern void bornclcs_gpu_clear();
extern int **bornclcs_gpu_compute_n(const int ago, const int inum_full, const int nall,
                                    double **host_x, int *host_type, double *sublo, double *subhi,
                                    tagint *tag, int **nspecial, tagint **special, const bool eflag,
                                    const bool vflag, const bool eatom, const bool vatom,
                                    int **ilist, int **jnum, bool &success, double *host_q,
                                    double *boxlo, double *prd, int *periodicity);
extern void bornclcs_gpu_compute(const int ago, const int inum_full, const int nall,
                                 double **host_x, int *host_type, int *ilist, int *numj,
                                 int **firstneigh, const bool eflag, const bool vflag,
                                 const bool eatom, const bool vatom, bool &success, double *host_q,
                                 const int nlocal, double *boxlo, double *prd);
extern double bornclcs_gpu_bytes();

// pair born/coul/long
extern int borncl_gpu_init(const int ntypes, double **cutsq, double **host_rhoinv,
                           double **host_born1, double **host_born2, double **host_born3,
                           double **host_a, double **host_c, double **host_d, double **sigma,
                           double **offset, double *special_lj, const int inum, const int nall,
                           const int max_nbors, const int maxspecial, const double cell_size,
                           int &gpu_mode, FILE *screen, double **host_cut_ljsq,
                           double host_cut_coulsq, double *host_special_coul, const double qqrd2e,
                           const double g_ewald);
extern void borncl_gpu_clear();
extern int **borncl_gpu_compute_n(const int ago, const int inum_full, const int nall,
                                  double **host_x, int *host_type, double *sublo, double *subhi,
                                  tagint *tag, int **nspecial, tagint **special, const bool eflag,
                                  const bool vflag, const bool eatom, const bool vatom, int **ilist,
                                  int **jnum, bool &success, double *host_q, double *boxlo,
                                  double *prd, int *periodicity);
extern void borncl_gpu_compute(const int ago, const int inum_full, const int nall, double **host_x,
                               int *host_type, int *ilist, int *numj, int **firstneigh,
                               const bool eflag, const bool vflag, const bool eatom,
                               const bool vatom, bool &success, double *host_q, const int nlocal,
                               double *boxlo, double *prd);
extern double borncl_gpu_bytes();

// pair born/coul/wolf/cs
extern int borncwcs_gpu_init(const int ntypes, double **cutsq, double **host_rhoinv,
                             double **host_born1, double **host_born2, double **host_born3,
                             double **host_a, double **host_c, double **host_d, double **sigma,
                             double **offset, double *special_lj, const int inum, const int nall,
                             const int max_nbors, const int maxspecial, const double cell_size,
                             int &gpu_mode, FILE *screen, double **host_cut_ljsq,
                             double host_cut_coulsq, double *host_special_coul, const double qqrd2e,
                             const double alf, const double e_shift, const double f_shift);
extern void borncwcs_gpu_clear();
extern int **borncwcs_gpu_compute_n(const int ago, const int inum_full, const int nall,
                                    double **host_x, int *host_type, double *sublo, double *subhi,
                                    tagint *tag, int **nspecial, tagint **special, const bool eflag,
                                    const bool vflag, const bool eatom, const bool vatom,
                                    int **ilist, int **jnum, bool &success, double *host_q,
                                    double *boxlo, double *prd, int *periodicity);
extern void borncwcs_gpu_compute(const int ago, const int inum_full, const int nall,
                                 double **host_x, int *host_type, int *ilist, int *numj,
                                 int **firstneigh, const bool eflag, const bool vflag,
                                 const bool eatom, const bool vatom, bool &success, double *host_q,
                                 const int nlocal, double *boxlo, double *prd);
extern double borncwcs_gpu_bytes();

// pair born/coul/wolf
extern int borncw_gpu_init(const int ntypes, double **cutsq, double **host_rhoinv,
                           double **host_born1, double **host_born2, double **host_born3,
                           double **host_a, double **host_c, double **host_d, double **sigma,
                           double **offset, double *special_lj, const int inum, const int nall,
                           const int max_nbors, const int maxspecial, const double cell_size,
                           int &gpu_mode, FILE *screen, double **host_cut_ljsq,
                           double host_cut_coulsq, double *host_special_coul, const double qqrd2e,
                           const double alf, const double e_shift, const double f_shift);
extern void borncw_gpu_clear();
extern int **borncw_gpu_compute_n(const int ago, const int inum_full, const int nall,
                                  double **host_x, int *host_type, double *sublo, double *subhi,
                                  tagint *tag, int **nspecial, tagint **special, const bool eflag,
                                  const bool vflag, const bool eatom, const bool vatom, int **ilist,
                                  int **jnum, bool &success, double *host_q, double *boxlo,
                                  double *prd, int *periodicity);
extern void borncw_gpu_compute(const int ago, const int inum_full, const int nall, double **host_x,
                               int *host_type, int *ilist, int *numj, int **firstneigh,
                               const bool eflag, const bool vflag, const bool eatom,
                               const bool vatom, bool &success, double *host_q, const int nlocal,
                               double *boxlo, double *prd);
extern double borncw_gpu_bytes();

// pair born
extern int born_gpu_init(const int ntypes, double **cutsq, double **host_rhoinv,
                         double **host_born1, double **host_born2, double **host_born3,
                         double **host_a, double **host_c, double **host_d, double **host_sigma,
                         double **offset, double *special_lj, const int inum, const int nall,
                         const int max_nbors, const int maxspecial, const double cell_size,
                         int &gpu_mode, FILE *screen);
extern void born_gpu_reinit(const int ntypes, double **host_rhoinv, double **host_born1,
                            double **host_born2, double **host_born3, double **host_a,
                            double **host_c, double **host_d, double **offset);
extern void born_gpu_clear();
extern int **born_gpu_compute_n(const int ago, const int inum_full, const int nall, double **host_x,
                                int *host_type, double *sublo, double *subhi, tagint *tag,
                                int **nspecial, tagint **special, const bool eflag,
                                const bool vflag, const bool eatom, const bool vatom, int **ilist,
                                int **jnum, bool &success, double *prd, int *periodicity);
extern void born_gpu_compute(const int ago, const int inum_full, const int nall, double **host_x,
                             int *host_type, int *ilist, int *numj, int **firstneigh,
                             const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                             bool &success);
extern double born_gpu_bytes();

// pair buck/coul/cut
extern int buckc_gpu_init(const int ntypes, double **cutsq, double **host_rhoinv,
                          double **host_buck1, double **host_buck2, double **host_a,
                          double **host_c, double **offset, double *special_lj, const int inum,
                          const int nall, const int max_nbors, const int maxspecial,
                          const double cell_size, int &gpu_mode, FILE *screen,
                          double **host_cut_ljsq, double **host_cut_coulsq,
                          double *host_special_coul, const double qqrd2e);
extern void buckc_gpu_clear();
extern int **buckc_gpu_compute_n(const int ago, const int inum_full, const int nall,
                                 double **host_x, int *host_type, double *sublo, double *subhi,
                                 tagint *tag, int **nspecial, tagint **special, const bool eflag,
                                 const bool vflag, const bool eatom, const bool vatom, int **ilist,
                                 int **jnum, bool &success, double *host_q, double *boxlo,
                                 double *prd, int *periodicity);
extern void buckc_gpu_compute(const int ago, const int inum_full, const int nall, double **host_x,
                              int *host_type, int *ilist, int *numj, int **firstneigh,
                              const bool eflag, const bool vflag, const bool eatom,
                              const bool vatom, bool &success, double *host_q, const int nlocal,
                              double *boxlo, double *prd);
extern double buckc_gpu_bytes();

// pair buck/coul/long
extern int buckcl_gpu_init(const int ntypes, double **cutsq, double **host_rhoinv,
                           double **host_buck1, double **host_buck2, double **host_a,
                           double **host_c, double **offset, double *special_lj, const int inum,
                           const int nall, const int max_nbors, const int maxspecial,
                           const double cell_size, int &gpu_mode, FILE *screen,
                           double **host_cut_ljsq, double host_cut_coulsq,
                           double *host_special_coul, const double qqrd2e, const double g_ewald);
extern void buckcl_gpu_clear();
extern int **buckcl_gpu_compute_n(const int ago, const int inum_full, const int nall,
                                  double **host_x, int *host_type, double *sublo, double *subhi,
                                  tagint *tag, int **nspecial, tagint **special, const bool eflag,
                                  const bool vflag, const bool eatom, const bool vatom, int **ilist,
                                  int **jnum, bool &success, double *host_q, double *boxlo,
                                  double *prd, int *periodicity);
extern void buckcl_gpu_compute(const int ago, const int inum_full, const int nall, double **host_x,
                               int *host_type, int *ilist, int *numj, int **firstneigh,
                               const bool eflag, const bool vflag, const bool eatom,
                               const bool vatom, bool &success, double *host_q, const int nlocal,
                               double *boxlo, double *prd);
extern double buckcl_gpu_bytes();

// pair buck
extern int buck_gpu_init(const int ntypes, double **cutsq, double **host_rhoinv,
                         double **host_buck1, double **host_buck2, double **host_a, double **host_c,
                         double **offset, double *special_lj, const int inum, const int nall,
                         const int max_nbors, const int maxspecial, const double cell_size,
                         int &gpu_mode, FILE *screen);
extern void buck_gpu_reinit(const int ntypes, double **cutsq, double **host_rhoinv,
                            double **host_buck1, double **host_buck2, double **host_a,
                            double **host_c, double **offset);
extern void buck_gpu_clear();
extern int **buck_gpu_compute_n(const int ago, const int inum_full, const int nall, double **host_x,
                                int *host_type, double *sublo, double *subhi, tagint *tag,
                                int **nspecial, tagint **special, const bool eflag,
                                const bool vflag, const bool eatom, const bool vatom, int **ilist,
                                int **jnum, bool &success, double *prd, int *periodicity);
extern void buck_gpu_compute(const int ago, const int inum_full, const int nall, double **host_x,
                             int *host_type, int *ilist, int *numj, int **firstneigh,
                             const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                             bool &success);
extern double buck_gpu_bytes();

// pair colloid
extern int colloid_gpu_init(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                            double **host_lj3, double **host_lj4, double **offset,
                            double *special_lj, double **host_a12, double **host_a1,
                            double **host_a2, double **host_d1, double **host_d2,
                            double **host_sigma3, double **host_sigma6, int **host_form,
                            const int nlocal, const int nall, const int max_nbors,
                            const int maxspecial, const double cell_size, int &gpu_mode,
                            FILE *screen);
extern void colloid_gpu_clear();
extern int **colloid_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                                   int *host_type, double *sublo, double *subhi, tagint *tag,
                                   int **nspecial, tagint **special, const bool eflag,
                                   const bool vflag, const bool eatom, const bool vatom,
                                   int **ilist, int **jnum, bool &success, double *prd,
                                   int *periodicity);
extern void colloid_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                                int *host_type, int *ilist, int *numj, int **firstneigh,
                                const bool eflag, const bool vflag, const bool eatom,
                                const bool vatom, bool &success);
extern double colloid_gpu_bytes();

// pair coul/cut
extern int coul_gpu_init(const int ntypes, double **host_scale, double **cutsq,
                         double *special_coul, const int nlocal, const int nall,
                         const int max_nbors, const int maxspecial, const double cell_size,
                         int &gpu_mode, FILE *screen, const double qqrd2e);
extern void coul_gpu_reinit(const int ntypes, double **host_scale);
extern void coul_gpu_clear();
extern int **coul_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                                int *host_type, double *sublo, double *subhi, tagint *tag,
                                int **nspecial, tagint **special, const bool eflag,
                                const bool vflag, const bool eatom, const bool vatom, int **ilist,
                                int **jnum, bool &success, double *host_q, double *boxlo,
                                double *prd, int *periodicity);
extern void coul_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                             int *host_type, int *ilist, int *numj, int **firstneigh,
                             const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                             bool &success, double *host_q, const int nlocal, double *boxlo,
                             double *prd);
extern double coul_gpu_bytes();

// pair coul/debye
extern int cdebye_gpu_init(const int ntypes, double **host_scale, double **cutsq,
                           double *special_coul, const int nlocal, const int nall,
                           const int max_nbors, const int maxspecial, const double cell_size,
                           int &gpu_mode, FILE *screen, const double qqrd2e, const double kappa);
extern void cdebye_gpu_reinit(const int ntypes, double **host_scale);
extern void cdebye_gpu_clear();
extern int **cdebye_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                                  int *host_type, double *sublo, double *subhi, tagint *tag,
                                  int **nspecial, tagint **special, const bool eflag,
                                  const bool vflag, const bool eatom, const bool vatom, int **ilist,
                                  int **jnum, bool &success, double *host_q, double *boxlo,
                                  double *prd, int *periodicity);
extern void cdebye_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                               int *host_type, int *ilist, int *numj, int **firstneigh,
                               const bool eflag, const bool vflag, const bool eatom,
                               const bool vatom, bool &success, double *host_q, const int nlocal,
                               double *boxlo, double *prd);
extern double cdebye_gpu_bytes();

// pair coul/dsf
extern int cdsf_gpu_init(const int ntypes, const int nlocal, const int nall, const int max_nbors,
                         const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen,
                         const double host_cut_coulsq, double *host_special_coul,
                         const double qqrd2e, const double e_shift, const double f_shift,
                         const double alpha);
extern void cdsf_gpu_clear();
extern int **cdsf_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                                int *host_type, double *sublo, double *subhi, tagint *tag,
                                int **nspecial, tagint **special, const bool eflag,
                                const bool vflag, const bool eatom, const bool vatom, int **ilist,
                                int **jnum, bool &success, double *host_q, double *boxlo,
                                double *prd, int *periodicity);
extern void cdsf_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                             int *host_type, int *ilist, int *numj, int **firstneigh,
                             const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                             bool &success, double *host_q, const int nlocal, double *boxlo,
                             double *prd);
extern double cdsf_gpu_bytes();

// pair coul/long/cs
extern int clcs_gpu_init(const int ntypes, double **scale, const int nlocal, const int nall,
                         const int max_nbors, const int maxspecial, const double cell_size,
                         int &gpu_mode, FILE *screen, double host_cut_coulsq,
                         double *host_special_coul, const double qqrd2e, const double g_ewald);
extern void clcs_gpu_reinit(const int ntypes, double **scale);
extern void clcs_gpu_clear();
extern int **clcs_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                                int *host_type, double *sublo, double *subhi, tagint *tag,
                                int **nspecial, tagint **special, const bool eflag,
                                const bool vflag, const bool eatom, const bool vatom, int **ilist,
                                int **jnum, bool &success, double *host_q, double *boxlo,
                                double *prd, int *periodicity);
extern void clcs_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                             int *host_type, int *ilist, int *numj, int **firstneigh,
                             const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                             bool &success, double *host_q, const int nlocal, double *boxlo,
                             double *prd);
extern double clcs_gpu_bytes();

// pair coul/long
extern int cl_gpu_init(const int ntypes, double **scale, const int nlocal, const int nall,
                       const int max_nbors, const int maxspecial, const double cell_size,
                       int &gpu_mode, FILE *screen, double host_cut_coulsq,
                       double *host_special_coul, const double qqrd2e, const double g_ewald);
extern void cl_gpu_reinit(const int ntypes, double **scale);
extern void cl_gpu_clear();
extern int **cl_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                              int *host_type, double *sublo, double *subhi, tagint *tag,
                              int **nspecial, tagint **special, const bool eflag, const bool vflag,
                              const bool eatom, const bool vatom, int **ilist, int **jnum,
                              bool &success, double *host_q, double *boxlo, double *prd,
                              int *periodicity);
extern void cl_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                           int *host_type, int *ilist, int *numj, int **firstneigh,
                           const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                           bool &success, double *host_q, const int nlocal, double *boxlo,
                           double *prd);
extern double cl_gpu_bytes();

// pair coul/slater/long
extern int csl_gpu_init(const int ntypes, double **scale, const int nlocal, const int nall,
                        const int max_nbors, const int maxspecial, const double cell_size,
                        int &gpu_mode, FILE *screen, double host_cut_coulsq,
                        double *host_special_coul, const double qqrd2e, const double g_ewald,
                        const double lamda);
extern void csl_gpu_reinit(const int ntypes, double **scale);
extern void csl_gpu_clear();
extern int **csl_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                               int *host_type, double *sublo, double *subhi, tagint *tag,
                               int **nspecial, tagint **special, const bool eflag, const bool vflag,
                               const bool eatom, const bool vatom, int **ilist, int **jnum,
                               bool &success, double *host_q, double *boxlo, double *prd,
                               int *periodicity);
extern void csl_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                            int *host_type, int *ilist, int *numj, int **firstneigh,
                            const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                            bool &success, double *host_q, const int nlocal, double *boxlo,
                            double *prd);
extern double csl_gpu_bytes();

// pair dpd/coul/slater/long
extern int dpd_coul_slater_long_gpu_init(
    const int ntypes, double **cutsq, double **host_a0, double **host_gamma, double **host_sigma,
    double **host_cut_dpd, double **host_cut_dpdsq, double **host_cut_slatersq, double *special_lj,
    const int inum, const int nall, const int max_nbors, const int maxspecial,
    const double cell_size, int &gpu_mode, FILE *screen, double *host_special_coul,
    const double qqrd2e, const double g_ewald, const double lamda);
extern void dpd_coul_slater_long_gpu_clear();
extern int **dpd_coul_slater_long_gpu_compute_n(
    const int ago, const int inum_full, const int nall, double **host_x, int *host_type,
    double *sublo, double *subhi, tagint *tag, int **nspecial, tagint **special, const bool eflag,
    const bool vflag, const bool eatom, const bool vatom, int **ilist, int **jnum, bool &success,
    double **host_v, const double dtinvsqrt, const int seed, const int timestep, double *boxlo,
    double *prd);
extern void dpd_coul_slater_long_gpu_compute(const int ago, const int inum_full, const int nall,
                                             double **host_x, int *host_type, int *ilist, int *numj,
                                             int **firstneigh, const bool eflag, const bool vflag,
                                             const bool eatom, const bool vatom, bool &success,
                                             tagint *tag, double **host_v, const double dtinvsqrt,
                                             const int seed, const int timestep, const int nlocal,
                                             double *boxlo, double *prd);
extern void dpd_coul_slater_long_gpu_get_extra_data(double *host_q);
extern double dpd_coul_slater_long_gpu_bytes();

// pair dpd
extern int dpd_gpu_init(const int ntypes, double **cutsq, double **host_a0, double **host_gamma,
                        double **host_sigma, double **host_cut, double *special_lj, const int inum,
                        const int nall, const int max_nbors, const int maxspecial,
                        const double cell_size, int &gpu_mode, FILE *screen);
extern void dpd_gpu_clear();
extern int **dpd_gpu_compute_n(const int ago, const int inum_full, const int nall, double **host_x,
                               int *host_type, double *sublo, double *subhi, tagint *tag,
                               int **nspecial, tagint **special, const bool eflag, const bool vflag,
                               const bool eatom, const bool vatom, int **ilist, int **jnum,
                               bool &success, double **host_v, const double dtinvsqrt,
                               const int seed, const int timestep, double *boxlo, double *prd);
extern void dpd_gpu_compute(const int ago, const int inum_full, const int nall, double **host_x,
                            int *host_type, int *ilist, int *numj, int **firstneigh,
                            const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                            bool &success, tagint *tag, double **host_v, const double dtinvsqrt,
                            const int seed, const int timestep, const int nlocal, double *boxlo,
                            double *prd);
extern double dpd_gpu_bytes();

// pair dpd/tstat
extern int dpd_tstat_gpu_init(const int ntypes, double **cutsq, double **host_a0,
                              double **host_gamma, double **host_sigma, double **host_cut,
                              double *special_lj, const int inum, const int nall,
                              const int max_nbors, const int maxspecial, const double cell_size,
                              int &gpu_mode, FILE *screen);
extern void dpd_tstat_gpu_clear();
extern int **dpd_tstat_gpu_compute_n(const int ago, const int inum_full, const int nall,
                                     double **host_x, int *host_type, double *sublo, double *subhi,
                                     tagint *tag, int **nspecial, tagint **special,
                                     const bool eflag, const bool vflag, const bool eatom,
                                     const bool vatom, int **ilist, int **jnum, bool &success,
                                     double **host_v, const double dtinvsqrt, const int seed,
                                     const int timestep, double *boxlo, double *prd);
extern void dpd_tstat_gpu_compute(const int ago, const int inum_full, const int nall,
                                  double **host_x, int *host_type, int *ilist, int *numj,
                                  int **firstneigh, const bool eflag, const bool vflag,
                                  const bool eatom, const bool vatom, bool &success, tagint *tag,
                                  double **host_v, const double dtinvsqrt, const int seed,
                                  const int timestep, const int nlocal, double *boxlo, double *prd);
extern void dpd_tstat_gpu_update_coeff(int ntypes, double **host_a0, double **host_gamma,
                                       double **host_sigma, double **host_cut);
extern double dpd_tstat_gpu_bytes();

// pair eam/alloy
extern int eam_alloy_gpu_init(const int ntypes, double host_cutforcesq, int **host_type2rhor,
                              int **host_type2z2r, int *host_type2frho, double ***host_rhor_spline,
                              double ***host_z2r_spline, double ***host_frho_spline,
                              double **host_cutsq, double rdr, double rdrho, double rhomax,
                              double rhomin, const int he_flag, int nrhor, int nrho, int nz2r,
                              int nfrho, int nr, const int nlocal, const int nall,
                              const int max_nbors, const int maxspecial, const double cell_size,
                              int &gpu_mode, FILE *screen, int &fp_size);
extern void eam_alloy_gpu_clear();
extern int **eam_alloy_gpu_compute_n(const int ago, const int inum_full, const int nall,
                                     double **host_x, int *host_type, double *sublo, double *subhi,
                                     tagint *tag, int **nspecial, tagint **special,
                                     const bool eflag, const bool vflag, const bool eatom,
                                     const bool vatom, int **ilist, int **jnum, bool &success,
                                     int &inum, void **fp_ptr, double *prd, int *periodicity);
extern void eam_alloy_gpu_compute(const int ago, const int inum_full, const int nlocal,
                                  const int nall, double **host_x, int *host_type, int *ilist,
                                  int *numj, int **firstneigh, const bool eflag, const bool vflag,
                                  const bool eatom, const bool vatom, bool &success, void **fp_ptr);
extern void eam_alloy_gpu_compute_force(int *ilist, const bool eflag, const bool vflag,
                                        const bool eatom, const bool vatom);
extern double eam_alloy_gpu_bytes();

// pair eam/fs
extern int eam_fs_gpu_init(const int ntypes, double host_cutforcesq, int **host_type2rhor,
                           int **host_type2z2r, int *host_type2frho, double ***host_rhor_spline,
                           double ***host_z2r_spline, double ***host_frho_spline,
                           double **host_cutsq, double rdr, double rdrho, double rhomax,
                           double rhomin, const int he_flag, int nrhor, int nrho, int nz2r,
                           int nfrho, int nr, const int nlocal, const int nall, const int max_nbors,
                           const int maxspecial, const double cell_size, int &gpu_mode,
                           FILE *screen, int &fp_size);
extern void eam_fs_gpu_clear();
extern int **eam_fs_gpu_compute_n(const int ago, const int inum_full, const int nall,
                                  double **host_x, int *host_type, double *sublo, double *subhi,
                                  tagint *tag, int **nspecial, tagint **special, const bool eflag,
                                  const bool vflag, const bool eatom, const bool vatom, int **ilist,
                                  int **jnum, bool &success, int &inum, void **fp_ptr, double *prd,
                                  int *periodicity);
extern void eam_fs_gpu_compute(const int ago, const int inum_full, const int nlocal, const int nall,
                               double **host_x, int *host_type, int *ilist, int *numj,
                               int **firstneigh, const bool eflag, const bool vflag,
                               const bool eatom, const bool vatom, bool &success, void **fp_ptr);
extern void eam_fs_gpu_compute_force(int *ilist, const bool eflag, const bool vflag,
                                     const bool eatom, const bool vatom);
extern double eam_fs_gpu_bytes();

// pair eam
extern int eam_gpu_init(const int ntypes, double host_cutforcesq, int **host_type2rhor,
                        int **host_type2z2r, int *host_type2frho, double ***host_rhor_spline,
                        double ***host_z2r_spline, double ***host_frho_spline, double **host_cutsq,
                        double rdr, double rdrho, double rhomax, double rhomin, const int he_flag,
                        int nrhor, int nrho, int nz2r, int nfrho, int nr, const int nlocal,
                        const int nall, const int max_nbors, const int maxspecial,
                        const double cell_size, int &gpu_mode, FILE *screen, int &fp_size);
extern void eam_gpu_clear();
extern int **eam_gpu_compute_n(const int ago, const int inum_full, const int nall, double **host_x,
                               int *host_type, double *sublo, double *subhi, tagint *tag,
                               int **nspecial, tagint **special, const bool eflag, const bool vflag,
                               const bool eatom, const bool vatom, int **ilist, int **jnum,
                               bool &success, int &inum, void **fp_ptr, double *prd,
                               int *periodicity);
extern void eam_gpu_compute(const int ago, const int inum_full, const int nlocal, const int nall,
                            double **host_x, int *host_type, int *ilist, int *numj,
                            int **firstneigh, const bool eflag, const bool vflag, const bool eatom,
                            const bool vatom, bool &success, void **fp_ptr);
extern void eam_gpu_compute_force(int *ilist, const bool eflag, const bool vflag, const bool eatom,
                                  const bool vatom);
extern double eam_gpu_bytes();

// pair eam/he
extern int eam_he_gpu_init(const int ntypes, double host_cutforcesq, int **host_type2rhor,
                           int **host_type2z2r, int *host_type2frho, double ***host_rhor_spline,
                           double ***host_z2r_spline, double ***host_frho_spline,
                           double **host_cutsq, double rdr, double rdrho, double rhomax,
                           double rhomin, const int he_flag, int nrhor, int nrho, int nz2r,
                           int nfrho, int nr, const int nlocal, const int nall, const int max_nbors,
                           const int maxspecial, const double cell_size, int &gpu_mode,
                           FILE *screen, int &fp_size);
extern void eam_he_gpu_clear();
extern int **eam_he_gpu_compute_n(const int ago, const int inum_full, const int nall,
                                  double **host_x, int *host_type, double *sublo, double *subhi,
                                  tagint *tag, int **nspecial, tagint **special, const bool eflag,
                                  const bool vflag, const bool eatom, const bool vatom, int **ilist,
                                  int **jnum, bool &success, int &inum, void **fp_ptr, double *prd,
                                  int *periodicity);
extern void eam_he_gpu_compute(const int ago, const int inum_full, const int nlocal, const int nall,
                               double **host_x, int *host_type, int *ilist, int *numj,
                               int **firstneigh, const bool eflag, const bool vflag,
                               const bool eatom, const bool vatom, bool &success, void **fp_ptr);
extern void eam_he_gpu_compute_force(int *ilist, const bool eflag, const bool vflag,
                                     const bool eatom, const bool vatom);
extern double eam_he_gpu_bytes();

// pair edpd
extern int edpd_gpu_init(const int ntypes, double **cutsq, double **host_a0, double **host_gamma,
                         double **host_cut, double **host_power, double **host_kappa,
                         double **host_powerT, double **host_cutT, double ***host_sc,
                         double ***host_kc, double *host_mass, double *special_lj,
                         const int power_flag, const int kappa_flag, const int inum, const int nall,
                         const int max_nbors, const int maxspecial, const double cell_size,
                         int &gpu_mode, FILE *screen);
extern void edpd_gpu_clear();
extern int **edpd_gpu_compute_n(const int ago, const int inum_full, const int nall, double **host_x,
                                int *host_type, double *sublo, double *subhi, tagint *tag,
                                int **nspecial, tagint **special, const bool eflag,
                                const bool vflag, const bool eatom, const bool vatom, int **ilist,
                                int **jnum, bool &success, double **host_v, const double dtinvsqrt,
                                const int seed, const int timestep, double *boxlo, double *prd);
extern void edpd_gpu_compute(const int ago, const int inum_full, const int nall, double **host_x,
                             int *host_type, int *ilist, int *numj, int **firstneigh,
                             const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                             bool &success, tagint *tag, double **host_v, const double dtinvsqrt,
                             const int seed, const int timestep, const int nlocal, double *boxlo,
                             double *prd);
extern void edpd_gpu_get_extra_data(double *host_T, double *host_cv);
extern void edpd_gpu_update_flux(void **flux_ptr);
extern double edpd_gpu_bytes();

// pair gauss
extern int gauss_gpu_init(const int ntypes, double **cutsq, double **host_a, double **b,
                          double **offset, double *special_lj, const int nlocal, const int nall,
                          const int max_nbors, const int maxspecial, const double cell_size,
                          int &gpu_mode, FILE *screen);
extern void gauss_gpu_reinit(const int ntypes, double **cutsq, double **host_a, double **b,
                             double **offset);
extern void gauss_gpu_clear();
extern int **gauss_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                                 int *host_type, double *sublo, double *subhi, tagint *tag,
                                 int **nspecial, tagint **special, const bool eflag,
                                 const bool vflag, const bool eatom, const bool vatom, int **ilist,
                                 int **jnum, bool &success, double *prd, int *periodicity);
extern void gauss_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                              int *host_type, int *ilist, int *numj, int **firstneigh,
                              const bool eflag, const bool vflag, const bool eatom,
                              const bool vatom, bool &success);
extern double gauss_gpu_bytes();

// pair gayberne
extern int gb_gpu_init(const int ntypes, const double gamma, const double upsilon, const double mu,
                       double **shape, double **well, double **cutsq, double **sigma,
                       double **epsilon, double *host_lshape, int **form, double **host_lj1,
                       double **host_lj2, double **host_lj3, double **host_lj4, double **offset,
                       double *special_lj, const int nlocal, const int nall, const int max_nbors,
                       const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen);
extern void gb_gpu_clear();
extern int **gb_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                              int *host_type, double *sublo, double *subhi, tagint *tag,
                              int **nspecial, tagint **special, const bool eflag, const bool vflag,
                              const bool eatom, const bool vatom, int **ilist, int **jnum,
                              bool &success, const int *ellipsoid, const void *bonus);
extern int *gb_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                           int *host_type, int *ilist, int *numj, int **firstneigh,
                           const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                           bool &success, const int *ellipsoid, const void *bonus);
extern double gb_gpu_bytes();

// pair hippo
extern int hippo_gpu_init(const int ntypes, const int max_amtype, const int max_amclass,
                          const double *host_pdamp, const double *host_thole,
                          const double *host_dirdamp, const int *host_amtype2class,
                          const double *host_special_repel, const double *host_special_disp,
                          const double *host_special_mpole, const double *host_special_polar_wscale,
                          const double *host_special_polar_piscale,
                          const double *host_special_polar_pscale, const double *host_sizpr,
                          const double *host_dmppr, const double *host_elepr,
                          const double *host_csix, const double *host_adisp,
                          const double *host_pcore, const double *host_palpha, const int nlocal,
                          const int nall, const int max_nbors, const int maxspecial,
                          const int maxspecial15, const double cell_size, int &gpu_mode,
                          FILE *screen, const double polar_dscale, const double polar_uscale);
extern void hippo_gpu_clear();
extern int **hippo_gpu_precompute(const int ago, const int inum_full, const int nall,
                                  double **host_x, int *host_type, int *host_amtype,
                                  int *host_amgroup, double **host_rpole, double **host_uind,
                                  double **host_uinp, double *host_pval, double *sublo,
                                  double *subhi, tagint *tag, int **nspecial, tagint **special,
                                  int *nspecial15, tagint **special15, const bool eflag_in,
                                  const bool vflag_in, const bool eatom, const bool vatom,
                                  int **ilist, int **jnum, bool &success, double *host_q,
                                  double *boxlo, double *prd);
extern void hippo_gpu_compute_repulsion(
    const int ago, const int inum_full, const int nall, double **host_x, int *host_type,
    int *host_amtype, int *host_amgroup, double **host_rpole, double *sublo, double *subhi,
    tagint *tag, int **nspecial, tagint **special, int *nspecial15, tagint **special15,
    const bool eflag, const bool vflag, const bool eatom, const bool vatom, int **ilist, int **jnum,
    bool &success, const double aewald, const double off2, double *host_q, double *boxlo,
    double *prd, double cut2, double c0, double c1, double c2, double c3, double c4, double c5,
    void **tep_ptr);
extern void hippo_gpu_compute_dispersion_real(int *host_amtype, int *host_amgroup,
                                              double **host_rpole, const double aewald,
                                              const double off2);
extern void hippo_gpu_compute_multipole_real(
    const int ago, const int inum, const int nall, double **host_x, int *host_type,
    int *host_amtype, int *host_amgroup, double **host_rpole, double *host_pval, double *sublo,
    double *subhi, tagint *tag, int **nspecial, tagint **special, int *nspecial15,
    tagint **special15, const bool eflag, const bool vflag, const bool eatom, const bool vatom,
    int **ilist, int **jnum, bool &success, const double aewald, const double felec,
    const double off2, double *host_q, double *boxlo, double *prd, void **tq_ptr);
extern void hippo_gpu_compute_udirect2b(int *host_amtype, int *host_amgroup, double **host_rpole,
                                        double **host_uind, double **host_uinp, double *host_pval,
                                        const double aewald, const double off2, void **fieldp_ptr);
extern void hippo_gpu_compute_umutual2b(int *host_amtype, int *host_amgroup, double **host_rpole,
                                        double **host_uind, double **host_uinp, double *host_pval,
                                        const double aewald, const double off2, void **fieldp_ptr);
extern void hippo_gpu_update_fieldp(void **fieldp_ptr);
extern void hippo_gpu_precompute_kspace(const int inum_full, const int bsorder,
                                        double ***host_thetai1, double ***host_thetai2,
                                        double ***host_thetai3, int **igrid, const int nzlo_out,
                                        const int nzhi_out, const int nylo_out, const int nyhi_out,
                                        const int nxlo_out, const int nxhi_out);
extern void hippo_gpu_fphi_uind(double ****host_grid_brick, void **host_fdip_phi1,
                                void **host_fdip_phi2, void **host_fdip_sum_phi);
extern void hippo_gpu_compute_polar_real(int *host_amtype, int *host_amgroup, double **host_rpole,
                                         double **host_uind, double **host_uinp, double *host_pval,
                                         const bool eflag, const bool vflag, const bool eatom,
                                         const bool vatom, const double aewald, const double felec,
                                         const double off2, void **tq_ptr);
extern double hippo_gpu_bytes();

// pair lj/96/cut, pair lj/class2
extern int lj96_gpu_init(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                         double **host_lj3, double **host_lj4, double **offset, double *special_lj,
                         const int nlocal, const int nall, const int max_nbors,
                         const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen);
extern void lj96_gpu_clear();
extern int **lj96_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                                int *host_type, double *sublo, double *subhi, tagint *tag,
                                int **nspecial, tagint **special, const bool eflag,
                                const bool vflag, const bool eatom, const bool vatom, int **ilist,
                                int **jnum, bool &success, double *prd, int *periodicity);
extern void lj96_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                             int *host_type, int *ilist, int *numj, int **firstneigh,
                             const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                             bool &success);
extern double lj96_gpu_bytes();

// pair lj/charmm/coul/charmm
extern int crm_gpu_init(const int ntypes, double cut_bothsq, double **host_lj1, double **host_lj2,
                        double **host_lj3, double **host_lj4, double *special_lj, const int nlocal,
                        const int nall, const int max_nbors, const int maxspecial,
                        const double cell_size, int &gpu_mode, FILE *screen, double host_cut_ljsq,
                        double host_cut_coulsq, double *host_special_coul, const double qqrd2e,
                        const double cut_lj_innersq, const double cut_coul_innersq,
                        const double denom_lj, const double denom_coul, double **epsilon,
                        double **sigma, const bool mix_arithmetic);
extern void crm_gpu_clear();
extern int **crm_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                               int *host_type, double *sublo, double *subhi, tagint *tag,
                               int **nspecial, tagint **special, const bool eflag, const bool vflag,
                               const bool eatom, const bool vatom, int **ilist, int **jnum,
                               bool &success, double *host_q, double *boxlo, double *prd,
                               int *periodicity);
extern void crm_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                            int *host_type, int *ilist, int *numj, int **firstneigh,
                            const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                            bool &success, double *host_q, const int nlocal, double *boxlo,
                            double *prd);
extern double crm_gpu_bytes();

// pair lj/charmm/coul/long
extern int crml_gpu_init(const int ntypes, double cut_bothsq, double **host_lj1, double **host_lj2,
                         double **host_lj3, double **host_lj4, double **offset, double *special_lj,
                         const int nlocal, const int nall, const int max_nbors,
                         const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen,
                         double host_cut_ljsq, double host_cut_coulsq, double *host_special_coul,
                         const double qqrd2e, const double g_ewald, const double cut_lj_innersq,
                         const double denom_lj, double **epsilon, double **sigma,
                         const bool mix_arithmetic);
extern void crml_gpu_clear();
extern int **crml_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                                int *host_type, double *sublo, double *subhi, tagint *tag,
                                int **nspecial, tagint **special, const bool eflag,
                                const bool vflag, const bool eatom, const bool vatom, int **ilist,
                                int **jnum, bool &success, double *host_q, double *boxlo,
                                double *prd, int *periodicity);
extern void crml_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                             int *host_type, int *ilist, int *numj, int **firstneigh,
                             const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                             bool &success, double *host_q, const int nlocal, double *boxlo,
                             double *prd);
extern double crml_gpu_bytes();

// pair lj/class2/coul/long
extern int c2cl_gpu_init(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                         double **host_lj3, double **host_lj4, double **offset, double *special_lj,
                         const int nlocal, const int nall, const int max_nbors,
                         const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen,
                         double **host_cut_ljsq, double host_cut_coulsq, double *host_special_coul,
                         const double qqrd2e, const double g_ewald);
extern void c2cl_gpu_clear();
extern int **c2cl_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                                int *host_type, double *sublo, double *subhi, tagint *tag,
                                int **nspecial, tagint **special, const bool eflag,
                                const bool vflag, const bool eatom, const bool vatom, int **ilist,
                                int **jnum, bool &success, double *host_q, double *boxlo,
                                double *prd, int *periodicity);
extern void c2cl_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                             int *host_type, int *ilist, int *numj, int **firstneigh,
                             const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                             bool &success, double *host_q, const int nlocal, double *boxlo,
                             double *prd);
extern double c2cl_gpu_bytes();

// pair lj/cubic
extern int ljcb_gpu_init(const int ntypes, double **cutsq, double **cut_inner_sq,
                         double **cut_inner, double **sigma, double **epsilon, double **host_lj1,
                         double **host_lj2, double **host_lj3, double **host_lj4,
                         double *special_lj, const int nlocal, const int nall, const int max_nbors,
                         const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen);
extern void ljcb_gpu_clear();
extern int **ljcb_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                                int *host_type, double *sublo, double *subhi, tagint *tag,
                                int **nspecial, tagint **special, const bool eflag,
                                const bool vflag, const bool eatom, const bool vatom, int **ilist,
                                int **jnum, bool &success, double *prd, int *periodicity);
extern void ljcb_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                             int *host_type, int *ilist, int *numj, int **firstneigh,
                             const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                             bool &success);
extern double ljcb_gpu_bytes();

// pair lj/cut/coul/cut
extern int ljc_gpu_init(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                        double **host_lj3, double **host_lj4, double **offset, double *special_lj,
                        const int nlocal, const int nall, const int max_nbors, const int maxspecial,
                        const double cell_size, int &gpu_mode, FILE *screen, double **host_cut_ljsq,
                        double **host_cut_coulsq, double *host_special_coul, const double qqrd2e);
extern void ljc_gpu_clear();
extern int **ljc_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                               int *host_type, double *sublo, double *subhi, tagint *tag,
                               int **nspecial, tagint **special, const bool eflag, const bool vflag,
                               const bool eatom, const bool vatom, int **ilist, int **jnum,
                               bool &success, double *host_q, double *boxlo, double *prd,
                               int *periodicity);
extern void ljc_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                            int *host_type, int *ilist, int *numj, int **firstneigh,
                            const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                            bool &success, double *host_q, const int nlocal, double *boxlo,
                            double *prd);
extern double ljc_gpu_bytes();

// pair lj/cut/coul/cut/soft
extern int ljcs_gpu_init(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                         double **host_lj3, double **host_lj4, double **offset, double **epsilon,
                         double *special_lj, const int nlocal, const int nall, const int max_nbors,
                         const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen,
                         double **host_cut_ljsq, double **host_cut_coulsq,
                         double *host_special_coul, const double qqrd2e);
extern void ljcs_gpu_clear();
extern int **ljcs_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                                int *host_type, double *sublo, double *subhi, tagint *tag,
                                int **nspecial, tagint **special, const bool eflag,
                                const bool vflag, const bool eatom, const bool vatom, int **ilist,
                                int **jnum, bool &success, double *host_q, double *boxlo,
                                double *prd, int *periodicity);
extern void ljcs_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                             int *host_type, int *ilist, int *numj, int **firstneigh,
                             const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                             bool &success, double *host_q, const int nlocal, double *boxlo,
                             double *prd);
extern double ljcs_gpu_bytes();

// pair lj/cut/coul/debye
extern int ljcd_gpu_init(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                         double **host_lj3, double **host_lj4, double **offset, double *special_lj,
                         const int nlocal, const int nall, const int max_nbors,
                         const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen,
                         double **host_cut_ljsq, double **host_cut_coulsq,
                         double *host_special_coul, const double qqrd2e, const double kappa);
extern void ljcd_gpu_clear();
extern int **ljcd_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                                int *host_type, double *sublo, double *subhi, tagint *tag,
                                int **nspecial, tagint **special, const bool eflag,
                                const bool vflag, const bool eatom, const bool vatom, int **ilist,
                                int **jnum, bool &success, double *host_q, double *boxlo,
                                double *prd, int *periodicity);
extern void ljcd_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                             int *host_type, int *ilist, int *numj, int **firstneigh,
                             const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                             bool &success, double *host_q, const int nlocal, double *boxlo,
                             double *prd);
extern double ljcd_gpu_bytes();

// pair lj/cut/coul/dsf
extern int ljd_gpu_init(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                        double **host_lj3, double **host_lj4, double **offset, double *special_lj,
                        const int nlocal, const int nall, const int max_nbors, const int maxspecial,
                        const double cell_size, int &gpu_mode, FILE *screen, double **host_cut_ljsq,
                        const double host_cut_coulsq, double *host_special_coul,
                        const double qqrd2e, const double e_shift, const double f_shift,
                        const double alpha);
extern void ljd_gpu_clear();
extern int **ljd_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                               int *host_type, double *sublo, double *subhi, tagint *tag,
                               int **nspecial, tagint **special, const bool eflag, const bool vflag,
                               const bool eatom, const bool vatom, int **ilist, int **jnum,
                               bool &success, double *host_q, double *boxlo, double *prd,
                               int *periodicity);
extern void ljd_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                            int *host_type, int *ilist, int *numj, int **firstneigh,
                            const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                            bool &success, double *host_q, const int nlocal, double *boxlo,
                            double *prd);
extern double ljd_gpu_bytes();

// pair lj/cut/coul/long
extern int ljcl_gpu_init(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                         double **host_lj3, double **host_lj4, double **offset, double *special_lj,
                         const int nlocal, const int nall, const int max_nbors,
                         const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen,
                         double **host_cut_ljsq, double host_cut_coulsq, double *host_special_coul,
                         const double qqrd2e, const double g_ewald);
extern void ljcl_gpu_reinit(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                            double **host_lj3, double **host_lj4, double **offset,
                            double **host_lj_cutsq);
extern void ljcl_gpu_clear();
extern int **ljcl_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                                int *host_type, double *sublo, double *subhi, tagint *tag,
                                int **nspecial, tagint **special, const bool eflag,
                                const bool vflag, const bool eatom, const bool vatom, int **ilist,
                                int **jnum, bool &success, double *host_q, double *boxlo,
                                double *prd, int *periodicity);
extern void ljcl_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                             int *host_type, int *ilist, int *numj, int **firstneigh,
                             const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                             bool &success, double *host_q, const int nlocal, double *boxlo,
                             double *prd);
extern double ljcl_gpu_bytes();

// pair lj/cut/coul/long/soft
extern int ljcls_gpu_init(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                          double **host_lj3, double **host_lj4, double **offset, double **epsilon,
                          double *special_lj, const int nlocal, const int nall, const int max_nbors,
                          const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen,
                          double **host_cut_ljsq, double host_cut_coulsq, double *host_special_coul,
                          const double qqrd2e, const double g_ewald);
extern void ljcls_gpu_reinit(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                             double **host_lj3, double **host_lj4, double **offset,
                             double **epsilon, double **host_lj_cutsq);
extern void ljcls_gpu_clear();
extern int **ljcls_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                                 int *host_type, double *sublo, double *subhi, tagint *tag,
                                 int **nspecial, tagint **special, const bool eflag,
                                 const bool vflag, const bool eatom, const bool vatom, int **ilist,
                                 int **jnum, bool &success, double *host_q, double *boxlo,
                                 double *prd, int *periodicity);
extern void ljcls_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                              int *host_type, int *ilist, int *numj, int **firstneigh,
                              const bool eflag, const bool vflag, const bool eatom,
                              const bool vatom, bool &success, double *host_q, const int nlocal,
                              double *boxlo, double *prd);
extern double ljcls_gpu_bytes();

// pair lj/cut/coul/msm
extern int ljcm_gpu_init(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                         double **host_lj3, double **host_lj4, double **host_gcons,
                         double **host_dgcons, double **offset, double *special_lj, const int inum,
                         const int nall, const int max_nbors, const int maxspecial,
                         const double cell_size, int &gpu_mode, FILE *screen,
                         double **host_cut_ljsq, double host_cut_coulsq, double *host_special_coul,
                         const int order, const double qqrd2e);
extern void ljcm_gpu_clear();
extern int **ljcm_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                                int *host_type, double *sublo, double *subhi, tagint *tag,
                                int **nspecial, tagint **special, const bool eflag,
                                const bool vflag, const bool eatom, const bool vatom, int **ilist,
                                int **jnum, bool &success, double *host_q, double *boxlo,
                                double *prd, int *periodicity);
extern void ljcm_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                             int *host_type, int *ilist, int *numj, int **firstneigh,
                             const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                             bool &success, double *host_q, const int nlocal, double *boxlo,
                             double *prd);
extern double ljcm_gpu_bytes();

// pair lj/cut/dipole/cut
extern int dpl_gpu_init(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                        double **host_lj3, double **host_lj4, double **offset, double *special_lj,
                        const int nlocal, const int nall, const int max_nbors, const int maxspecial,
                        const double cell_size, int &gpu_mode, FILE *screen, double **host_cut_ljsq,
                        double **host_cut_coulsq, double *host_special_coul, const double qqrd2e);
extern void dpl_gpu_clear();
extern int **dpl_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                               int *host_type, double *sublo, double *subhi, tagint *tag,
                               int **nspecial, tagint **special, const bool eflag, const bool vflag,
                               const bool eatom, const bool vatom, int **ilist, int **jnum,
                               bool &success, double *host_q, double **host_mu, double *boxlo,
                               double *prd);
extern void dpl_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                            int *host_type, int *ilist, int *numj, int **firstneigh,
                            const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                            bool &success, double *host_q, double **host_mu, const int nlocal,
                            double *boxlo, double *prd);
extern double dpl_gpu_bytes();

// pair lj/cut/dipole/long
extern int dplj_gpu_init(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                         double **host_lj3, double **host_lj4, double **offset, double *special_lj,
                         const int nlocal, const int nall, const int max_nbors,
                         const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen,
                         double **host_cut_ljsq, const double host_cut_coulsq,
                         double *host_special_coul, const double qqrd2e, const double g_ewald);
extern void dplj_gpu_clear();
extern int **dplj_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                                int *host_type, double *sublo, double *subhi, tagint *tag,
                                int **nspecial, tagint **special, const bool eflag,
                                const bool vflag, const bool eatom, const bool vatom, int **ilist,
                                int **jnum, bool &success, double *host_q, double **host_mu,
                                double *boxlo, double *prd);
extern void dplj_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                             int *host_type, int *ilist, int *numj, int **firstneigh,
                             const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                             bool &success, double *host_q, double **host_mu, const int nlocal,
                             double *boxlo, double *prd);
extern double dplj_gpu_bytes();

// pair lj/cut
extern int ljl_gpu_init(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                        double **host_lj3, double **host_lj4, double **offset, double *special_lj,
                        const int nlocal, const int nall, const int max_nbors, const int maxspecial,
                        const double cell_size, int &gpu_mode, FILE *screen);
extern void ljl_gpu_reinit(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                           double **host_lj3, double **host_lj4, double **offset);
extern void ljl_gpu_clear();
extern int **ljl_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                               int *host_type, double *sublo, double *subhi, tagint *tag,
                               int **nspecial, tagint **special, const bool eflag, const bool vflag,
                               const bool eatom, const bool vatom, int **ilist, int **jnum,
                               bool &success, double *prd, int *periodicity);
extern void ljl_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                            int *host_type, int *ilist, int *numj, int **firstneigh,
                            const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                            bool &success);
extern double ljl_gpu_bytes();

// pair lj/cut/tip4p/long
extern int ljtip4p_long_gpu_init(const int ntypes, double **cutsq, double **host_lj1,
                                 double **host_lj2, double **host_lj3, double **host_lj4,
                                 double **offset, double *special_lj, const int nlocal,
                                 const int tH, const int tO, const double alpha, const double qdist,
                                 const int nall, const int max_nbors, const int maxspecial,
                                 const double cell_size, int &gpu_mode, FILE *screen,
                                 double **host_cut_ljsq, const double host_cut_coulsq,
                                 const double host_cut_coulsqplus, double *host_special_coul,
                                 const double qqrd2e, const double g_ewald, int map_size,
                                 int max_same);
extern void ljtip4p_long_gpu_clear();
extern int **ljtip4p_long_gpu_compute_n(const int ago, const int inum, const int nall,
                                        double **host_x, int *host_type, double *sublo,
                                        double *subhi, tagint *tag, int *map_array, int map_size,
                                        int *sametag, int max_same, int **nspecial,
                                        tagint **special, const bool eflag, const bool vflag,
                                        const bool eatom, const bool vatom, int **ilist, int **jnum,
                                        bool &success, double *host_q, double *boxlo, double *prd,
                                        int *periodicity);
extern void ljtip4p_long_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                                     int *host_type, int *ilist, int *numj, int **firstneigh,
                                     const bool eflag, const bool vflag, const bool eatom,
                                     const bool vatom, bool &success, double *host_q,
                                     const int nlocal, double *boxlo, double *prd);
extern double ljtip4p_long_gpu_bytes();
extern void ljtip4p_long_copy_molecule_data(int, tagint *, int *, int, int *, int, int);

// pair lj/expand/coul/long
extern int ljecl_gpu_init(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                          double **host_lj3, double **host_lj4, double **offset, double **shift,
                          double *special_lj, const int nlocal, const int nall, const int max_nbors,
                          const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen,
                          double **host_cut_ljsq, double host_cut_coulsq, double *host_special_coul,
                          const double qqrd2e, const double g_ewald);
extern void ljecl_gpu_reinit(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                             double **host_lj3, double **host_lj4, double **offset, double **shift,
                             double **host_lj_cutsq);
extern void ljecl_gpu_clear();
extern int **ljecl_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                                 int *host_type, double *sublo, double *subhi, tagint *tag,
                                 int **nspecial, tagint **special, const bool eflag,
                                 const bool vflag, const bool eatom, const bool vatom, int **ilist,
                                 int **jnum, bool &success, double *host_q, double *boxlo,
                                 double *prd, int *periodicity);
extern void ljecl_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                              int *host_type, int *ilist, int *numj, int **firstneigh,
                              const bool eflag, const bool vflag, const bool eatom,
                              const bool vatom, bool &success, double *host_q, const int nlocal,
                              double *boxlo, double *prd);
extern double ljecl_gpu_bytes();

// pair lj/expand
extern int lje_gpu_init(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                        double **host_lj3, double **host_lj4, double **offset, double **shift,
                        double *special_lj, const int nlocal, const int nall, const int max_nbors,
                        const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen);
extern void lje_gpu_reinit(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                           double **host_lj3, double **host_lj4, double **offset, double **shift);
extern void lje_gpu_clear();
extern int **lje_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                               int *host_type, double *sublo, double *subhi, tagint *tag,
                               int **nspecial, tagint **special, const bool eflag, const bool vflag,
                               const bool eatom, const bool vatom, int **ilist, int **jnum,
                               bool &success, double *prd, int *periodicity);
extern void lje_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                            int *host_type, int *ilist, int *numj, int **firstneigh,
                            const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                            bool &success);
extern double lje_gpu_bytes();

// pair lj/gromacs
extern int ljgrm_gpu_init(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                          double **host_lj3, double **host_lj4, double *special_lj, const int inum,
                          const int nall, const int max_nbors, const int maxspecial,
                          const double cell_size, int &gpu_mode, FILE *screen, double **host_ljsw1,
                          double **host_ljsw2, double **host_ljsw3, double **host_ljsw4,
                          double **host_ljsw5, double **cut_inner, double **cut_innersq);
extern void ljgrm_gpu_clear();
extern int **ljgrm_gpu_compute_n(const int ago, const int inum_full, const int nall,
                                 double **host_x, int *host_type, double *sublo, double *subhi,
                                 tagint *tag, int **nspecial, tagint **special, const bool eflag,
                                 const bool vflag, const bool eatom, const bool vatom, int **ilist,
                                 int **jnum, bool &success, double *prd, int *periodicity);
extern void ljgrm_gpu_compute(const int ago, const int inum_full, const int nall, double **host_x,
                              int *host_type, int *ilist, int *numj, int **firstneigh,
                              const bool eflag, const bool vflag, const bool eatom,
                              const bool vatom, bool &success);
extern double ljgrm_gpu_bytes();

// pair lj/sf/dipole/sf
extern int dplsf_gpu_init(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                          double **host_lj3, double **host_lj4, double *special_lj,
                          const int nlocal, const int nall, const int max_nbors,
                          const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen,
                          double **host_cut_ljsq, double **host_cut_coulsq,
                          double *host_special_coul, const double qqrd2e);
extern void dplsf_gpu_clear();
extern int **dplsf_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                                 int *host_type, double *sublo, double *subhi, tagint *tag,
                                 int **nspecial, tagint **special, const bool eflag,
                                 const bool vflag, const bool eatom, const bool vatom, int **ilist,
                                 int **jnum, bool &success, double *host_q, double **host_mu,
                                 double *boxlo, double *prd);
extern void dplsf_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                              int *host_type, int *ilist, int *numj, int **firstneigh,
                              const bool eflag, const bool vflag, const bool eatom,
                              const bool vatom, bool &success, double *host_q, double **host_mu,
                              const int nlocal, double *boxlo, double *prd);
extern double dplsf_gpu_bytes();

// pair lj/smooth
extern int ljsmt_gpu_init(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                          double **host_lj3, double **host_lj4, double **offset, double *special_lj,
                          const int nlocal, const int nall, const int max_nbors,
                          const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen,
                          double **host_ljsw0, double **host_ljsw1, double **host_ljsw2,
                          double **host_ljsw3, double **host_ljsw4, double **cut_inner,
                          double **cut_innersq);
extern void ljsmt_gpu_reinit(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                             double **host_lj3, double **host_lj4, double **offset,
                             double **host_ljsw0, double **host_ljsw1, double **host_ljsw2,
                             double **host_ljsw3, double **host_ljsw4, double **cut_inner,
                             double **cut_innersq);
extern void ljsmt_gpu_clear();
extern int **ljsmt_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                                 int *host_type, double *sublo, double *subhi, tagint *tag,
                                 int **nspecial, tagint **special, const bool eflag,
                                 const bool vflag, const bool eatom, const bool vatom, int **ilist,
                                 int **jnum, bool &success, double *prd, int *periodicity);
extern void ljsmt_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                              int *host_type, int *ilist, int *numj, int **firstneigh,
                              const bool eflag, const bool vflag, const bool eatom,
                              const bool vatom, bool &success);
extern double ljsmt_gpu_bytes();

// pair lj/spica/coul/long
extern int spical_gpu_init(const int ntypes, double **cutsq, int **lj_type, double **host_lj1,
                           double **host_lj2, double **host_lj3, double **host_lj4, double **offset,
                           double *special_lj, const int nlocal, const int nall,
                           const int max_nbors, const int maxspecial, const double cell_size,
                           int &gpu_mode, FILE *screen, double **host_cut_ljsq,
                           double host_cut_coulsq, double *host_special_coul, const double qqrd2e,
                           const double g_ewald);
extern void spical_gpu_clear();
extern int **spical_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                                  int *host_type, double *sublo, double *subhi, tagint *tag,
                                  int **nspecial, tagint **special, const bool eflag,
                                  const bool vflag, const bool eatom, const bool vatom, int **ilist,
                                  int **jnum, bool &success, double *host_q, double *boxlo,
                                  double *prd, int *periodicity);
extern void spical_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                               int *host_type, int *ilist, int *numj, int **firstneigh,
                               const bool eflag, const bool vflag, const bool eatom,
                               const bool vatom, bool &success, double *host_q, const int nlocal,
                               double *boxlo, double *prd);
extern double spical_gpu_bytes();

// pair lj/spica
extern int spica_gpu_init(const int ntypes, double **cutsq, int **cg_types, double **host_lj1,
                          double **host_lj2, double **host_lj3, double **host_lj4, double **offset,
                          double *special_lj, const int nlocal, const int nall, const int max_nbors,
                          const int maxspecial, const double cell_size, int &gpu_mode,
                          FILE *screen);
extern void spica_gpu_clear();
extern int **spica_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                                 int *host_type, double *sublo, double *subhi, tagint *tag,
                                 int **nspecial, tagint **special, const bool eflag,
                                 const bool vflag, const bool eatom, const bool vatom, int **ilist,
                                 int **jnum, bool &success, double *prd, int *periodicity);
extern void spica_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                              int *host_type, int *ilist, int *numj, int **firstneigh,
                              const bool eflag, const bool vflag, const bool eatom,
                              const bool vatom, bool &success);
extern double spica_gpu_bytes();

// pair mdpd
extern int mdpd_gpu_init(const int ntypes, double **cutsq, double **host_A_att, double **host_B_rep,
                         double **host_gamma, double **host_sigma, double **host_cut,
                         double **host_cut_r, double *special_lj, const int inum, const int nall,
                         const int max_nbors, const int maxspecial, const double cell_size,
                         int &gpu_mode, FILE *screen);
extern void mdpd_gpu_clear();
extern int **mdpd_gpu_compute_n(const int ago, const int inum_full, const int nall, double **host_x,
                                int *host_type, double *sublo, double *subhi, tagint *tag,
                                int **nspecial, tagint **special, const bool eflag,
                                const bool vflag, const bool eatom, const bool vatom, int **ilist,
                                int **jnum, bool &success, double **host_v, const double dtinvsqrt,
                                const int seed, const int timestep, double *boxlo, double *prd);
extern void mdpd_gpu_compute(const int ago, const int inum_full, const int nall, double **host_x,
                             int *host_type, int *ilist, int *numj, int **firstneigh,
                             const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                             bool &success, tagint *tag, double **host_v, const double dtinvsqrt,
                             const int seed, const int timestep, const int nlocal, double *boxlo,
                             double *prd);
extern void mdpd_gpu_get_extra_data(double *host_rho);
extern double mdpd_gpu_bytes();

// pair mie/cut
extern int mie_gpu_init(const int ntypes, double **cutsq, double **host_mie1, double **host_mie2,
                        double **host_mie3, double **host_mie4, double **host_gamA,
                        double **host_gamR, double **offset, double *special_lj, const int nlocal,
                        const int nall, const int max_nbors, const int maxspecial,
                        const double cell_size, int &gpu_mode, FILE *screen);
extern void mie_gpu_clear();
extern int **mie_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                               int *host_type, double *sublo, double *subhi, tagint *tag,
                               int **nspecial, tagint **special, const bool eflag, const bool vflag,
                               const bool eatom, const bool vatom, int **ilist, int **jnum,
                               bool &success, double *prd, int *periodicity);
extern void mie_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                            int *host_type, int *ilist, int *numj, int **firstneigh,
                            const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                            bool &success);
extern double mie_gpu_bytes();

// pair morse
extern int mor_gpu_init(const int ntypes, double **cutsq, double **host_morse1, double **host_r0,
                        double **host_alpha, double **host_d0, double **offset, double *special_lj,
                        const int nlocal, const int nall, const int max_nbors, const int maxspecial,
                        const double cell_size, int &gpu_mode, FILE *screen);
extern void mor_gpu_clear();
extern int **mor_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                               int *host_type, double *sublo, double *subhi, tagint *tag,
                               int **nspecial, tagint **special, const bool eflag, const bool vflag,
                               const bool eatom, const bool vatom, int **ilist, int **jnum,
                               bool &success, double *prd, int *periodicity);
extern void mor_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                            int *host_type, int *ilist, int *numj, int **firstneigh,
                            const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                            bool &success);
extern double mor_gpu_bytes();

// pair resquared
extern int re_gpu_init(const int ntypes, double **shape, double **well, double **cutsq,
                       double **sigma, double **epsilon, int **form, double **host_lj1,
                       double **host_lj2, double **host_lj3, double **host_lj4, double **offset,
                       double *special_lj, const int nlocal, const int nall, const int max_nbors,
                       const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen);
extern void re_gpu_clear();
extern int **re_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                              int *host_type, double *sublo, double *subhi, tagint *tag,
                              int **nspecial, tagint **special, const bool eflag, const bool vflag,
                              const bool eatom, const bool vatom, int **ilist, int **jnum,
                              bool &success, const int *ellipsoid, const void *bonus);
extern int *re_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                           int *host_type, int *ilist, int *numj, int **firstneigh,
                           const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                           bool &success, const int *ellipsoid, const void *bonus);
extern double re_gpu_bytes();

// pair soft
extern int soft_gpu_init(const int ntypes, double **cutsq, double **prefactor, double **cut,
                         double *special_lj, const int nlocal, const int nall, const int max_nbors,
                         const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen);
extern void soft_gpu_reinit(const int ntypes, double **cutsq, double **host_prefactor,
                            double **host_cut);
extern void soft_gpu_clear();
extern int **soft_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                                int *host_type, double *sublo, double *subhi, tagint *tag,
                                int **nspecial, tagint **special, const bool eflag,
                                const bool vflag, const bool eatom, const bool vatom, int **ilist,
                                int **jnum, bool &success, double *prd, int *periodicity);
extern void soft_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                             int *host_type, int *ilist, int *numj, int **firstneigh,
                             const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                             bool &success);
extern double soft_gpu_bytes();

// pair sph/heatconduction
extern int sph_heatconduction_gpu_init(const int ntypes, double **cutsq, double **host_cut,
                                       double **host_alpha, double *host_mass, const int dimension,
                                       double *special_lj, const int inum, const int nall,
                                       const int max_nbors, const int maxspecial,
                                       const double cell_size, int &gpu_mode, FILE *screen);
extern void sph_heatconduction_gpu_clear();
extern int **sph_heatconduction_gpu_compute_n(const int ago, const int inum_full, const int nall,
                                              double **host_x, int *host_type, double *sublo,
                                              double *subhi, tagint *host_tag, int **nspecial,
                                              tagint **special, const bool eflag, const bool vflag,
                                              const bool eatom, const bool vatom, int **ilist,
                                              int **jnum, bool &success, double **host_v);
extern void sph_heatconduction_gpu_compute(const int ago, const int inum_full, const int nall,
                                           double **host_x, int *host_type, int *ilist, int *numj,
                                           int **firstneigh, const bool eflag, const bool vflag,
                                           const bool eatom, const bool vatom, bool &success,
                                           tagint *host_tag, double **host_v);
extern void sph_heatconduction_gpu_get_extra_data(double *host_rho, double *host_esph);
extern void sph_heatconduction_gpu_update_dE(void **dE_ptr);
extern double sph_heatconduction_gpu_bytes();

// pair sph/lj
extern int sph_lj_gpu_init(const int ntypes, double **cutsq, double **host_cut,
                           double **host_viscosity, double *host_mass, const int dimension,
                           double *special_lj, const int inum, const int nall, const int max_nbors,
                           const int maxspecial, const double cell_size, int &gpu_mode,
                           FILE *screen);
extern void sph_lj_gpu_clear();
extern int **sph_lj_gpu_compute_n(const int ago, const int inum_full, const int nall,
                                  double **host_x, int *host_type, double *sublo, double *subhi,
                                  tagint *host_tag, int **nspecial, tagint **special,
                                  const bool eflag, const bool vflag, const bool eatom,
                                  const bool vatom, int **ilist, int **jnum, bool &success,
                                  double **host_v);
extern void sph_lj_gpu_compute(const int ago, const int inum_full, const int nall, double **host_x,
                               int *host_type, int *ilist, int *numj, int **firstneigh,
                               const bool eflag, const bool vflag, const bool eatom,
                               const bool vatom, bool &success, tagint *host_tag, double **host_v);
extern void sph_lj_gpu_get_extra_data(double *host_rho, double *host_esph, double *host_cv);
extern void sph_lj_gpu_update_drhoE(void **drhoE_ptr);
extern double sph_lj_gpu_bytes();

// pair sph/taitwater
extern int sph_taitwater_gpu_init(const int ntypes, double **cutsq, double **host_cut,
                                  double **host_viscosity, double *host_mass, double *host_rho0,
                                  double *host_soundspeed, double *host_B, const int dimension,
                                  double *special_lj, const int inum, const int nall,
                                  const int max_nbors, const int maxspecial, const double cell_size,
                                  int &gpu_mode, FILE *screen);
extern void sph_taitwater_gpu_clear();
extern int **sph_taitwater_gpu_compute_n(const int ago, const int inum_full, const int nall,
                                         double **host_x, int *host_type, double *sublo,
                                         double *subhi, tagint *tag, int **nspecial,
                                         tagint **special, const bool eflag, const bool vflag,
                                         const bool eatom, const bool vatom, int **ilist,
                                         int **jnum, bool &success, double **host_v);
extern void sph_taitwater_gpu_compute(const int ago, const int inum_full, const int nall,
                                      double **host_x, int *host_type, int *ilist, int *numj,
                                      int **firstneigh, const bool eflag, const bool vflag,
                                      const bool eatom, const bool vatom, bool &success,
                                      tagint *tag, double **host_v);
extern void sph_taitwater_gpu_get_extra_data(double *host_rho);
extern void sph_taitwater_gpu_update_drhoE(void **drhoE_ptr);
extern double sph_taitwater_gpu_bytes();

// pair sw
extern int sw_gpu_init(const int ntypes, const int inum, const int nall, const int max_nbors,
                       const double cell_size, int &gpu_mode, FILE *screen, double **ncutsq,
                       double **ncut, double **sigma, double **powerp, double **powerq,
                       double **sigma_gamma, double **c1, double **c2, double **c3, double **c4,
                       double **c5, double **c6, double ***lambda_epsilon, double ***costheta,
                       const int *map, int ***e2param);
extern void sw_gpu_clear();
extern int **sw_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                              int *host_type, double *sublo, double *subhi, tagint *tag,
                              int **nspecial, tagint **special, const bool eflag, const bool vflag,
                              const bool eatom, const bool vatom, int **ilist, int **jnum,
                              bool &success);
extern void sw_gpu_compute(const int ago, const int nloc, const int nall, const int ln,
                           double **host_x, int *host_type, int *ilist, int *numj, int **firstneigh,
                           const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                           bool &success);
extern double sw_gpu_bytes();

// pair table
extern int table_gpu_init(const int ntypes, double **cutsq, double ***host_table_coeffs,
                          double **host_table_data, double *special_lj, const int nlocal,
                          const int nall, const int max_nbors, const int maxspecial,
                          const double cell_size, int &gpu_mode, FILE *screen, int tabstyle,
                          int ntables, int tablength);
extern void table_gpu_clear();
extern int **table_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                                 int *host_type, double *sublo, double *subhi, tagint *tag,
                                 int **nspecial, tagint **special, const bool eflag,
                                 const bool vflag, const bool eatom, const bool vatom, int **ilist,
                                 int **jnum, bool &success, double *prd, int *periodicity);
extern void table_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                              int *host_type, int *ilist, int *numj, int **firstneigh,
                              const bool eflag, const bool vflag, const bool eatom,
                              const bool vatom, bool &success);
extern double table_gpu_bytes();

// pair tersoff
extern int tersoff_gpu_init(const int ntypes, const int inum, const int nall, const int max_nbors,
                            const double cell_size, int &gpu_mode, FILE *screen, int *host_map,
                            const int nelements, int ***host_elem3param, const int nparams,
                            const double *ts_lam1, const double *ts_lam2, const double *ts_lam3,
                            const double *ts_powermint, const double *ts_biga,
                            const double *ts_bigb, const double *ts_bigr, const double *ts_bigd,
                            const double *ts_c1, const double *ts_c2, const double *ts_c3,
                            const double *ts_c4, const double *ts_c, const double *ts_d,
                            const double *ts_h, const double *ts_gamma, const double *ts_beta,
                            const double *ts_powern, const double *ts_cutsq);
extern void tersoff_gpu_clear();
extern int **tersoff_gpu_compute_n(const int ago, const int inum_full, const int nall,
                                   double **host_x, int *host_type, double *sublo, double *subhi,
                                   tagint *tag, int **nspecial, tagint **special, const bool eflag,
                                   const bool vflag, const bool eatom, const bool vatom,
                                   int **ilist, int **jnum, bool &success);
extern void tersoff_gpu_compute(const int ago, const int nlocal, const int nall, const int nlist,
                                double **host_x, int *host_type, int *ilist, int *numj,
                                int **firstneigh, const bool eflag, const bool vflag,
                                const bool eatom, const bool vatom, bool &success);
extern double tersoff_gpu_bytes();

// pair tersoff/mod
extern int tersoff_mod_gpu_init(
    const int ntypes, const int inum, const int nall, const int max_nbors, const double cell_size,
    int &gpu_mode, FILE *screen, int *host_map, const int nelements, int ***host_elem3param,
    const int nparams, const double *ts_lam1, const double *ts_lam2, const double *ts_lam3,
    const double *ts_powermint, const double *ts_biga, const double *ts_bigb, const double *ts_bigr,
    const double *ts_bigd, const double *ts_c1, const double *ts_c2, const double *ts_c3,
    const double *ts_c4, const double *ts_c5, const double *ts_h, const double *ts_beta,
    const double *ts_powern, const double *ts_powern_del, const double *ts_ca1,
    const double *ts_cutsq);
extern void tersoff_mod_gpu_clear();
extern int **tersoff_mod_gpu_compute_n(const int ago, const int inum_full, const int nall,
                                       double **host_x, int *host_type, double *sublo,
                                       double *subhi, tagint *tag, int **nspecial, tagint **special,
                                       const bool eflag, const bool vflag, const bool eatom,
                                       const bool vatom, int **ilist, int **jnum, bool &success);
extern void tersoff_mod_gpu_compute(const int ago, const int nlocal, const int nall,
                                    const int nlist, double **host_x, int *host_type, int *ilist,
                                    int *numj, int **firstneigh, const bool eflag, const bool vflag,
                                    const bool eatom, const bool vatom, bool &success);
extern double tersoff_mod_gpu_bytes();

// pair tersoff/zbl
extern int tersoff_zbl_gpu_init(
    const int ntypes, const int inum, const int nall, const int max_nbors, const double cell_size,
    int &gpu_mode, FILE *screen, int *host_map, const int nelements, int ***host_elem3param,
    const int nparams, const double *ts_lam1, const double *ts_lam2, const double *ts_lam3,
    const double *ts_powermint, const double *ts_biga, const double *ts_bigb, const double *ts_bigr,
    const double *ts_bigd, const double *ts_c1, const double *ts_c2, const double *ts_c3,
    const double *ts_c4, const double *ts_c, const double *ts_d, const double *ts_h,
    const double *ts_gamma, const double *ts_beta, const double *ts_powern, const double *ts_Z_i,
    const double *ts_Z_j, const double *ts_ZBLcut, const double *ts_ZBLexpscale,
    const double global_e, const double global_a_0, const double global_epsilon_0,
    const double *ts_cutsq);
extern void tersoff_zbl_gpu_clear();
extern int **tersoff_zbl_gpu_compute_n(const int ago, const int inum_full, const int nall,
                                       double **host_x, int *host_type, double *sublo,
                                       double *subhi, tagint *tag, int **nspecial, tagint **special,
                                       const bool eflag, const bool vflag, const bool eatom,
                                       const bool vatom, int **ilist, int **jnum, bool &success);
extern void tersoff_zbl_gpu_compute(const int ago, const int nlocal, const int nall,
                                    const int nlist, double **host_x, int *host_type, int *ilist,
                                    int *numj, int **firstneigh, const bool eflag, const bool vflag,
                                    const bool eatom, const bool vatom, bool &success);
extern double tersoff_zbl_gpu_bytes();

// pair ufm
extern int ufml_gpu_init(const int ntypes, double **cutsq, double **host_uf1, double **host_uf2,
                         double **host_uf3, double **offset, double *special_lj, const int nlocal,
                         const int nall, const int max_nbors, const int maxspecial,
                         const double cell_size, int &gpu_mode, FILE *screen);
extern void ufml_gpu_reinit(const int ntypes, double **cutsq, double **host_uf1, double **host_uf2,
                            double **host_uf3, double **offset);
extern void ufml_gpu_clear();
extern int **ufml_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                                int *host_type, double *sublo, double *subhi, tagint *tag,
                                int **nspecial, tagint **special, const bool eflag,
                                const bool vflag, const bool eatom, const bool vatom, int **ilist,
                                int **jnum, bool &success, double *prd, int *periodicity);
extern void ufml_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                             int *host_type, int *ilist, int *numj, int **firstneigh,
                             const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                             bool &success);
extern double ufml_gpu_bytes();

// pair vashishta
extern int vashishta_gpu_init(const int ntypes, const int inum, const int nall, const int max_nbors,
                              const double cell_size, int &gpu_mode, FILE *screen, int *host_map,
                              const int nelements, int ***host_elem3param, const int nparams,
                              const double *cutsq, const double *r0, const double *gamma,
                              const double *eta, const double *lam1inv, const double *lam4inv,
                              const double *zizj, const double *mbigd, const double *dvrc,
                              const double *big6w, const double *heta, const double *bigh,
                              const double *bigw, const double *c0, const double *costheta,
                              const double *bigb, const double *big2b, const double *bigc);
extern void vashishta_gpu_clear();
extern int **vashishta_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                                     int *host_type, double *sublo, double *subhi, tagint *tag,
                                     int **nspecial, tagint **special, const bool eflag,
                                     const bool vflag, const bool eatom, const bool vatom,
                                     int **ilist, int **jnum, bool &success);
extern void vashishta_gpu_compute(const int ago, const int nloc, const int nall, const int ln,
                                  double **host_x, int *host_type, int *ilist, int *numj,
                                  int **firstneigh, const bool eflag, const bool vflag,
                                  const bool eatom, const bool vatom, bool &success);
extern double vashishta_gpu_bytes();

// pair yukawa/colloid
extern int ykcolloid_gpu_init(const int ntypes, double **cutsq, double **host_a,
                              double **host_offset, double *special_lj, const int inum,
                              const int nall, const int max_nbors, const int maxspecial,
                              const double cell_size, int &gpu_mode, FILE *screen,
                              const double kappa);
extern void ykcolloid_gpu_clear();
extern int **ykcolloid_gpu_compute_n(const int ago, const int inum_full, const int nall,
                                     double **host_x, int *host_type, double *sublo, double *subhi,
                                     tagint *tag, int **nspecial, tagint **special,
                                     const bool eflag, const bool vflag, const bool eatom,
                                     const bool vatom, int **ilist, int **jnum, bool &success,
                                     double *host_rad, double *prd, int *periodicity);
extern void ykcolloid_gpu_compute(const int ago, const int inum_full, const int nall,
                                  double **host_x, int *host_type, int *ilist, int *numj,
                                  int **firstneigh, const bool eflag, const bool vflag,
                                  const bool eatom, const bool vatom, bool &success,
                                  double *host_rad);
extern double ykcolloid_gpu_bytes();

// pair yukawa
extern int yukawa_gpu_init(const int ntypes, double **cutsq, double kappa, double **host_a,
                           double **offset, double *special_lj, const int inum, const int nall,
                           const int max_nbors, const int maxspecial, const double cell_size,
                           int &gpu_mode, FILE *screen);
extern void yukawa_gpu_clear();
extern int **yukawa_gpu_compute_n(const int ago, const int inum_full, const int nall,
                                  double **host_x, int *host_type, double *sublo, double *subhi,
                                  tagint *tag, int **nspecial, tagint **special, const bool eflag,
                                  const bool vflag, const bool eatom, const bool vatom, int **ilist,
                                  int **jnum, bool &success, double *prd, int *periodicity);
extern void yukawa_gpu_compute(const int ago, const int inum_full, const int nall, double **host_x,
                               int *host_type, int *ilist, int *numj, int **firstneigh,
                               const bool eflag, const bool vflag, const bool eatom,
                               const bool vatom, bool &success);
extern double yukawa_gpu_bytes();

// pair zbl
extern int zbl_gpu_init(const int ntypes, double **cutsq, double **host_sw1, double **host_sw2,
                        double **host_sw3, double **host_sw4, double **host_sw5, double **host_d1a,
                        double **host_d2a, double **host_d3a, double **host_d4a, double **host_zze,
                        double cut_globalsq, double cut_innersq, double cut_inner, const int inum,
                        const int nall, const int max_nbors, const int maxspecial,
                        const double cell_size, int &gpu_mode, FILE *screen);
extern void zbl_gpu_clear();
extern int **zbl_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                               int *host_type, double *sublo, double *subhi, tagint *tag,
                               int **nspecial, tagint **special, const bool eflag, const bool vflag,
                               const bool eatom, const bool vatom, int **ilist, int **jnum,
                               bool &success, double *prd, int *periodicity);
extern void zbl_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                            int *host_type, int *ilist, int *numj, int **firstneigh,
                            const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                            bool &success);
extern double zbl_gpu_bytes();

// kspace ewald
extern void ewald_gpu_init(const int nlocal, const int nall, FILE *screen, int &success);
extern void ewald_gpu_setup(const int kmax, const int kcount, int *kxvecs, int *kyvecs,
                            int *kzvecs, double *ug, double **eg, double **vg, double *unitk,
                            int &success);
extern int ewald_gpu_structure(const int ago, const int nlocal, const int nall, double **host_x,
                               int *host_type, double *host_q, double *host_sfacrl,
                               double *host_sfacim, bool &success);
extern void ewald_gpu_compute(double *host_sfacrl_all, double *host_sfacim_all,
                              const double qscale, const int slabflag, const int eflag_atom,
                              const int vflag_atom, double *host_eatom, double **host_vatom,
                              bool &success);
extern void ewald_gpu_clear(const double cpu_time);
extern double ewald_gpu_bytes();

// pppm need two versions for single and double precision FFTs

extern float *pppm_gpu_init_f(const int nlocal, const int nall, FILE *screen, const int order,
                             const int nxlo_out, const int nylo_out, const int nzlo_out,
                             const int nxhi_out, const int nyhi_out, const int nzhi_out,
                             float **rho_coeff, float **_vd_brick, const double slab_volfactor,
                             const int nx_pppm, const int ny_pppm, const int nz_pppm,
                             const bool split, const bool respa, int &success);
extern void pppm_gpu_clear_f(const double poisson_time);
extern int pppm_gpu_spread_f(const int ago, const int nlocal, const int nall, double **host_x,
                             int *host_type, bool &success, double *host_q, double *boxlo,
                             const double delxinv, const double delyinv, const double delzinv);
extern void pppm_gpu_interp_f(const float qqrd2e_scale);
extern double pppm_gpu_bytes_f();
extern void pppm_gpu_forces_f(double **f);

extern double *pppm_gpu_init_d(const int nlocal, const int nall, FILE *screen, const int order,
                               const int nxlo_out, const int nylo_out, const int nzlo_out,
                               const int nxhi_out, const int nyhi_out, const int nzhi_out,
                               double **rho_coeff, double **_vd_brick, const double slab_volfactor,
                               const int nx_pppm, const int ny_pppm, const int nz_pppm,
                               const bool split, const bool respa, int &success);
extern void pppm_gpu_clear_d(const double poisson_time);
extern int pppm_gpu_spread_d(const int ago, const int nlocal, const int nall, double **host_x,
                             int *host_type, bool &success, double *host_q, double *boxlo,
                             const double delxinv, const double delyinv, const double delzinv);
extern void pppm_gpu_interp_d(const double qqrd2e_scale);
extern double pppm_gpu_bytes_d();
extern void pppm_gpu_forces_d(double **f);

}    // namespace LAMMPS_GPU

#endif    // LMP_LAMMPS_GPU_H
