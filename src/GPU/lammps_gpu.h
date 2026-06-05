// clang-format off
/* ----------------------------------------------------------------------
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
   the corresponding lib/gpu/*_ext.cpp or lal_device.cpp.

   Prerequisites: mpi.h (for MPI_Comm) and either lmptype.h (LAMMPS side)
   or lal_precision.h (lib/gpu side) must be included before this header
   for the tagint typedef; <string> must be visible for std::string.        */

#ifndef LMP_LAMMPS_GPU_H
#define LMP_LAMMPS_GPU_H

#include <cstdio>
#include <string>
#include <mpi.h>

/* tagint: re-declaration is legal in C++ when it resolves to the same type.
   On the LAMMPS side lmptype.h already defined it; on the lib/gpu side
   lal_precision.h did.  Guard ensures we don't collide.                    */
#ifndef LAMMPS_GPU_TAGINT_DEFINED
#define LAMMPS_GPU_TAGINT_DEFINED
#ifdef LAMMPS_BIGBIG
#include <cstdint>
typedef int64_t tagint;
#else
typedef int tagint;
#endif
#endif

namespace LAMMPS_GPU {

// Device-level functions (lal_device.cpp)
bool lmp_has_compatible_gpu_device();
bool lmp_gpu_requires_host_neighbor();
std::string lmp_gpu_device_info();
int lmp_init_device(MPI_Comm world, MPI_Comm replica, const int ngpu,
                    const int first_gpu_id, const int gpu_mode,
                    const int t_per_atom,
                    const double user_cell_size, char *opencl_config,
                    const int ocl_platform, char *device_type_flags,
                    const int block_pair);
void lmp_clear_device();
double lmp_gpu_forces(double **f, double **tor, double *eatom, double **vatom,
                      double *virial, double &ecoul, int &error_flag);
double lmp_gpu_update_bin_size(const double subx, const double suby, const double subz,
                               const int nlocal, const double cut);
bool lmp_gpu_config(const std::string &category, const std::string &setting);

// Per-style entry points

int amoeba_gpu_init(const int ntypes, const int max_amtype, const int max_amclass,
                    const double *host_pdamp, const double *host_thole,
                    const double *host_dirdamp, const int* host_amtype2class,
                    const double *host_special_hal, const double *host_special_repel,
                    const double *host_special_disp, const double *host_special_mpole,
                    const double *host_special_polar_wscale,
                    const double *host_special_polar_piscale,
                    const double *host_special_polar_pscale,
                    const double *host_csix, const double *host_adisp,
                    const int nlocal, const int nall, const int max_nbors,
                    const int maxspecial, const int maxspecial15,
                    const double cell_size, int &gpu_mode, FILE *screen,
                    const double polar_dscale, const double polar_uscale);
void amoeba_gpu_clear();
int** amoeba_gpu_precompute(const int ago, const int inum_full, const int nall,
                            double **host_x, int *host_type, int *host_amtype,
                            int *host_amgroup, double **host_rpole,
                            double **host_uind, double **host_uinp, double *host_pval,
                            double *sublo, double *subhi, tagint *tag,
                            int **nspecial, tagint **special,
                            int *nspecial15, tagint **special15,
                            const bool eflag_in, const bool vflag_in,
                            const bool eatom, const bool vatom, int **ilist, int **jnum, bool &success, double *host_q, double *boxlo, double *prd);
void amoeba_gpu_compute_multipole_real(const int ago, const int inum, const int nall,
              double **host_x, int *host_type, int *host_amtype, int *host_amgroup,
              double **host_rpole, double *sublo, double *subhi, tagint *tag,
              int **nspecial, tagint **special, int* nspecial15, tagint** special15,
              const bool eflag, const bool vflag, const bool eatom, const bool vatom,
              int **ilist, int **jnum, bool &success, const double aewald, const double felec, const double off2,
              double *host_q, double *boxlo, double *prd, void **tq_ptr);
void amoeba_gpu_compute_udirect2b(int *host_amtype, int *host_amgroup,
              double **host_rpole, double **host_uind, double **host_uinp,
              const double aewald, const double off2, void **fieldp_ptr);
void amoeba_gpu_compute_umutual2b(int *host_amtype, int *host_amgroup,
              double **host_rpole, double **host_uind, double **host_uinp,
              const double aewald, const double off2, void **fieldp_ptr);
void amoeba_gpu_update_fieldp(void **fieldp_ptr);
void amoeba_gpu_precompute_kspace(const int inum_full, const int bsorder,
              double ***host_thetai1, double ***host_thetai2,
              double ***host_thetai3, int** igrid,
              const int nzlo_out, const int nzhi_out,
              const int nylo_out, const int nyhi_out,
              const int nxlo_out, const int nxhi_out);
void amoeba_gpu_fphi_uind(double ****host_grid_brick, void **host_fdip_phi1,
                          void **host_fdip_phi2, void **host_fdip_sum_phi);
void amoeba_gpu_fphi_mpole(double ***host_grid_brick, void **host_fdip_sum_phi,
                           const double felec);
void amoeba_gpu_compute_polar_real(int *host_amtype, int *host_amgroup,
              double **host_rpole, double **host_uind, double **host_uinp,
              const bool eflag, const bool vflag, const bool eatom, const bool vatom,
              const double aewald, const double felec, const double off2,
              void **tq_ptr);
double amoeba_gpu_bytes();

int beck_gpu_init(const int ntypes, double **cutsq, double **host_aa, double **alpha, double **beta,
                  double **AA, double **BB, double *special_lj, const int nlocal, const int nall,
                  const int max_nbors, const int maxspecial, const double cell_size, int &gpu_mode,
                  FILE *screen);
void beck_gpu_clear();
int **beck_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                         int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                         tagint **special, const bool eflag, const bool vflag, const bool eatom,
                         const bool vatom, int **ilist, int **jnum,
                         bool &success, double *prd, int *periodicity);
void beck_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                      int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                      const bool vflag, const bool eatom, const bool vatom, bool &success);
double beck_gpu_bytes();

int bornclcs_gpu_init(const int ntypes, double **cutsq, double **host_rhoinv, double **host_born1,
                      double **host_born2, double **host_born3, double **host_a, double **host_c,
                      double **host_d, double **sigma, double **offset, double *special_lj,
                      const int inum, const int nall, const int max_nbors, const int maxspecial,
                      const double cell_size, int &gpu_mode, FILE *screen, double **host_cut_ljsq,
                      double host_cut_coulsq, double *host_special_coul, const double qqrd2e,
                      const double g_ewald);
void bornclcs_gpu_clear();
int **bornclcs_gpu_compute_n(const int ago, const int inum_full, const int nall, double **host_x,
                             int *host_type, double *sublo, double *subhi, tagint *tag,
                             int **nspecial, tagint **special, const bool eflag, const bool vflag,
                             const bool eatom, const bool vatom, int **ilist,
                             int **jnum, bool &success, double *host_q,
                             double *boxlo, double *prd, int* periodicity);
void bornclcs_gpu_compute(const int ago, const int inum_full, const int nall, double **host_x,
                          int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                          const bool vflag, const bool eatom, const bool vatom, bool &success, double *host_q, const int nlocal,
                          double *boxlo, double *prd);
double bornclcs_gpu_bytes();

int borncl_gpu_init(const int ntypes, double **cutsq, double **host_rhoinv, double **host_born1,
                    double **host_born2, double **host_born3, double **host_a, double **host_c,
                    double **host_d, double **sigma, double **offset, double *special_lj,
                    const int inum, const int nall, const int max_nbors, const int maxspecial,
                    const double cell_size, int &gpu_mode, FILE *screen, double **host_cut_ljsq,
                    double host_cut_coulsq, double *host_special_coul, const double qqrd2e,
                    const double g_ewald);
void borncl_gpu_clear();
int **borncl_gpu_compute_n(const int ago, const int inum_full, const int nall, double **host_x,
                           int *host_type, double *sublo, double *subhi, tagint *tag,
                           int **nspecial, tagint **special, const bool eflag, const bool vflag,
                           const bool eatom, const bool vatom, int **ilist,
                           int **jnum, bool &success, double *host_q,
                           double *boxlo, double *prd, int* periodicity);
void borncl_gpu_compute(const int ago, const int inum_full, const int nall, double **host_x,
                        int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                        const bool vflag, const bool eatom, const bool vatom, bool &success, double *host_q, const int nlocal,
                        double *boxlo, double *prd);
double borncl_gpu_bytes();

int borncwcs_gpu_init(const int ntypes, double **cutsq, double **host_rhoinv, double **host_born1,
                      double **host_born2, double **host_born3, double **host_a, double **host_c,
                      double **host_d, double **sigma, double **offset, double *special_lj,
                      const int inum, const int nall, const int max_nbors, const int maxspecial,
                      const double cell_size, int &gpu_mode, FILE *screen, double **host_cut_ljsq,
                      double host_cut_coulsq, double *host_special_coul, const double qqrd2e,
                      const double alf, const double e_shift, const double f_shift);
void borncwcs_gpu_clear();
int **borncwcs_gpu_compute_n(const int ago, const int inum_full, const int nall, double **host_x,
                             int *host_type, double *sublo, double *subhi, tagint *tag,
                             int **nspecial, tagint **special, const bool eflag, const bool vflag,
                             const bool eatom, const bool vatom, int **ilist,
                             int **jnum, bool &success, double *host_q,
                             double *boxlo, double *prd, int* periodicity);
void borncwcs_gpu_compute(const int ago, const int inum_full, const int nall, double **host_x,
                          int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                          const bool vflag, const bool eatom, const bool vatom, bool &success, double *host_q, const int nlocal,
                          double *boxlo, double *prd);
double borncwcs_gpu_bytes();

int borncw_gpu_init(const int ntypes, double **cutsq, double **host_rhoinv, double **host_born1,
                    double **host_born2, double **host_born3, double **host_a, double **host_c,
                    double **host_d, double **sigma, double **offset, double *special_lj,
                    const int inum, const int nall, const int max_nbors, const int maxspecial,
                    const double cell_size, int &gpu_mode, FILE *screen, double **host_cut_ljsq,
                    double host_cut_coulsq, double *host_special_coul, const double qqrd2e,
                    const double alf, const double e_shift, const double f_shift);
void borncw_gpu_clear();
int **borncw_gpu_compute_n(const int ago, const int inum_full, const int nall, double **host_x,
                           int *host_type, double *sublo, double *subhi, tagint *tag,
                           int **nspecial, tagint **special, const bool eflag, const bool vflag,
                           const bool eatom, const bool vatom, int **ilist,
                           int **jnum, bool &success, double *host_q,
                           double *boxlo, double *prd, int* periodicity);
void borncw_gpu_compute(const int ago, const int inum_full, const int nall, double **host_x,
                        int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                        const bool vflag, const bool eatom, const bool vatom, bool &success, double *host_q, const int nlocal,
                        double *boxlo, double *prd);
double borncw_gpu_bytes();

int born_gpu_init(const int ntypes, double **cutsq, double **host_rhoinv, double **host_born1,
                  double **host_born2, double **host_born3, double **host_a, double **host_c,
                  double **host_d, double **host_sigma, double **offset, double *special_lj,
                  const int inum, const int nall, const int max_nbors, const int maxspecial,
                  const double cell_size, int &gpu_mode, FILE *screen);
void born_gpu_reinit(const int ntypes, double **host_rhoinv, double **host_born1,
                     double **host_born2, double **host_born3, double **host_a, double **host_c,
                     double **host_d, double **offset);
void born_gpu_clear();
int **born_gpu_compute_n(const int ago, const int inum_full, const int nall, double **host_x,
                         int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                         tagint **special, const bool eflag, const bool vflag, const bool eatom,
                         const bool vatom, int **ilist, int **jnum,
                         bool &success, double *prd, int *periodicity);
void born_gpu_compute(const int ago, const int inum_full, const int nall, double **host_x,
                      int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                      const bool vflag, const bool eatom, const bool vatom, bool &success);
double born_gpu_bytes();

int buckc_gpu_init(const int ntypes, double **cutsq, double **host_rhoinv, double **host_buck1,
                   double **host_buck2, double **host_a, double **host_c, double **offset,
                   double *special_lj, const int inum, const int nall, const int max_nbors,
                   const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen,
                   double **host_cut_ljsq, double **host_cut_coulsq, double *host_special_coul,
                   const double qqrd2e);
void buckc_gpu_clear();
int **buckc_gpu_compute_n(const int ago, const int inum_full, const int nall, double **host_x,
                          int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                          tagint **special, const bool eflag, const bool vflag, const bool eatom,
                          const bool vatom, int **ilist, int **jnum,
                          bool &success, double *host_q, double *boxlo,
                          double *prd, int* periodicity);
void buckc_gpu_compute(const int ago, const int inum_full, const int nall, double **host_x,
                       int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                       const bool vflag, const bool eatom, const bool vatom, bool &success, double *host_q, const int nlocal,
                       double *boxlo, double *prd);
double buckc_gpu_bytes();

int buckcl_gpu_init(const int ntypes, double **cutsq, double **host_rhoinv, double **host_buck1,
                    double **host_buck2, double **host_a, double **host_c, double **offset,
                    double *special_lj, const int inum, const int nall, const int max_nbors,
                    const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen,
                    double **host_cut_ljsq, double host_cut_coulsq, double *host_special_coul,
                    const double qqrd2e, const double g_ewald);
void buckcl_gpu_clear();
int **buckcl_gpu_compute_n(const int ago, const int inum_full, const int nall, double **host_x,
                           int *host_type, double *sublo, double *subhi, tagint *tag,
                           int **nspecial, tagint **special, const bool eflag, const bool vflag,
                           const bool eatom, const bool vatom, int **ilist,
                           int **jnum, bool &success, double *host_q,
                           double *boxlo, double *prd, int *periodicity);
void buckcl_gpu_compute(const int ago, const int inum_full, const int nall, double **host_x,
                        int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                        const bool vflag, const bool eatom, const bool vatom, bool &success, double *host_q, const int nlocal,
                        double *boxlo, double *prd);
double buckcl_gpu_bytes();

int buck_gpu_init(const int ntypes, double **cutsq, double **host_rhoinv, double **host_buck1,
                  double **host_buck2, double **host_a, double **host_c, double **offset,
                  double *special_lj, const int inum, const int nall, const int max_nbors,
                  const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen);
void buck_gpu_reinit(const int ntypes, double **cutsq, double **host_rhoinv, double **host_buck1,
                     double **host_buck2, double **host_a, double **host_c, double **offset);
void buck_gpu_clear();
int **buck_gpu_compute_n(const int ago, const int inum_full, const int nall, double **host_x,
                         int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                         tagint **special, const bool eflag, const bool vflag, const bool eatom,
                         const bool vatom, int **ilist, int **jnum,
                         bool &success, double *prd, int *periodicity);
void buck_gpu_compute(const int ago, const int inum_full, const int nall, double **host_x,
                      int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                      const bool vflag, const bool eatom, const bool vatom, bool &success);
double buck_gpu_bytes();

int colloid_gpu_init(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                     double **host_lj3, double **host_lj4, double **offset, double *special_lj,
                     double **host_a12, double **host_a1, double **host_a2, double **host_d1,
                     double **host_d2, double **host_sigma3, double **host_sigma6, int **host_form,
                     const int nlocal, const int nall, const int max_nbors, const int maxspecial,
                     const double cell_size, int &gpu_mode, FILE *screen);
void colloid_gpu_clear();
int **colloid_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                            int *host_type, double *sublo, double *subhi, tagint *tag,
                            int **nspecial, tagint **special, const bool eflag, const bool vflag,
                            const bool eatom, const bool vatom, int **ilist,
                            int **jnum, bool &success, double *prd, int* periodicity);
void colloid_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                         int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                         const bool vflag, const bool eatom, const bool vatom, bool &success);
double colloid_gpu_bytes();

int coul_gpu_init(const int ntypes, double **host_scale, double **cutsq, double *special_coul,
                  const int nlocal, const int nall, const int max_nbors, const int maxspecial,
                  const double cell_size, int &gpu_mode, FILE *screen, const double qqrd2e);
void coul_gpu_reinit(const int ntypes, double **host_scale);
void coul_gpu_clear();
int **coul_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                         int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                         tagint **special, const bool eflag, const bool vflag, const bool eatom,
                         const bool vatom, int **ilist, int **jnum,
                         bool &success, double *host_q, double *boxlo,
                         double *prd, int *periodicity);
void coul_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                      int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                      const bool vflag, const bool eatom, const bool vatom, bool &success, double *host_q, const int nlocal,
                      double *boxlo, double *prd);
double coul_gpu_bytes();

int cdebye_gpu_init(const int ntypes, double **host_scale, double **cutsq, double *special_coul,
                    const int nlocal, const int nall, const int max_nbors, const int maxspecial,
                    const double cell_size, int &gpu_mode, FILE *screen, const double qqrd2e,
                    const double kappa);
void cdebye_gpu_reinit(const int ntypes, double **host_scale);
void cdebye_gpu_clear();
int **cdebye_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                           int *host_type, double *sublo, double *subhi, tagint *tag,
                           int **nspecial, tagint **special, const bool eflag, const bool vflag,
                           const bool eatom, const bool vatom, int **ilist,
                           int **jnum, bool &success, double *host_q,
                           double *boxlo, double *prd, int *periodicity);
void cdebye_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                        int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                        const bool vflag, const bool eatom, const bool vatom, bool &success, double *host_q, const int nlocal,
                        double *boxlo, double *prd);
double cdebye_gpu_bytes();

int cdsf_gpu_init(const int ntypes, const int nlocal, const int nall, const int max_nbors,
                  const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen,
                  const double host_cut_coulsq, double *host_special_coul, const double qqrd2e,
                  const double e_shift, const double f_shift, const double alpha);
void cdsf_gpu_clear();
int **cdsf_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                         int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                         tagint **special, const bool eflag, const bool vflag, const bool eatom,
                         const bool vatom, int **ilist, int **jnum,
                         bool &success, double *host_q, double *boxlo,
                         double *prd, int *periodicity);
void cdsf_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                      int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                      const bool vflag, const bool eatom, const bool vatom, bool &success, double *host_q, const int nlocal,
                      double *boxlo, double *prd);
double cdsf_gpu_bytes();

int clcs_gpu_init(const int ntypes, double **scale, const int nlocal, const int nall,
                  const int max_nbors, const int maxspecial, const double cell_size, int &gpu_mode,
                  FILE *screen, double host_cut_coulsq, double *host_special_coul,
                  const double qqrd2e, const double g_ewald);
void clcs_gpu_reinit(const int ntypes, double **scale);
void clcs_gpu_clear();
int **clcs_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                         int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                         tagint **special, const bool eflag, const bool vflag, const bool eatom,
                         const bool vatom, int **ilist, int **jnum,
                         bool &success, double *host_q, double *boxlo,
                         double *prd, int* periodicity);
void clcs_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                      int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                      const bool vflag, const bool eatom, const bool vatom, bool &success, double *host_q, const int nlocal,
                      double *boxlo, double *prd);
double clcs_gpu_bytes();

int cl_gpu_init(const int ntypes, double **scale, const int nlocal, const int nall,
                const int max_nbors, const int maxspecial, const double cell_size, int &gpu_mode,
                FILE *screen, double host_cut_coulsq, double *host_special_coul,
                const double qqrd2e, const double g_ewald);
void cl_gpu_reinit(const int ntypes, double **scale);
void cl_gpu_clear();
int **cl_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                       int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                       tagint **special, const bool eflag, const bool vflag, const bool eatom,
                       const bool vatom, int **ilist, int **jnum,
                       bool &success, double *host_q, double *boxlo,
                       double *prd, int *periodicity);
void cl_gpu_compute(const int ago, const int inum, const int nall, double **host_x, int *host_type,
                    int *ilist, int *numj, int **firstneigh, const bool eflag, const bool vflag,
                    const bool eatom, const bool vatom, bool &success, double *host_q, const int nlocal, double *boxlo, double *prd);
double cl_gpu_bytes();

int csl_gpu_init(const int ntypes, double **scale, const int nlocal, const int nall,
                const int max_nbors, const int maxspecial, const double cell_size, int &gpu_mode,
                FILE *screen, double host_cut_coulsq, double *host_special_coul,
                const double qqrd2e, const double g_ewald, const double lamda);
void csl_gpu_reinit(const int ntypes, double **scale);
void csl_gpu_clear();
int **csl_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                       int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                       tagint **special, const bool eflag, const bool vflag, const bool eatom,
                       const bool vatom, int **ilist, int **jnum,
                       bool &success, double *host_q, double *boxlo,
                       double *prd, int *periodicity);
void csl_gpu_compute(const int ago, const int inum, const int nall, double **host_x, int *host_type,
                    int *ilist, int *numj, int **firstneigh, const bool eflag, const bool vflag,
                    const bool eatom, const bool vatom, bool &success, double *host_q, const int nlocal, double *boxlo, double *prd);
double csl_gpu_bytes();

int dpd_coul_slater_long_gpu_init(const int ntypes, double **cutsq, double **host_a0,
                                  double **host_gamma, double **host_sigma, double **host_cut_dpd,
                                  double **host_cut_dpdsq, double **host_cut_slatersq,
                                  double *special_lj, const int inum, const int nall,
                                  const int max_nbors, const int maxspecial,
                                  const double cell_size, int &gpu_mode, FILE *screen,
                                  double *host_special_coul, const double qqrd2e,
                                  const double g_ewald, const double lamda);
void dpd_coul_slater_long_gpu_clear();
int **dpd_coul_slater_long_gpu_compute_n(const int ago, const int inum_full, const int nall,
                                         double **host_x, int *host_type, double *sublo,
                                         double *subhi, tagint *tag, int **nspecial,
                                         tagint **special, const bool eflag, const bool vflag,
                                         const bool eatom, const bool vatom, int **ilist, int **jnum, bool &success, double **host_v, const double dtinvsqrt,
                                         const int seed, const int timestep, double *boxlo,
                                         double *prd);
void dpd_coul_slater_long_gpu_compute(const int ago, const int inum_full, const int nall,
                                      double **host_x, int *host_type, int *ilist, int *numj,
                                      int **firstneigh, const bool eflag, const bool vflag,
                                      const bool eatom, const bool vatom, bool &success, tagint *tag,
                                      double **host_v, const double dtinvsqrt, const int seed,
                                      const int timestep, const int nlocal, double *boxlo,
                                      double *prd);
void dpd_coul_slater_long_gpu_get_extra_data(double *host_q);
double dpd_coul_slater_long_gpu_bytes();
static constexpr double EPSILON = 1.0e-10;
#if !defined(_USE_UNIFORM_SARU_LCG) && !defined(_USE_UNIFORM_SARU_TEA8) && \
    !defined(_USE_GAUSSIAN_SARU_LCG)
#define _USE_UNIFORM_SARU_LCG
#endif
#define LCGA 0x4beb5d59    // Full period 32 bit LCG
#define LCGC 0x2600e1f7
#define oWeylPeriod 0xda879add    // Prime period 3666320093
#define oWeylOffset 0x8009d14b
#define TWO_N32 0.232830643653869628906250e-9f /* 2^-32 */
#ifdef _USE_UNIFORM_SARU_LCG
#define numtyp double
#define SQRT3 (numtyp) 1.7320508075688772935274463
#define saru(seed1, seed2, seed, timestep, randnum)                                \
  {                                                                                \
    unsigned int seed3 = seed + timestep;                                          \
    seed3 ^= (seed1 << 7) ^ (seed2 >> 6);                                          \
    seed2 += (seed1 >> 4) ^ (seed3 >> 15);                                         \
    seed1 ^= (seed2 << 9) + (seed3 << 8);                                          \
    seed3 ^= 0xA5366B4D * ((seed2 >> 11) ^ (seed1 << 1));                          \
    seed2 += 0x72BE1579 * ((seed1 << 4) ^ (seed3 >> 16));                          \
    seed1 ^= 0x3F38A6ED * ((seed3 >> 5) ^ (((signed int) seed2) >> 22));           \
    seed2 += seed1 * seed3;                                                        \
    seed1 += seed3 ^ (seed2 >> 2);                                                 \
    seed2 ^= ((signed int) seed2) >> 17;                                           \
    unsigned int state = 0x79dedea3 * (seed1 ^ (((signed int) seed1) >> 14));      \
    unsigned int wstate = (state + seed2) ^ (((signed int) state) >> 8);           \
    state = state + (wstate * (wstate ^ 0xdddf97f5));                              \
    wstate = 0xABCB96F7 + (wstate >> 1);                                           \
    state = LCGA * state + LCGC;                                                   \
    wstate = wstate + oWeylOffset + ((((signed int) wstate) >> 31) & oWeylPeriod); \
    unsigned int v = (state ^ (state >> 26)) + wstate;                             \
    unsigned int s = (signed int) ((v ^ (v >> 20)) * 0x6957f5a7);                  \
    randnum = SQRT3 * (s * TWO_N32 * (numtyp) 2.0 - (numtyp) 1.0);                 \
  }
#endif
#ifdef _USE_UNIFORM_SARU_TEA8
#define numtyp double
#define SQRT3 (numtyp) 1.7320508075688772935274463
#define k0 0xA341316C
#define k1 0xC8013EA4
#define k2 0xAD90777D
#define k3 0x7E95761E
#define delta 0x9e3779b9
#define rounds 8
#define saru(seed1, seed2, seed, timestep, randnum)                           \
  {                                                                           \
    unsigned int seed3 = seed + timestep;                                     \
    seed3 ^= (seed1 << 7) ^ (seed2 >> 6);                                     \
    seed2 += (seed1 >> 4) ^ (seed3 >> 15);                                    \
    seed1 ^= (seed2 << 9) + (seed3 << 8);                                     \
    seed3 ^= 0xA5366B4D * ((seed2 >> 11) ^ (seed1 << 1));                     \
    seed2 += 0x72BE1579 * ((seed1 << 4) ^ (seed3 >> 16));                     \
    seed1 ^= 0x3F38A6ED * ((seed3 >> 5) ^ (((signed int) seed2) >> 22));      \
    seed2 += seed1 * seed3;                                                   \
    seed1 += seed3 ^ (seed2 >> 2);                                            \
    seed2 ^= ((signed int) seed2) >> 17;                                      \
    unsigned int state = 0x79dedea3 * (seed1 ^ (((signed int) seed1) >> 14)); \
    unsigned int wstate = (state + seed2) ^ (((signed int) state) >> 8);      \
    state = state + (wstate * (wstate ^ 0xdddf97f5));                         \
    wstate = 0xABCB96F7 + (wstate >> 1);                                      \
    unsigned int sum = 0;                                                     \
    for (int i = 0; i < rounds; i++) {                                        \
      sum += delta;                                                           \
      state += ((wstate << 4) + k0) ^ (wstate + sum) ^ ((wstate >> 5) + k1);  \
      wstate += ((state << 4) + k2) ^ (state + sum) ^ ((state >> 5) + k3);    \
    }                                                                         \
    unsigned int v = (state ^ (state >> 26)) + wstate;                        \
    unsigned int s = (signed int) ((v ^ (v >> 20)) * 0x6957f5a7);             \
    randnum = SQRT3 * (s * TWO_N32 * (numtyp) 2.0 - (numtyp) 1.0);            \
  }
#endif
#ifdef _USE_GAUSSIAN_SARU_LCG
#define numtyp double
#define saru(seed1, seed2, seed, timestep, randnum)                                  \
  {                                                                                  \
    unsigned int seed3 = seed + timestep;                                            \
    seed3 ^= (seed1 << 7) ^ (seed2 >> 6);                                            \
    seed2 += (seed1 >> 4) ^ (seed3 >> 15);                                           \
    seed1 ^= (seed2 << 9) + (seed3 << 8);                                            \
    seed3 ^= 0xA5366B4D * ((seed2 >> 11) ^ (seed1 << 1));                            \
    seed2 += 0x72BE1579 * ((seed1 << 4) ^ (seed3 >> 16));                            \
    seed1 ^= 0x3F38A6ED * ((seed3 >> 5) ^ (((signed int) seed2) >> 22));             \
    seed2 += seed1 * seed3;                                                          \
    seed1 += seed3 ^ (seed2 >> 2);                                                   \
    seed2 ^= ((signed int) seed2) >> 17;                                             \
    unsigned int state = 0x12345678;                                                 \
    unsigned int wstate = 12345678;                                                  \
    state = 0x79dedea3 * (seed1 ^ (((signed int) seed1) >> 14));                     \
    wstate = (state + seed2) ^ (((signed int) state) >> 8);                          \
    state = state + (wstate * (wstate ^ 0xdddf97f5));                                \
    wstate = 0xABCB96F7 + (wstate >> 1);                                             \
    unsigned int v, s;                                                               \
    numtyp r1, r2, rsq;                                                              \
    while (1) {                                                                      \
      state = LCGA * state + LCGC;                                                   \
      wstate = wstate + oWeylOffset + ((((signed int) wstate) >> 31) & oWeylPeriod); \
      v = (state ^ (state >> 26)) + wstate;                                          \
      s = (signed int) ((v ^ (v >> 20)) * 0x6957f5a7);                               \
      r1 = s * TWO_N32 * (numtyp) 2.0 - (numtyp) 1.0;                                \
      state = LCGA * state + LCGC;                                                   \
      wstate = wstate + oWeylOffset + ((((signed int) wstate) >> 31) & oWeylPeriod); \
      v = (state ^ (state >> 26)) + wstate;                                          \
      s = (signed int) ((v ^ (v >> 20)) * 0x6957f5a7);                               \
      r2 = s * TWO_N32 * (numtyp) 2.0 - (numtyp) 1.0;                                \
      rsq = r1 * r1 + r2 * r2;                                                       \
      if (rsq < (numtyp) 1.0) break;                                                 \
    }                                                                                \
    numtyp fac = sqrt((numtyp) -2.0 * log(rsq) / rsq);                               \
    randnum = r2 * fac;                                                              \
  }
#endif

int dpd_gpu_init(const int ntypes, double **cutsq, double **host_a0, double **host_gamma,
                 double **host_sigma, double **host_cut, double *special_lj, const int inum,
                 const int nall, const int max_nbors, const int maxspecial, const double cell_size,
                 int &gpu_mode, FILE *screen);
void dpd_gpu_clear();
int **dpd_gpu_compute_n(const int ago, const int inum_full, const int nall, double **host_x,
                        int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                        tagint **special, const bool eflag, const bool vflag, const bool eatom,
                        const bool vatom, int **ilist, int **jnum,
                        bool &success, double **host_v,
                        const double dtinvsqrt, const int seed, const int timestep, double *boxlo,
                        double *prd);
void dpd_gpu_compute(const int ago, const int inum_full, const int nall, double **host_x,
                     int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                     const bool vflag, const bool eatom, const bool vatom, bool &success, tagint *tag, double **host_v,
                     const double dtinvsqrt, const int seed, const int timestep, const int nlocal,
                     double *boxlo, double *prd);
double dpd_gpu_bytes();
static constexpr double EPSILON = 1.0e-10;
#if !defined(_USE_UNIFORM_SARU_LCG) && !defined(_USE_UNIFORM_SARU_TEA8) && \
    !defined(_USE_GAUSSIAN_SARU_LCG)
#define _USE_UNIFORM_SARU_LCG
#endif
#define LCGA 0x4beb5d59    // Full period 32 bit LCG
#define LCGC 0x2600e1f7
#define oWeylPeriod 0xda879add    // Prime period 3666320093
#define oWeylOffset 0x8009d14b
#define TWO_N32 0.232830643653869628906250e-9f /* 2^-32 */
#ifdef _USE_UNIFORM_SARU_LCG
#define numtyp double
#define SQRT3 (numtyp) 1.7320508075688772935274463
#define saru(seed1, seed2, seed, timestep, randnum)                                \
  {                                                                                \
    unsigned int seed3 = seed + timestep;                                          \
    seed3 ^= (seed1 << 7) ^ (seed2 >> 6);                                          \
    seed2 += (seed1 >> 4) ^ (seed3 >> 15);                                         \
    seed1 ^= (seed2 << 9) + (seed3 << 8);                                          \
    seed3 ^= 0xA5366B4D * ((seed2 >> 11) ^ (seed1 << 1));                          \
    seed2 += 0x72BE1579 * ((seed1 << 4) ^ (seed3 >> 16));                          \
    seed1 ^= 0x3F38A6ED * ((seed3 >> 5) ^ (((signed int) seed2) >> 22));           \
    seed2 += seed1 * seed3;                                                        \
    seed1 += seed3 ^ (seed2 >> 2);                                                 \
    seed2 ^= ((signed int) seed2) >> 17;                                           \
    unsigned int state = 0x79dedea3 * (seed1 ^ (((signed int) seed1) >> 14));      \
    unsigned int wstate = (state + seed2) ^ (((signed int) state) >> 8);           \
    state = state + (wstate * (wstate ^ 0xdddf97f5));                              \
    wstate = 0xABCB96F7 + (wstate >> 1);                                           \
    state = LCGA * state + LCGC;                                                   \
    wstate = wstate + oWeylOffset + ((((signed int) wstate) >> 31) & oWeylPeriod); \
    unsigned int v = (state ^ (state >> 26)) + wstate;                             \
    unsigned int s = (signed int) ((v ^ (v >> 20)) * 0x6957f5a7);                  \
    randnum = SQRT3 * (s * TWO_N32 * (numtyp) 2.0 - (numtyp) 1.0);                 \
  }
#endif
#ifdef _USE_UNIFORM_SARU_TEA8
#define numtyp double
#define SQRT3 (numtyp) 1.7320508075688772935274463
#define k0 0xA341316C
#define k1 0xC8013EA4
#define k2 0xAD90777D
#define k3 0x7E95761E
#define delta 0x9e3779b9
#define rounds 8
#define saru(seed1, seed2, seed, timestep, randnum)                           \
  {                                                                           \
    unsigned int seed3 = seed + timestep;                                     \
    seed3 ^= (seed1 << 7) ^ (seed2 >> 6);                                     \
    seed2 += (seed1 >> 4) ^ (seed3 >> 15);                                    \
    seed1 ^= (seed2 << 9) + (seed3 << 8);                                     \
    seed3 ^= 0xA5366B4D * ((seed2 >> 11) ^ (seed1 << 1));                     \
    seed2 += 0x72BE1579 * ((seed1 << 4) ^ (seed3 >> 16));                     \
    seed1 ^= 0x3F38A6ED * ((seed3 >> 5) ^ (((signed int) seed2) >> 22));      \
    seed2 += seed1 * seed3;                                                   \
    seed1 += seed3 ^ (seed2 >> 2);                                            \
    seed2 ^= ((signed int) seed2) >> 17;                                      \
    unsigned int state = 0x79dedea3 * (seed1 ^ (((signed int) seed1) >> 14)); \
    unsigned int wstate = (state + seed2) ^ (((signed int) state) >> 8);      \
    state = state + (wstate * (wstate ^ 0xdddf97f5));                         \
    wstate = 0xABCB96F7 + (wstate >> 1);                                      \
    unsigned int sum = 0;                                                     \
    for (int i = 0; i < rounds; i++) {                                        \
      sum += delta;                                                           \
      state += ((wstate << 4) + k0) ^ (wstate + sum) ^ ((wstate >> 5) + k1);  \
      wstate += ((state << 4) + k2) ^ (state + sum) ^ ((state >> 5) + k3);    \
    }                                                                         \
    unsigned int v = (state ^ (state >> 26)) + wstate;                        \
    unsigned int s = (signed int) ((v ^ (v >> 20)) * 0x6957f5a7);             \
    randnum = SQRT3 * (s * TWO_N32 * (numtyp) 2.0 - (numtyp) 1.0);            \
  }
#endif
#ifdef _USE_GAUSSIAN_SARU_LCG
#define numtyp double
#define saru(seed1, seed2, seed, timestep, randnum)                                  \
  {                                                                                  \
    unsigned int seed3 = seed + timestep;                                            \
    seed3 ^= (seed1 << 7) ^ (seed2 >> 6);                                            \
    seed2 += (seed1 >> 4) ^ (seed3 >> 15);                                           \
    seed1 ^= (seed2 << 9) + (seed3 << 8);                                            \
    seed3 ^= 0xA5366B4D * ((seed2 >> 11) ^ (seed1 << 1));                            \
    seed2 += 0x72BE1579 * ((seed1 << 4) ^ (seed3 >> 16));                            \
    seed1 ^= 0x3F38A6ED * ((seed3 >> 5) ^ (((signed int) seed2) >> 22));             \
    seed2 += seed1 * seed3;                                                          \
    seed1 += seed3 ^ (seed2 >> 2);                                                   \
    seed2 ^= ((signed int) seed2) >> 17;                                             \
    unsigned int state = 0x12345678;                                                 \
    unsigned int wstate = 12345678;                                                  \
    state = 0x79dedea3 * (seed1 ^ (((signed int) seed1) >> 14));                     \
    wstate = (state + seed2) ^ (((signed int) state) >> 8);                          \
    state = state + (wstate * (wstate ^ 0xdddf97f5));                                \
    wstate = 0xABCB96F7 + (wstate >> 1);                                             \
    unsigned int v, s;                                                               \
    numtyp r1, r2, rsq;                                                              \
    while (1) {                                                                      \
      state = LCGA * state + LCGC;                                                   \
      wstate = wstate + oWeylOffset + ((((signed int) wstate) >> 31) & oWeylPeriod); \
      v = (state ^ (state >> 26)) + wstate;                                          \
      s = (signed int) ((v ^ (v >> 20)) * 0x6957f5a7);                               \
      r1 = s * TWO_N32 * (numtyp) 2.0 - (numtyp) 1.0;                                \
      state = LCGA * state + LCGC;                                                   \
      wstate = wstate + oWeylOffset + ((((signed int) wstate) >> 31) & oWeylPeriod); \
      v = (state ^ (state >> 26)) + wstate;                                          \
      s = (signed int) ((v ^ (v >> 20)) * 0x6957f5a7);                               \
      r2 = s * TWO_N32 * (numtyp) 2.0 - (numtyp) 1.0;                                \
      rsq = r1 * r1 + r2 * r2;                                                       \
      if (rsq < (numtyp) 1.0) break;                                                 \
    }                                                                                \
    numtyp fac = sqrt((numtyp) -2.0 * log(rsq) / rsq);                               \
    randnum = r2 * fac;                                                              \
  }
#endif

int dpd_tstat_gpu_init(const int ntypes, double **cutsq, double **host_a0, double **host_gamma,
                       double **host_sigma, double **host_cut, double *special_lj, const int inum,
                       const int nall, const int max_nbors, const int maxspecial,
                       const double cell_size, int &gpu_mode, FILE *screen);
void dpd_tstat_gpu_clear();
int **dpd_tstat_gpu_compute_n(const int ago, const int inum_full, const int nall, double **host_x,
                              int *host_type, double *sublo, double *subhi, tagint *tag,
                              int **nspecial, tagint **special, const bool eflag, const bool vflag,
                              const bool eatom, const bool vatom, int **ilist,
                              int **jnum, bool &success, double **host_v,
                              const double dtinvsqrt, const int seed, const int timestep,
                              double *boxlo, double *prd);
void dpd_tstat_gpu_compute(const int ago, const int inum_full, const int nall, double **host_x,
                           int *host_type, int *ilist, int *numj, int **firstneigh,
                           const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                           bool &success, tagint *tag,
                           double **host_v, const double dtinvsqrt, const int seed,
                           const int timestep, const int nlocal, double *boxlo, double *prd);
void dpd_tstat_gpu_update_coeff(int ntypes, double **host_a0, double **host_gamma,
                                double **host_sigma, double **host_cut);
double dpd_tstat_gpu_bytes();
static constexpr double EPSILON = 1.0e-10;
#if !defined(_USE_UNIFORM_SARU_LCG) && !defined(_USE_UNIFORM_SARU_TEA8) && \
    !defined(_USE_GAUSSIAN_SARU_LCG)
#define _USE_UNIFORM_SARU_LCG
#endif
#define LCGA 0x4beb5d59    // Full period 32 bit LCG
#define LCGC 0x2600e1f7
#define oWeylPeriod 0xda879add    // Prime period 3666320093
#define oWeylOffset 0x8009d14b
#define TWO_N32 0.232830643653869628906250e-9f /* 2^-32 */
#ifdef _USE_UNIFORM_SARU_LCG
#define numtyp double
#define SQRT3 (numtyp) 1.7320508075688772935274463
#define saru(seed1, seed2, seed, timestep, randnum)                                \
  {                                                                                \
    unsigned int seed3 = seed + timestep;                                          \
    seed3 ^= (seed1 << 7) ^ (seed2 >> 6);                                          \
    seed2 += (seed1 >> 4) ^ (seed3 >> 15);                                         \
    seed1 ^= (seed2 << 9) + (seed3 << 8);                                          \
    seed3 ^= 0xA5366B4D * ((seed2 >> 11) ^ (seed1 << 1));                          \
    seed2 += 0x72BE1579 * ((seed1 << 4) ^ (seed3 >> 16));                          \
    seed1 ^= 0x3F38A6ED * ((seed3 >> 5) ^ (((signed int) seed2) >> 22));           \
    seed2 += seed1 * seed3;                                                        \
    seed1 += seed3 ^ (seed2 >> 2);                                                 \
    seed2 ^= ((signed int) seed2) >> 17;                                           \
    unsigned int state = 0x79dedea3 * (seed1 ^ (((signed int) seed1) >> 14));      \
    unsigned int wstate = (state + seed2) ^ (((signed int) state) >> 8);           \
    state = state + (wstate * (wstate ^ 0xdddf97f5));                              \
    wstate = 0xABCB96F7 + (wstate >> 1);                                           \
    state = LCGA * state + LCGC;                                                   \
    wstate = wstate + oWeylOffset + ((((signed int) wstate) >> 31) & oWeylPeriod); \
    unsigned int v = (state ^ (state >> 26)) + wstate;                             \
    unsigned int s = (signed int) ((v ^ (v >> 20)) * 0x6957f5a7);                  \
    randnum = SQRT3 * (s * TWO_N32 * (numtyp) 2.0 - (numtyp) 1.0);                 \
  }
#endif
#ifdef _USE_UNIFORM_SARU_TEA8
#define numtyp double
#define SQRT3 (numtyp) 1.7320508075688772935274463
#define k0 0xA341316C
#define k1 0xC8013EA4
#define k2 0xAD90777D
#define k3 0x7E95761E
#define delta 0x9e3779b9
#define rounds 8
#define saru(seed1, seed2, seed, timestep, randnum)                           \
  {                                                                           \
    unsigned int seed3 = seed + timestep;                                     \
    seed3 ^= (seed1 << 7) ^ (seed2 >> 6);                                     \
    seed2 += (seed1 >> 4) ^ (seed3 >> 15);                                    \
    seed1 ^= (seed2 << 9) + (seed3 << 8);                                     \
    seed3 ^= 0xA5366B4D * ((seed2 >> 11) ^ (seed1 << 1));                     \
    seed2 += 0x72BE1579 * ((seed1 << 4) ^ (seed3 >> 16));                     \
    seed1 ^= 0x3F38A6ED * ((seed3 >> 5) ^ (((signed int) seed2) >> 22));      \
    seed2 += seed1 * seed3;                                                   \
    seed1 += seed3 ^ (seed2 >> 2);                                            \
    seed2 ^= ((signed int) seed2) >> 17;                                      \
    unsigned int state = 0x79dedea3 * (seed1 ^ (((signed int) seed1) >> 14)); \
    unsigned int wstate = (state + seed2) ^ (((signed int) state) >> 8);      \
    state = state + (wstate * (wstate ^ 0xdddf97f5));                         \
    wstate = 0xABCB96F7 + (wstate >> 1);                                      \
    unsigned int sum = 0;                                                     \
    for (int i = 0; i < rounds; i++) {                                        \
      sum += delta;                                                           \
      state += ((wstate << 4) + k0) ^ (wstate + sum) ^ ((wstate >> 5) + k1);  \
      wstate += ((state << 4) + k2) ^ (state + sum) ^ ((state >> 5) + k3);    \
    }                                                                         \
    unsigned int v = (state ^ (state >> 26)) + wstate;                        \
    unsigned int s = (signed int) ((v ^ (v >> 20)) * 0x6957f5a7);             \
    randnum = SQRT3 * (s * TWO_N32 * (numtyp) 2.0 - (numtyp) 1.0);            \
  }
#endif
#ifdef _USE_GAUSSIAN_SARU_LCG
#define numtyp double
#define saru(seed1, seed2, seed, timestep, randnum)                                  \
  {                                                                                  \
    unsigned int seed3 = seed + timestep;                                            \
    seed3 ^= (seed1 << 7) ^ (seed2 >> 6);                                            \
    seed2 += (seed1 >> 4) ^ (seed3 >> 15);                                           \
    seed1 ^= (seed2 << 9) + (seed3 << 8);                                            \
    seed3 ^= 0xA5366B4D * ((seed2 >> 11) ^ (seed1 << 1));                            \
    seed2 += 0x72BE1579 * ((seed1 << 4) ^ (seed3 >> 16));                            \
    seed1 ^= 0x3F38A6ED * ((seed3 >> 5) ^ (((signed int) seed2) >> 22));             \
    seed2 += seed1 * seed3;                                                          \
    seed1 += seed3 ^ (seed2 >> 2);                                                   \
    seed2 ^= ((signed int) seed2) >> 17;                                             \
    unsigned int state = 0x12345678;                                                 \
    unsigned int wstate = 12345678;                                                  \
    state = 0x79dedea3 * (seed1 ^ (((signed int) seed1) >> 14));                     \
    wstate = (state + seed2) ^ (((signed int) state) >> 8);                          \
    state = state + (wstate * (wstate ^ 0xdddf97f5));                                \
    wstate = 0xABCB96F7 + (wstate >> 1);                                             \
    unsigned int v, s;                                                               \
    numtyp r1, r2, rsq;                                                              \
    while (1) {                                                                      \
      state = LCGA * state + LCGC;                                                   \
      wstate = wstate + oWeylOffset + ((((signed int) wstate) >> 31) & oWeylPeriod); \
      v = (state ^ (state >> 26)) + wstate;                                          \
      s = (signed int) ((v ^ (v >> 20)) * 0x6957f5a7);                               \
      r1 = s * TWO_N32 * (numtyp) 2.0 - (numtyp) 1.0;                                \
      state = LCGA * state + LCGC;                                                   \
      wstate = wstate + oWeylOffset + ((((signed int) wstate) >> 31) & oWeylPeriod); \
      v = (state ^ (state >> 26)) + wstate;                                          \
      s = (signed int) ((v ^ (v >> 20)) * 0x6957f5a7);                               \
      r2 = s * TWO_N32 * (numtyp) 2.0 - (numtyp) 1.0;                                \
      rsq = r1 * r1 + r2 * r2;                                                       \
      if (rsq < (numtyp) 1.0) break;                                                 \
    }                                                                                \
    numtyp fac = sqrt((numtyp) -2.0 * log(rsq) / rsq);                               \
    randnum = r2 * fac;                                                              \
  }
#endif

int eam_alloy_gpu_init(const int ntypes, double host_cutforcesq, int **host_type2rhor,
                       int **host_type2z2r, int *host_type2frho, double ***host_rhor_spline,
                       double ***host_z2r_spline, double ***host_frho_spline, double **host_cutsq,
                       double rdr, double rdrho, double rhomax, int nrhor, int nrho, int nz2r,
                       int nfrho, int nr, const int nlocal, const int nall, const int max_nbors,
                       const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen,
                       int &fp_size);
void eam_alloy_gpu_clear();
int **eam_alloy_gpu_compute_n(const int ago, const int inum_full, const int nall, double **host_x,
                              int *host_type, double *sublo, double *subhi, tagint *tag,
                              int **nspecial, tagint **special, const bool eflag, const bool vflag,
                              const bool eatom, const bool vatom, int **ilist,
                              int **jnum, bool &success, int &inum,
                              void **fp_ptr, double *prd, int *periodicity);
void eam_alloy_gpu_compute(const int ago, const int inum_full, const int nlocal, const int nall,
                           double **host_x, int *host_type, int *ilist, int *numj, int **firstneigh,
                           const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                           bool &success, void **fp_ptr);
void eam_alloy_gpu_compute_force(int *ilist, const bool eflag, const bool vflag, const bool eatom,
                                 const bool vatom);
double eam_alloy_gpu_bytes();

int eam_fs_gpu_init(const int ntypes, double host_cutforcesq, int **host_type2rhor,
                    int **host_type2z2r, int *host_type2frho, double ***host_rhor_spline,
                    double ***host_z2r_spline, double ***host_frho_spline, double **host_cutsq,
                    double rdr, double rdrho, double rhomax, int nrhor, int nrho, int nz2r,
                    int nfrho, int nr, const int nlocal, const int nall, const int max_nbors,
                    const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen,
                    int &fp_size);
void eam_fs_gpu_clear();
int **eam_fs_gpu_compute_n(const int ago, const int inum_full, const int nall, double **host_x,
                           int *host_type, double *sublo, double *subhi, tagint *tag,
                           int **nspecial, tagint **special, const bool eflag, const bool vflag,
                           const bool eatom, const bool vatom, int **ilist,
                           int **jnum, bool &success, int &inum,
                           void **fp_ptr, double *prd, int *periodicity);
void eam_fs_gpu_compute(const int ago, const int inum_full, const int nlocal, const int nall,
                        double **host_x, int *host_type, int *ilist, int *numj, int **firstneigh,
                        const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                        bool &success, void **fp_ptr);
void eam_fs_gpu_compute_force(int *ilist, const bool eflag, const bool vflag, const bool eatom,
                              const bool vatom);
double eam_fs_gpu_bytes();

int eam_gpu_init(const int ntypes, double host_cutforcesq, int **host_type2rhor,
                 int **host_type2z2r, int *host_type2frho, double ***host_rhor_spline,
                 double ***host_z2r_spline, double ***host_frho_spline, double **host_cutsq,
                 double rdr, double rdrho, double rhomax, int nrhor, int nrho, int nz2r, int nfrho,
                 int nr, const int nlocal, const int nall, const int max_nbors,
                 const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen,
                 int &fp_size);
void eam_gpu_clear();
int **eam_gpu_compute_n(const int ago, const int inum_full, const int nall, double **host_x,
                        int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                        tagint **special, const bool eflag, const bool vflag, const bool eatom,
                        const bool vatom, int **ilist, int **jnum,
                        bool &success, int &inum, void **fp_ptr,
                        double *prd, int *periodicity);
void eam_gpu_compute(const int ago, const int inum_full, const int nlocal, const int nall,
                     double **host_x, int *host_type, int *ilist, int *numj, int **firstneigh,
                     const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                     bool &success, void **fp_ptr);
void eam_gpu_compute_force(int *ilist, const bool eflag, const bool vflag, const bool eatom,
                           const bool vatom);
double eam_gpu_bytes();

int edpd_gpu_init(const int ntypes, double **cutsq, double **host_a0, double **host_gamma,
                  double **host_cut, double **host_power, double **host_kappa,
                  double **host_powerT, double** host_cutT, double*** host_sc, double ***host_kc,
                  double *host_mass, double *special_lj, const int power_flag, const int kappa_flag,
                  const int inum, const int nall, const int max_nbors,
                  const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen);
void edpd_gpu_clear();
int **edpd_gpu_compute_n(const int ago, const int inum_full, const int nall, double **host_x,
                        int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                        tagint **special, const bool eflag, const bool vflag, const bool eatom,
                        const bool vatom, int **ilist, int **jnum,
                        bool &success, double **host_v,
                        const double dtinvsqrt, const int seed, const int timestep, double *boxlo,
                        double *prd);
void edpd_gpu_compute(const int ago, const int inum_full, const int nall, double **host_x,
                     int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                     const bool vflag, const bool eatom, const bool vatom, bool &success, tagint *tag, double **host_v,
                     const double dtinvsqrt, const int seed, const int timestep, const int nlocal,
                     double *boxlo, double *prd);
void edpd_gpu_get_extra_data(double *host_T, double *host_cv);
void edpd_gpu_update_flux(void **flux_ptr);
double edpd_gpu_bytes();

int gauss_gpu_init(const int ntypes, double **cutsq, double **host_a, double **b, double **offset,
                   double *special_lj, const int nlocal, const int nall, const int max_nbors,
                   const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen);
void gauss_gpu_reinit(const int ntypes, double **cutsq, double **host_a, double **b,
                      double **offset);
void gauss_gpu_clear();
int **gauss_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                          int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                          tagint **special, const bool eflag, const bool vflag, const bool eatom,
                          const bool vatom, int **ilist, int **jnum,
                          bool &success, double *prd, int *periodicity);
void gauss_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                       int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                       const bool vflag, const bool eatom, const bool vatom, bool &success);
double gauss_gpu_bytes();

int gb_gpu_init(const int ntypes, const double gamma, const double upsilon,
                const double mu, double **shape, double **well, double **cutsq,
                double **sigma, double **epsilon, double *host_lshape,
                int **form, double **host_lj1, double **host_lj2,
                double **host_lj3, double **host_lj4, double **offset,
                double *special_lj, const int nlocal, const int nall,
                const int max_nbors, const int maxspecial,
                const double cell_size, int &gpu_mode, FILE *screen);
void gb_gpu_clear();
int **gb_gpu_compute_n(const int ago, const int inum, const int nall,
                       double **host_x, int *host_type, double *sublo,
                       double *subhi, tagint *tag, int **nspecial,
                       tagint **special, const bool eflag, const bool vflag,
                       const bool eatom, const bool vatom, int **ilist, int **jnum, bool &success, const int *ellipsoid,
                       const void *bonus);
int *gb_gpu_compute(const int ago, const int inum, const int nall,
                    double **host_x, int *host_type, int *ilist, int *numj,
                    int **firstneigh, const bool eflag, const bool vflag,
                    const bool eatom, const bool vatom, bool &success, const int *ellipsoid,
                    const void *bonus);
double gb_gpu_bytes();
enum { SPHERE_SPHERE, SPHERE_ELLIPSE, ELLIPSE_SPHERE, ELLIPSE_ELLIPSE };

int hippo_gpu_init(const int ntypes, const int max_amtype, const int max_amclass,
                    const double *host_pdamp, const double *host_thole,
                    const double *host_dirdamp, const int* host_amtype2class,
                    const double *host_special_repel, const double *host_special_disp,
                    const double *host_special_mpole,
                    const double *host_special_polar_wscale,
                    const double *host_special_polar_piscale,
                    const double *host_special_polar_pscale,
                    const double *host_sizpr, const double *host_dmppr, const double *host_elepr,
                    const double *host_csix, const double *host_adisp,
                    const double *host_pcore, const double *host_palpha,
                    const int nlocal, const int nall, const int max_nbors,
                    const int maxspecial, const int maxspecial15,
                    const double cell_size, int &gpu_mode, FILE *screen,
                    const double polar_dscale, const double polar_uscale);
void hippo_gpu_clear();
int** hippo_gpu_precompute(const int ago, const int inum_full, const int nall,
                            double **host_x, int *host_type, int *host_amtype,
                            int *host_amgroup, double **host_rpole,
                            double **host_uind, double **host_uinp, double *host_pval,
                            double *sublo, double *subhi, tagint *tag,
                            int **nspecial, tagint **special,
                            int *nspecial15, tagint **special15,
                            const bool eflag_in, const bool vflag_in,
                            const bool eatom, const bool vatom, int **ilist, int **jnum, bool &success, double *host_q, double *boxlo, double *prd);
void hippo_gpu_compute_repulsion(const int ago, const int inum_full,
                           const int nall, double **host_x, int *host_type,
                           int *host_amtype, int *host_amgroup, double **host_rpole,
                           double *sublo, double *subhi, tagint *tag, int **nspecial,
                           tagint **special, int *nspecial15, tagint** special15,
                           const bool eflag, const bool vflag, const bool eatom,
                           const bool vatom, int **ilist, int **jnum, bool &success, const double aewald, const double off2,
                           double *host_q, double *boxlo, double *prd,
                           double cut2, double c0, double c1, double c2,
                           double c3, double c4, double c5, void **tep_ptr);
void hippo_gpu_compute_dispersion_real(int *host_amtype, int *host_amgroup, double **host_rpole,
                                        const double aewald, const double off2);
void hippo_gpu_compute_multipole_real(const int ago, const int inum, const int nall,
              double **host_x, int *host_type, int *host_amtype, int *host_amgroup,
              double **host_rpole, double *host_pval, double *sublo, double *subhi, tagint *tag,
              int **nspecial, tagint **special, int* nspecial15, tagint** special15,
              const bool eflag, const bool vflag, const bool eatom, const bool vatom,
              int **ilist, int **jnum, bool &success, const double aewald, const double felec, const double off2,
              double *host_q, double *boxlo, double *prd, void **tq_ptr);
void hippo_gpu_compute_udirect2b(int *host_amtype, int *host_amgroup,
              double **host_rpole, double **host_uind, double **host_uinp,
              double *host_pval, const double aewald, const double off2, void **fieldp_ptr);
void hippo_gpu_compute_umutual2b(int *host_amtype, int *host_amgroup,
              double **host_rpole, double **host_uind, double **host_uinp, double *host_pval,
              const double aewald, const double off2, void **fieldp_ptr);
void hippo_gpu_update_fieldp(void **fieldp_ptr);
void hippo_gpu_precompute_kspace(const int inum_full, const int bsorder,
              double ***host_thetai1, double ***host_thetai2,
              double ***host_thetai3, int** igrid,
              const int nzlo_out, const int nzhi_out,
              const int nylo_out, const int nyhi_out,
              const int nxlo_out, const int nxhi_out);
void hippo_gpu_fphi_uind(double ****host_grid_brick, void **host_fdip_phi1,
                          void **host_fdip_phi2, void **host_fdip_sum_phi);
void hippo_gpu_compute_polar_real(int *host_amtype, int *host_amgroup,
              double **host_rpole, double **host_uind, double **host_uinp, double *host_pval,
              const bool eflag, const bool vflag, const bool eatom, const bool vatom,
              const double aewald, const double felec, const double off2,
              void **tq_ptr);
double hippo_gpu_bytes();

int lj96_gpu_init(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                  double **host_lj3, double **host_lj4, double **offset, double *special_lj,
                  const int nlocal, const int nall, const int max_nbors, const int maxspecial,
                  const double cell_size, int &gpu_mode, FILE *screen);
void lj96_gpu_clear();
int **lj96_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                         int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                         tagint **special, const bool eflag, const bool vflag, const bool eatom,
                         const bool vatom, int **ilist, int **jnum,
                         bool &success, double *prd, int *periodicity);
void lj96_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                      int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                      const bool vflag, const bool eatom, const bool vatom, bool &success);
double lj96_gpu_bytes();

int crm_gpu_init(const int ntypes, double cut_bothsq, double **host_lj1, double **host_lj2,
                 double **host_lj3, double **host_lj4, double *special_lj, const int nlocal,
                 const int nall, const int max_nbors, const int maxspecial, const double cell_size,
                 int &gpu_mode, FILE *screen, double host_cut_ljsq, double host_cut_coulsq,
                 double *host_special_coul, const double qqrd2e, const double cut_lj_innersq,
                 const double cut_coul_innersq, const double denom_lj, const double denom_coul,
                 double **epsilon, double **sigma, const bool mix_arithmetic);
void crm_gpu_clear();
int **crm_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                        int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                        tagint **special, const bool eflag, const bool vflag, const bool eatom,
                        const bool vatom, int **ilist, int **jnum,
                        bool &success, double *host_q, double *boxlo,
                        double *prd, int* periodicity);
void crm_gpu_compute(const int ago, const int inum, const int nall, double **host_x, int *host_type,
                     int *ilist, int *numj, int **firstneigh, const bool eflag, const bool vflag,
                     const bool eatom, const bool vatom, bool &success, double *host_q, const int nlocal, double *boxlo, double *prd);
double crm_gpu_bytes();

int crml_gpu_init(const int ntypes, double cut_bothsq, double **host_lj1, double **host_lj2,
                  double **host_lj3, double **host_lj4, double **offset, double *special_lj,
                  const int nlocal, const int nall, const int max_nbors, const int maxspecial,
                  const double cell_size, int &gpu_mode, FILE *screen, double host_cut_ljsq,
                  double host_cut_coulsq, double *host_special_coul, const double qqrd2e,
                  const double g_ewald, const double cut_lj_innersq, const double denom_lj,
                  double **epsilon, double **sigma, const bool mix_arithmetic);
void crml_gpu_clear();
int **crml_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                         int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                         tagint **special, const bool eflag, const bool vflag, const bool eatom,
                         const bool vatom, int **ilist, int **jnum,
                         bool &success, double *host_q, double *boxlo,
                         double *prd, int* periodicity);
void crml_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                      int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                      const bool vflag, const bool eatom, const bool vatom, bool &success, double *host_q, const int nlocal,
                      double *boxlo, double *prd);
double crml_gpu_bytes();

int c2cl_gpu_init(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                  double **host_lj3, double **host_lj4, double **offset, double *special_lj,
                  const int nlocal, const int nall, const int max_nbors, const int maxspecial,
                  const double cell_size, int &gpu_mode, FILE *screen, double **host_cut_ljsq,
                  double host_cut_coulsq, double *host_special_coul, const double qqrd2e,
                  const double g_ewald);
void c2cl_gpu_clear();
int **c2cl_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                         int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                         tagint **special, const bool eflag, const bool vflag, const bool eatom,
                         const bool vatom, int **ilist, int **jnum,
                         bool &success, double *host_q, double *boxlo,
                         double *prd, int *periodicity);
void c2cl_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                      int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                      const bool vflag, const bool eatom, const bool vatom, bool &success, double *host_q, const int nlocal,
                      double *boxlo, double *prd);
double c2cl_gpu_bytes();

int lj96_gpu_init(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                  double **host_lj3, double **host_lj4, double **offset, double *special_lj,
                  const int nlocal, const int nall, const int max_nbors, const int maxspecial,
                  const double cell_size, int &gpu_mode, FILE *screen);
void lj96_gpu_clear();
int **lj96_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                         int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                         tagint **special, const bool eflag, const bool vflag, const bool eatom,
                         const bool vatom, int **ilist, int **jnum,
                         bool &success, double *prd, int *periodicity);
void lj96_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                      int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                      const bool vflag, const bool eatom, const bool vatom, bool &success);
double lj96_gpu_bytes();

int ljcb_gpu_init(const int ntypes, double **cutsq, double **cut_inner_sq, double **cut_inner,
                  double **sigma, double **epsilon, double **host_lj1, double **host_lj2,
                  double **host_lj3, double **host_lj4, double *special_lj, const int nlocal,
                  const int nall, const int max_nbors, const int maxspecial, const double cell_size,
                  int &gpu_mode, FILE *screen);
void ljcb_gpu_clear();
int **ljcb_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                         int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                         tagint **special, const bool eflag, const bool vflag, const bool eatom,
                         const bool vatom, int **ilist, int **jnum,
                         bool &success, double *prd, int *periodicity);
void ljcb_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                      int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                      const bool vflag, const bool eatom, const bool vatom, bool &success);
double ljcb_gpu_bytes();

int ljc_gpu_init(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                 double **host_lj3, double **host_lj4, double **offset, double *special_lj,
                 const int nlocal, const int nall, const int max_nbors, const int maxspecial,
                 const double cell_size, int &gpu_mode, FILE *screen, double **host_cut_ljsq,
                 double **host_cut_coulsq, double *host_special_coul, const double qqrd2e);
void ljc_gpu_clear();
int **ljc_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                        int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                        tagint **special, const bool eflag, const bool vflag, const bool eatom,
                        const bool vatom, int **ilist, int **jnum,
                        bool &success, double *host_q, double *boxlo,
                        double *prd, int *periodicity);
void ljc_gpu_compute(const int ago, const int inum, const int nall, double **host_x, int *host_type,
                     int *ilist, int *numj, int **firstneigh, const bool eflag, const bool vflag,
                     const bool eatom, const bool vatom, bool &success, double *host_q, const int nlocal, double *boxlo, double *prd);
double ljc_gpu_bytes();

int ljcs_gpu_init(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                 double **host_lj3, double **host_lj4, double **offset, double **epsilon, double *special_lj,
                 const int nlocal, const int nall, const int max_nbors, const int maxspecial,
                 const double cell_size, int &gpu_mode, FILE *screen, double **host_cut_ljsq,
                 double **host_cut_coulsq, double *host_special_coul, const double qqrd2e);
void ljcs_gpu_clear();
int **ljcs_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                        int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                        tagint **special, const bool eflag, const bool vflag, const bool eatom,
                        const bool vatom, int **ilist, int **jnum,
                        bool &success, double *host_q, double *boxlo,
                        double *prd, int* periodicity);
void ljcs_gpu_compute(const int ago, const int inum, const int nall, double **host_x, int *host_type,
                     int *ilist, int *numj, int **firstneigh, const bool eflag, const bool vflag,
                     const bool eatom, const bool vatom, bool &success, double *host_q, const int nlocal, double *boxlo, double *prd);
double ljcs_gpu_bytes();

int ljcd_gpu_init(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                  double **host_lj3, double **host_lj4, double **offset, double *special_lj,
                  const int nlocal, const int nall, const int max_nbors, const int maxspecial,
                  const double cell_size, int &gpu_mode, FILE *screen, double **host_cut_ljsq,
                  double **host_cut_coulsq, double *host_special_coul, const double qqrd2e,
                  const double kappa);
void ljcd_gpu_clear();
int **ljcd_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                         int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                         tagint **special, const bool eflag, const bool vflag, const bool eatom,
                         const bool vatom, int **ilist, int **jnum,
                         bool &success, double *host_q, double *boxlo,
                         double *prd, int *periodicity);
void ljcd_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                      int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                      const bool vflag, const bool eatom, const bool vatom, bool &success, double *host_q, const int nlocal,
                      double *boxlo, double *prd);
double ljcd_gpu_bytes();

int ljd_gpu_init(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                 double **host_lj3, double **host_lj4, double **offset, double *special_lj,
                 const int nlocal, const int nall, const int max_nbors, const int maxspecial,
                 const double cell_size, int &gpu_mode, FILE *screen, double **host_cut_ljsq,
                 const double host_cut_coulsq, double *host_special_coul, const double qqrd2e,
                 const double e_shift, const double f_shift, const double alpha);
void ljd_gpu_clear();
int **ljd_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                        int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                        tagint **special, const bool eflag, const bool vflag, const bool eatom,
                        const bool vatom, int **ilist, int **jnum,
                        bool &success, double *host_q, double *boxlo,
                        double *prd, int *periodicity);
void ljd_gpu_compute(const int ago, const int inum, const int nall, double **host_x, int *host_type,
                     int *ilist, int *numj, int **firstneigh, const bool eflag, const bool vflag,
                     const bool eatom, const bool vatom, bool &success, double *host_q, const int nlocal, double *boxlo, double *prd);
double ljd_gpu_bytes();

int ljcl_gpu_init(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                  double **host_lj3, double **host_lj4, double **offset, double *special_lj,
                  const int nlocal, const int nall, const int max_nbors, const int maxspecial,
                  const double cell_size, int &gpu_mode, FILE *screen, double **host_cut_ljsq,
                  double host_cut_coulsq, double *host_special_coul, const double qqrd2e,
                  const double g_ewald);
void ljcl_gpu_reinit(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                     double **host_lj3, double **host_lj4, double **offset, double **host_lj_cutsq);
void ljcl_gpu_clear();
int **ljcl_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                         int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                         tagint **special, const bool eflag, const bool vflag, const bool eatom,
                         const bool vatom, int **ilist, int **jnum,
                         bool &success, double *host_q, double *boxlo,
                         double *prd, int* periodicity);
void ljcl_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                      int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                      const bool vflag, const bool eatom, const bool vatom, bool &success, double *host_q, const int nlocal,
                      double *boxlo, double *prd);
double ljcl_gpu_bytes();

int ljcls_gpu_init(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                  double **host_lj3, double **host_lj4, double **offset,  double **epsilon, double *special_lj,
                  const int nlocal, const int nall, const int max_nbors, const int maxspecial,
                  const double cell_size, int &gpu_mode, FILE *screen, double **host_cut_ljsq,
                  double host_cut_coulsq, double *host_special_coul, const double qqrd2e,
                  const double g_ewald);
void ljcls_gpu_reinit(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                     double **host_lj3, double **host_lj4, double **offset, double **epsilon,
                     double **host_lj_cutsq);
void ljcls_gpu_clear();
int **ljcls_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                         int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                         tagint **special, const bool eflag, const bool vflag, const bool eatom,
                         const bool vatom, int **ilist, int **jnum,
                         bool &success, double *host_q, double *boxlo,
                         double *prd, int *periodicity);
void ljcls_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                      int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                      const bool vflag, const bool eatom, const bool vatom, bool &success, double *host_q, const int nlocal,
                      double *boxlo, double *prd);
double ljcls_gpu_bytes();

int ljcm_gpu_init(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                  double **host_lj3, double **host_lj4, double **host_gcons, double **host_dgcons,
                  double **offset, double *special_lj, const int inum, const int nall,
                  const int max_nbors, const int maxspecial, const double cell_size, int &gpu_mode,
                  FILE *screen, double **host_cut_ljsq, double host_cut_coulsq,
                  double *host_special_coul, const int order, const double qqrd2e);
void ljcm_gpu_clear();
int **ljcm_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                         int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                         tagint **special, const bool eflag, const bool vflag, const bool eatom,
                         const bool vatom, int **ilist, int **jnum,
                         bool &success, double *host_q, double *boxlo,
                         double *prd, int *periodicity);
void ljcm_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                      int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                      const bool vflag, const bool eatom, const bool vatom, bool &success, double *host_q, const int nlocal,
                      double *boxlo, double *prd);
double ljcm_gpu_bytes();

int dpl_gpu_init(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                 double **host_lj3, double **host_lj4, double **offset, double *special_lj,
                 const int nlocal, const int nall, const int max_nbors, const int maxspecial,
                 const double cell_size, int &gpu_mode, FILE *screen, double **host_cut_ljsq,
                 double **host_cut_coulsq, double *host_special_coul, const double qqrd2e);
void dpl_gpu_clear();
int **dpl_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                        int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                        tagint **special, const bool eflag, const bool vflag, const bool eatom,
                        const bool vatom, int **ilist, int **jnum,
                        bool &success, double *host_q, double **host_mu,
                        double *boxlo, double *prd);
void dpl_gpu_compute(const int ago, const int inum, const int nall, double **host_x, int *host_type,
                     int *ilist, int *numj, int **firstneigh, const bool eflag, const bool vflag,
                     const bool eatom, const bool vatom, bool &success, double *host_q, double **host_mu, const int nlocal,
                     double *boxlo, double *prd);
double dpl_gpu_bytes();

int dplj_gpu_init(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                  double **host_lj3, double **host_lj4, double **offset, double *special_lj,
                  const int nlocal, const int nall, const int max_nbors, const int maxspecial,
                  const double cell_size, int &gpu_mode, FILE *screen, double **host_cut_ljsq,
                  const double host_cut_coulsq, double *host_special_coul, const double qqrd2e,
                  const double g_ewald);
void dplj_gpu_clear();
int **dplj_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                         int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                         tagint **special, const bool eflag, const bool vflag, const bool eatom,
                         const bool vatom, int **ilist, int **jnum,
                         bool &success, double *host_q, double **host_mu,
                         double *boxlo, double *prd);
void dplj_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                      int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                      const bool vflag, const bool eatom, const bool vatom, bool &success, double *host_q, double **host_mu,
                      const int nlocal, double *boxlo, double *prd);
double dplj_gpu_bytes();

int ljl_gpu_init(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                 double **host_lj3, double **host_lj4, double **offset, double *special_lj,
                 const int nlocal, const int nall, const int max_nbors, const int maxspecial,
                 const double cell_size, int &gpu_mode, FILE *screen);
void ljl_gpu_reinit(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                    double **host_lj3, double **host_lj4, double **offset);
void ljl_gpu_clear();
int **ljl_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                        int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                        tagint **special, const bool eflag, const bool vflag, const bool eatom,
                        const bool vatom, int **ilist, int **jnum,
                        bool &success, double *prd, int *periodicity);
void ljl_gpu_compute(const int ago, const int inum, const int nall, double **host_x, int *host_type,
                     int *ilist, int *numj, int **firstneigh, const bool eflag, const bool vflag,
                     const bool eatom, const bool vatom, bool &success);
double ljl_gpu_bytes();

int ljtip4p_long_gpu_init(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                          double **host_lj3, double **host_lj4, double **offset, double *special_lj,
                          const int nlocal, const int tH, const int tO, const double alpha,
                          const double qdist, const int nall, const int max_nbors,
                          const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen,
                          double **host_cut_ljsq, const double host_cut_coulsq,
                          const double host_cut_coulsqplus, double *host_special_coul,
                          const double qqrd2e, const double g_ewald, int map_size, int max_same);
void ljtip4p_long_gpu_clear();
int **ljtip4p_long_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                                 int *host_type, double *sublo, double *subhi, tagint *tag,
                                 int *map_array, int map_size, int *sametag, int max_same,
                                 int **nspecial, tagint **special, const bool eflag,
                                 const bool vflag, const bool eatom, const bool vatom,
                                 int **ilist, int **jnum, bool &success, double *host_q, double *boxlo, double *prd,
                                 int *periodicity);
void ljtip4p_long_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                              int *host_type, int *ilist, int *numj, int **firstneigh,
                              const bool eflag, const bool vflag, const bool eatom,
                              const bool vatom, bool &success, double *host_q, const int nlocal, double *boxlo,
                              double *prd);
double ljtip4p_long_gpu_bytes();
void ljtip4p_long_copy_molecule_data(int, tagint *, int *, int, int *, int, int);

int ljecl_gpu_init(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                   double **host_lj3, double **host_lj4, double **offset, double **shift,
                   double *special_lj, const int nlocal, const int nall, const int max_nbors,
                   const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen,
                   double **host_cut_ljsq, double host_cut_coulsq, double *host_special_coul,
                   const double qqrd2e, const double g_ewald);
void ljecl_gpu_reinit(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                      double **host_lj3, double **host_lj4, double **offset, double **shift,
                      double **host_lj_cutsq);
void ljecl_gpu_clear();
int **ljecl_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                          int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                          tagint **special, const bool eflag, const bool vflag, const bool eatom,
                          const bool vatom, int **ilist, int **jnum,
                          bool &success, double *host_q, double *boxlo,
                          double *prd, int *periodicity);
void ljecl_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                       int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                       const bool vflag, const bool eatom, const bool vatom, bool &success, double *host_q, const int nlocal,
                       double *boxlo, double *prd);
double ljecl_gpu_bytes();

int lje_gpu_init(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                 double **host_lj3, double **host_lj4, double **offset, double **shift,
                 double *special_lj, const int nlocal, const int nall, const int max_nbors,
                 const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen);
void lje_gpu_reinit(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                    double **host_lj3, double **host_lj4, double **offset, double **shift);
void lje_gpu_clear();
int **lje_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                        int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                        tagint **special, const bool eflag, const bool vflag, const bool eatom,
                        const bool vatom, int **ilist, int **jnum,
                        bool &success, double *prd, int *periodicity);
void lje_gpu_compute(const int ago, const int inum, const int nall, double **host_x, int *host_type,
                     int *ilist, int *numj, int **firstneigh, const bool eflag, const bool vflag,
                     const bool eatom, const bool vatom, bool &success);
double lje_gpu_bytes();

int ljgrm_gpu_init(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                   double **host_lj3, double **host_lj4, double *special_lj, const int inum,
                   const int nall, const int max_nbors, const int maxspecial,
                   const double cell_size, int &gpu_mode, FILE *screen, double **host_ljsw1,
                   double **host_ljsw2, double **host_ljsw3, double **host_ljsw4,
                   double **host_ljsw5, double **cut_inner, double **cut_innersq);
void ljgrm_gpu_clear();
int **ljgrm_gpu_compute_n(const int ago, const int inum_full, const int nall, double **host_x,
                          int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                          tagint **special, const bool eflag, const bool vflag, const bool eatom,
                          const bool vatom, int **ilist, int **jnum,
                          bool &success, double *prd, int *periodicity);
void ljgrm_gpu_compute(const int ago, const int inum_full, const int nall, double **host_x,
                       int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                       const bool vflag, const bool eatom, const bool vatom, bool &success);
double ljgrm_gpu_bytes();

int dplsf_gpu_init(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                   double **host_lj3, double **host_lj4, double *special_lj, const int nlocal,
                   const int nall, const int max_nbors, const int maxspecial,
                   const double cell_size, int &gpu_mode, FILE *screen, double **host_cut_ljsq,
                   double **host_cut_coulsq, double *host_special_coul, const double qqrd2e);
void dplsf_gpu_clear();
int **dplsf_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                          int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                          tagint **special, const bool eflag, const bool vflag, const bool eatom,
                          const bool vatom, int **ilist, int **jnum,
                          bool &success, double *host_q, double **host_mu,
                          double *boxlo, double *prd);
void dplsf_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                       int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                       const bool vflag, const bool eatom, const bool vatom, bool &success, double *host_q, double **host_mu,
                       const int nlocal, double *boxlo, double *prd);
double dplsf_gpu_bytes();

int ljsmt_gpu_init(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                   double **host_lj3, double **host_lj4, double **offset, double *special_lj,
                   const int nlocal, const int nall, const int max_nbors, const int maxspecial,
                   const double cell_size, int &gpu_mode, FILE *screen, double **host_ljsw0,
                   double **host_ljsw1, double **host_ljsw2, double **host_ljsw3,
                   double **host_ljsw4, double **cut_inner, double **cut_innersq);
void ljsmt_gpu_reinit(const int ntypes, double **cutsq, double **host_lj1, double **host_lj2,
                      double **host_lj3, double **host_lj4, double **offset, double **host_ljsw0,
                      double **host_ljsw1, double **host_ljsw2, double **host_ljsw3,
                      double **host_ljsw4, double **cut_inner, double **cut_innersq);
void ljsmt_gpu_clear();
int **ljsmt_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                          int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                          tagint **special, const bool eflag, const bool vflag, const bool eatom,
                          const bool vatom, int **ilist, int **jnum,
                          bool &success, double *prd, int *periodicity);
void ljsmt_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                       int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                       const bool vflag, const bool eatom, const bool vatom, bool &success);
double ljsmt_gpu_bytes();

int spical_gpu_init(const int ntypes, double **cutsq, int **lj_type, double **host_lj1,
                    double **host_lj2, double **host_lj3, double **host_lj4, double **offset,
                    double *special_lj, const int nlocal, const int nall, const int max_nbors,
                    const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen,
                    double **host_cut_ljsq, double host_cut_coulsq, double *host_special_coul,
                    const double qqrd2e, const double g_ewald);
void spical_gpu_clear();
int **spical_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                           int *host_type, double *sublo, double *subhi, tagint *tag,
                           int **nspecial, tagint **special, const bool eflag, const bool vflag,
                           const bool eatom, const bool vatom, int **ilist,
                           int **jnum, bool &success, double *host_q,
                           double *boxlo, double *prd, int *periodicity);
void spical_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                        int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                        const bool vflag, const bool eatom, const bool vatom, bool &success, double *host_q, const int nlocal,
                        double *boxlo, double *prd);
double spical_gpu_bytes();
#include "lj_spica_common.h"
using namespace LJSPICAParms;

int spica_gpu_init(const int ntypes, double **cutsq, int **cg_types, double **host_lj1,
                   double **host_lj2, double **host_lj3, double **host_lj4, double **offset,
                   double *special_lj, const int nlocal, const int nall, const int max_nbors,
                   const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen);
void spica_gpu_clear();
int **spica_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                          int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                          tagint **special, const bool eflag, const bool vflag, const bool eatom,
                          const bool vatom, int **ilist, int **jnum,
                          bool &success, double *prd, int *periodicity);
void spica_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                       int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                       const bool vflag, const bool eatom, const bool vatom, bool &success);
double spica_gpu_bytes();
#include "lj_spica_common.h"
using namespace LJSPICAParms;

int mdpd_gpu_init(const int ntypes, double **cutsq, double **host_A_att, double **host_B_rep,
                  double **host_gamma, double **host_sigma, double **host_cut, double **host_cut_r,
                  double *special_lj, const int inum, const int nall, const int max_nbors,
                  const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen);
void mdpd_gpu_clear();
int **mdpd_gpu_compute_n(const int ago, const int inum_full, const int nall, double **host_x,
                         int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                         tagint **special, const bool eflag, const bool vflag, const bool eatom,
                         const bool vatom, int **ilist, int **jnum,
                         bool &success, double **host_v,
                         const double dtinvsqrt, const int seed, const int timestep, double *boxlo,
                         double *prd);
void mdpd_gpu_compute(const int ago, const int inum_full, const int nall, double **host_x,
                      int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                      const bool vflag, const bool eatom, const bool vatom, bool &success, tagint *tag, double **host_v,
                      const double dtinvsqrt, const int seed, const int timestep, const int nlocal,
                      double *boxlo, double *prd);
void mdpd_gpu_get_extra_data(double *host_rho);
double mdpd_gpu_bytes();

int mie_gpu_init(const int ntypes, double **cutsq, double **host_mie1, double **host_mie2,
                 double **host_mie3, double **host_mie4, double **host_gamA, double **host_gamR,
                 double **offset, double *special_lj, const int nlocal, const int nall,
                 const int max_nbors, const int maxspecial, const double cell_size, int &gpu_mode,
                 FILE *screen);
void mie_gpu_clear();
int **mie_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                        int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                        tagint **special, const bool eflag, const bool vflag, const bool eatom,
                        const bool vatom, int **ilist, int **jnum,
                        bool &success, double *prd, int *periodicity);
void mie_gpu_compute(const int ago, const int inum, const int nall, double **host_x, int *host_type,
                     int *ilist, int *numj, int **firstneigh, const bool eflag, const bool vflag,
                     const bool eatom, const bool vatom, bool &success);
double mie_gpu_bytes();

int mor_gpu_init(const int ntypes, double **cutsq, double **host_morse1, double **host_r0,
                 double **host_alpha, double **host_d0, double **offset, double *special_lj,
                 const int nlocal, const int nall, const int max_nbors, const int maxspecial,
                 const double cell_size, int &gpu_mode, FILE *screen);
void mor_gpu_clear();
int **mor_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                        int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                        tagint **special, const bool eflag, const bool vflag, const bool eatom,
                        const bool vatom, int **ilist, int **jnum,
                        bool &success, double *prd, int *periodicity);
void mor_gpu_compute(const int ago, const int inum, const int nall, double **host_x, int *host_type,
                     int *ilist, int *numj, int **firstneigh, const bool eflag, const bool vflag,
                     const bool eatom, const bool vatom, bool &success);
double mor_gpu_bytes();

int re_gpu_init(const int ntypes, double **shape, double **well, double **cutsq,
                double **sigma, double **epsilon, int **form,
                double **host_lj1, double **host_lj2, double **host_lj3,
                double **host_lj4, double **offset, double *special_lj,
                const int nlocal, const int nall, const int max_nbors,
                const int maxspecial, const double cell_size, int &gpu_mode,
                FILE *screen);
void re_gpu_clear();
int **re_gpu_compute_n(const int ago, const int inum, const int nall,
                       double **host_x, int *host_type, double *sublo,
                       double *subhi, tagint *tag, int **nspecial,
                       tagint **special, const bool eflag, const bool vflag,
                       const bool eatom, const bool vatom, int **ilist, int **jnum, bool &success, const int *ellipsoid,
                       const void *bonus);
int *re_gpu_compute(const int ago, const int inum, const int nall,
                    double **host_x, int *host_type, int *ilist, int *numj,
                    int **firstneigh, const bool eflag, const bool vflag,
                    const bool eatom, const bool vatom, bool &success, const int *ellipsoid,
                    const void *bonus);
double re_gpu_bytes();
enum { SPHERE_SPHERE, SPHERE_ELLIPSE, ELLIPSE_SPHERE, ELLIPSE_ELLIPSE };

int soft_gpu_init(const int ntypes, double **cutsq, double **prefactor, double **cut,
                  double *special_lj, const int nlocal, const int nall, const int max_nbors,
                  const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen);
void soft_gpu_reinit(const int ntypes, double **cutsq, double **host_prefactor, double **host_cut);
void soft_gpu_clear();
int **soft_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                         int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                         tagint **special, const bool eflag, const bool vflag, const bool eatom,
                         const bool vatom, int **ilist, int **jnum,
                         bool &success, double *prd, int *periodicity);
void soft_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                      int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                      const bool vflag, const bool eatom, const bool vatom, bool &success);
double soft_gpu_bytes();
using namespace MathConst;

int sph_heatconduction_gpu_init(const int ntypes, double **cutsq, double** host_cut,
                                double **host_alpha, double* host_mass,
                                const int dimension, double *special_lj,
                                const int inum, const int nall,
                                const int max_nbors,  const int maxspecial,
                                const double cell_size, int &gpu_mode, FILE *screen);
void sph_heatconduction_gpu_clear();
int **sph_heatconduction_gpu_compute_n(const int ago, const int inum_full, const int nall,
                           double **host_x, int *host_type, double *sublo,
                           double *subhi, tagint *host_tag, int **nspecial,
                           tagint **special, const bool eflag, const bool vflag,
                           const bool eatom, const bool vatom, int **ilist, int **jnum, bool &success,
                           double **host_v);
void sph_heatconduction_gpu_compute(const int ago, const int inum_full, const int nall,
                        double **host_x, int *host_type, int *ilist, int *numj,
                        int **firstneigh, const bool eflag, const bool vflag,
                        const bool eatom, const bool vatom, bool &success, tagint *host_tag,
                        double **host_v);
void sph_heatconduction_gpu_get_extra_data(double *host_rho, double *host_esph);
void sph_heatconduction_gpu_update_dE(void **dE_ptr);
double sph_heatconduction_gpu_bytes();

int sph_lj_gpu_init(const int ntypes, double **cutsq, double** host_cut,
                    double **host_viscosity, double* host_mass,
                     const int dimension, double *special_lj,
                    const int inum, const int nall,
                    const int max_nbors,  const int maxspecial,
                    const double cell_size, int &gpu_mode, FILE *screen);
void sph_lj_gpu_clear();
int **sph_lj_gpu_compute_n(const int ago, const int inum_full, const int nall,
                           double **host_x, int *host_type, double *sublo,
                           double *subhi, tagint *host_tag, int **nspecial,
                           tagint **special, const bool eflag, const bool vflag,
                           const bool eatom, const bool vatom, int **ilist, int **jnum, bool &success,
                           double **host_v);
void sph_lj_gpu_compute(const int ago, const int inum_full, const int nall,
                        double **host_x, int *host_type, int *ilist, int *numj,
                        int **firstneigh, const bool eflag, const bool vflag,
                        const bool eatom, const bool vatom, bool &success, tagint *host_tag,
                        double **host_v);
void sph_lj_gpu_get_extra_data(double *host_rho, double *host_esph,
                               double *host_cv);
void sph_lj_gpu_update_drhoE(void **drhoE_ptr);
double sph_lj_gpu_bytes();

int sph_taitwater_gpu_init(const int ntypes, double **cutsq, double** host_cut,
                           double **host_viscosity, double* host_mass, double* host_rho0,
                           double* host_soundspeed, double* host_B, const int dimension,
                           double *special_lj, const int inum, const int nall,
                           const int max_nbors,  const int maxspecial,
                           const double cell_size, int &gpu_mode, FILE *screen);
void sph_taitwater_gpu_clear();
int **sph_taitwater_gpu_compute_n(const int ago, const int inum_full, const int nall,
                         double **host_x, int *host_type, double *sublo,
                         double *subhi, tagint *tag, int **nspecial,
                         tagint **special, const bool eflag, const bool vflag,
                         const bool eatom, const bool vatom, int **ilist, int **jnum, bool &success,
                         double **host_v);
void sph_taitwater_gpu_compute(const int ago, const int inum_full, const int nall,
                        double **host_x, int *host_type, int *ilist, int *numj,
                        int **firstneigh, const bool eflag, const bool vflag,
                        const bool eatom, const bool vatom, bool &success, tagint *tag,
                        double **host_v);
void sph_taitwater_gpu_get_extra_data(double *host_rho);
void sph_taitwater_gpu_update_drhoE(void **drhoE_ptr);
double sph_taitwater_gpu_bytes();

int sw_gpu_init(const int ntypes, const int inum, const int nall, const int max_nbors,
                const double cell_size, int &gpu_mode, FILE *screen, double **ncutsq, double **ncut,
                double **sigma, double **powerp, double **powerq, double **sigma_gamma, double **c1,
                double **c2, double **c3, double **c4, double **c5, double **c6,
                double ***lambda_epsilon, double ***costheta, const int *map, int ***e2param);
void sw_gpu_clear();
int **sw_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                       int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                       tagint **special, const bool eflag, const bool vflag, const bool eatom,
                       const bool vatom, int **ilist, int **jnum,
                       bool &success);
void sw_gpu_compute(const int ago, const int nloc, const int nall, const int ln, double **host_x,
                    int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                    const bool vflag, const bool eatom, const bool vatom, bool &success);
double sw_gpu_bytes();

int table_gpu_init(const int ntypes, double **cutsq, double ***host_table_coeffs,
                   double **host_table_data, double *special_lj, const int nlocal, const int nall,
                   const int max_nbors, const int maxspecial, const double cell_size, int &gpu_mode,
                   FILE *screen, int tabstyle, int ntables, int tablength);
void table_gpu_clear();
int **table_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                          int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                          tagint **special, const bool eflag, const bool vflag, const bool eatom,
                          const bool vatom, int **ilist, int **jnum,
                          bool &success, double *prd, int *periodicity);
void table_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                       int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                       const bool vflag, const bool eatom, const bool vatom, bool &success);
double table_gpu_bytes();

int tersoff_gpu_init(const int ntypes, const int inum, const int nall, const int max_nbors,
                     const double cell_size, int &gpu_mode, FILE *screen, int *host_map,
                     const int nelements, int ***host_elem3param, const int nparams,
                     const double *ts_lam1, const double *ts_lam2, const double *ts_lam3,
                     const double *ts_powermint, const double *ts_biga, const double *ts_bigb,
                     const double *ts_bigr, const double *ts_bigd, const double *ts_c1,
                     const double *ts_c2, const double *ts_c3, const double *ts_c4,
                     const double *ts_c, const double *ts_d, const double *ts_h,
                     const double *ts_gamma, const double *ts_beta, const double *ts_powern,
                     const double *ts_cutsq);
void tersoff_gpu_clear();
int **tersoff_gpu_compute_n(const int ago, const int inum_full, const int nall, double **host_x,
                            int *host_type, double *sublo, double *subhi, tagint *tag,
                            int **nspecial, tagint **special, const bool eflag, const bool vflag,
                            const bool eatom, const bool vatom, int **ilist,
                            int **jnum, bool &success);
void tersoff_gpu_compute(const int ago, const int nlocal, const int nall, const int nlist,
                         double **host_x, int *host_type, int *ilist, int *numj, int **firstneigh,
                         const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                         bool &success);
double tersoff_gpu_bytes();

int tersoff_mod_gpu_init(const int ntypes, const int inum, const int nall, const int max_nbors,
                         const double cell_size, int &gpu_mode, FILE *screen, int *host_map,
                         const int nelements, int ***host_elem3param, const int nparams,
                         const double *ts_lam1, const double *ts_lam2, const double *ts_lam3,
                         const double *ts_powermint, const double *ts_biga, const double *ts_bigb,
                         const double *ts_bigr, const double *ts_bigd, const double *ts_c1,
                         const double *ts_c2, const double *ts_c3, const double *ts_c4,
                         const double *ts_c5, const double *ts_h, const double *ts_beta,
                         const double *ts_powern, const double *ts_powern_del, const double *ts_ca1,
                         const double *ts_cutsq);
void tersoff_mod_gpu_clear();
int **tersoff_mod_gpu_compute_n(const int ago, const int inum_full, const int nall, double **host_x,
                                int *host_type, double *sublo, double *subhi, tagint *tag,
                                int **nspecial, tagint **special, const bool eflag,
                                const bool vflag, const bool eatom, const bool vatom,
                                int **ilist, int **jnum, bool &success);
void tersoff_mod_gpu_compute(const int ago, const int nlocal, const int nall, const int nlist,
                             double **host_x, int *host_type, int *ilist, int *numj,
                             int **firstneigh, const bool eflag, const bool vflag, const bool eatom,
                             const bool vatom, bool &success);
double tersoff_mod_gpu_bytes();

int tersoff_zbl_gpu_init(const int ntypes, const int inum, const int nall, const int max_nbors,
                         const double cell_size, int &gpu_mode, FILE *screen, int *host_map,
                         const int nelements, int ***host_elem3param, const int nparams,
                         const double *ts_lam1, const double *ts_lam2, const double *ts_lam3,
                         const double *ts_powermint, const double *ts_biga, const double *ts_bigb,
                         const double *ts_bigr, const double *ts_bigd, const double *ts_c1,
                         const double *ts_c2, const double *ts_c3, const double *ts_c4,
                         const double *ts_c, const double *ts_d, const double *ts_h,
                         const double *ts_gamma, const double *ts_beta, const double *ts_powern,
                         const double *ts_Z_i, const double *ts_Z_j, const double *ts_ZBLcut,
                         const double *ts_ZBLexpscale, const double global_e,
                         const double global_a_0, const double global_epsilon_0,
                         const double *ts_cutsq);
void tersoff_zbl_gpu_clear();
int **tersoff_zbl_gpu_compute_n(const int ago, const int inum_full, const int nall, double **host_x,
                                int *host_type, double *sublo, double *subhi, tagint *tag,
                                int **nspecial, tagint **special, const bool eflag,
                                const bool vflag, const bool eatom, const bool vatom,
                                int **ilist, int **jnum, bool &success);
void tersoff_zbl_gpu_compute(const int ago, const int nlocal, const int nall, const int nlist,
                             double **host_x, int *host_type, int *ilist, int *numj,
                             int **firstneigh, const bool eflag, const bool vflag, const bool eatom,
                             const bool vatom, bool &success);
double tersoff_zbl_gpu_bytes();

int ufml_gpu_init(const int ntypes, double **cutsq, double **host_uf1, double **host_uf2,
                  double **host_uf3, double **offset, double *special_lj, const int nlocal,
                  const int nall, const int max_nbors, const int maxspecial, const double cell_size,
                  int &gpu_mode, FILE *screen);
void ufml_gpu_reinit(const int ntypes, double **cutsq, double **host_uf1, double **host_uf2,
                     double **host_uf3, double **offset);
void ufml_gpu_clear();
int **ufml_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                         int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                         tagint **special, const bool eflag, const bool vflag, const bool eatom,
                         const bool vatom, int **ilist, int **jnum,
                         bool &success, double *prd, int *periodicity);
void ufml_gpu_compute(const int ago, const int inum, const int nall, double **host_x,
                      int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                      const bool vflag, const bool eatom, const bool vatom, bool &success);
double ufml_gpu_bytes();

int vashishta_gpu_init(const int ntypes, const int inum, const int nall, const int max_nbors,
                       const double cell_size, int &gpu_mode, FILE *screen, int *host_map,
                       const int nelements, int ***host_elem3param, const int nparams,
                       const double *cutsq, const double *r0, const double *gamma,
                       const double *eta, const double *lam1inv, const double *lam4inv,
                       const double *zizj, const double *mbigd, const double *dvrc,
                       const double *big6w, const double *heta, const double *bigh,
                       const double *bigw, const double *c0, const double *costheta,
                       const double *bigb, const double *big2b, const double *bigc);
void vashishta_gpu_clear();
int **vashishta_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                              int *host_type, double *sublo, double *subhi, tagint *tag,
                              int **nspecial, tagint **special, const bool eflag, const bool vflag,
                              const bool eatom, const bool vatom, int **ilist,
                              int **jnum, bool &success);
void vashishta_gpu_compute(const int ago, const int nloc, const int nall, const int ln,
                           double **host_x, int *host_type, int *ilist, int *numj, int **firstneigh,
                           const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                           bool &success);
double vashishta_gpu_bytes();

int ykcolloid_gpu_init(const int ntypes, double **cutsq, double **host_a, double **host_offset,
                       double *special_lj, const int inum, const int nall, const int max_nbors,
                       const int maxspecial, const double cell_size, int &gpu_mode, FILE *screen,
                       const double kappa);
void ykcolloid_gpu_clear();
int **ykcolloid_gpu_compute_n(const int ago, const int inum_full, const int nall, double **host_x,
                              int *host_type, double *sublo, double *subhi, tagint *tag,
                              int **nspecial, tagint **special, const bool eflag, const bool vflag,
                              const bool eatom, const bool vatom, int **ilist,
                              int **jnum, bool &success, double *host_rad, double *prd,
                              int *periodicity);
void ykcolloid_gpu_compute(const int ago, const int inum_full, const int nall, double **host_x,
                           int *host_type, int *ilist, int *numj, int **firstneigh,
                           const bool eflag, const bool vflag, const bool eatom, const bool vatom,
                           bool &success, double *host_rad);
double ykcolloid_gpu_bytes();

int yukawa_gpu_init(const int ntypes, double **cutsq, double kappa, double **host_a,
                    double **offset, double *special_lj, const int inum, const int nall,
                    const int max_nbors, const int maxspecial, const double cell_size,
                    int &gpu_mode, FILE *screen);
void yukawa_gpu_clear();
int **yukawa_gpu_compute_n(const int ago, const int inum_full, const int nall, double **host_x,
                           int *host_type, double *sublo, double *subhi, tagint *tag,
                           int **nspecial, tagint **special, const bool eflag, const bool vflag,
                           const bool eatom, const bool vatom, int **ilist,
                           int **jnum, bool &success, double *prd, int *periodicity);
void yukawa_gpu_compute(const int ago, const int inum_full, const int nall, double **host_x,
                        int *host_type, int *ilist, int *numj, int **firstneigh, const bool eflag,
                        const bool vflag, const bool eatom, const bool vatom, bool &success);
double yukawa_gpu_bytes();

int zbl_gpu_init(const int ntypes, double **cutsq, double **host_sw1, double **host_sw2,
                 double **host_sw3, double **host_sw4, double **host_sw5, double **host_d1a,
                 double **host_d2a, double **host_d3a, double **host_d4a, double **host_zze,
                 double cut_globalsq, double cut_innersq, double cut_inner, const int inum,
                 const int nall, const int max_nbors, const int maxspecial, const double cell_size,
                 int &gpu_mode, FILE *screen);
void zbl_gpu_clear();
int **zbl_gpu_compute_n(const int ago, const int inum, const int nall, double **host_x,
                        int *host_type, double *sublo, double *subhi, tagint *tag, int **nspecial,
                        tagint **special, const bool eflag, const bool vflag, const bool eatom,
                        const bool vatom, int **ilist, int **jnum,
                        bool &success, double *prd, int *periodicity);
void zbl_gpu_compute(const int ago, const int inum, const int nall, double **host_x, int *host_type,
                     int *ilist, int *numj, int **firstneigh, const bool eflag, const bool vflag,
                     const bool eatom, const bool vatom, bool &success);
double zbl_gpu_bytes();

} // namespace LAMMPS_GPU

#endif // LMP_LAMMPS_GPU_H
