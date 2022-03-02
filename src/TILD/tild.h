/* ------------------------------------------------------------------
    TILD - Theoretically Informed Langevan Dynamics
    This replicates the TILD coe done by the Riggleman group, 
    previously known as Dynamical Mean Field Theory. 
    
    Copyright (2019) Christian Tabedzki and Zachariah Vicars.
    tabedzki@seas.upenn.edu zvicars@seas.upenn.edu
-------------------------------------------------------------------- */

#ifdef KSPACE_CLASS
// clang-format off
KSpaceStyle(tild,TILD)
// clang-format on
#else

#ifndef LMP_TILD_H
#define LMP_TILD_H

#if defined(FFT_FFTW3)
#define LMP_FFT_LIB "FFTW3"
#elif defined(FFT_MKL)
#define LMP_FFT_LIB "MKL FFT"
#elif defined(FFT_CUFFT)
#define LMP_FFT_LIB "cuFFT"
#else
#define LMP_FFT_LIB "KISS FFT"
#endif

#ifdef FFT_SINGLE
typedef float FFT_SCALAR;
#define LMP_FFT_PREC "single"
#define MPI_FFT_SCALAR MPI_FLOAT
#else

typedef double FFT_SCALAR;
#define LMP_FFT_PREC "double"
#define MPI_FFT_SCALAR MPI_DOUBLE
#endif

#include "kspace.h"
#include <algorithm>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace LAMMPS_NS {

class Interaction;

class TILD : public KSpace {
 public:
  TILD(class LAMMPS *);
  virtual ~TILD();
  virtual void init();
  virtual void setup();
  virtual void settings(int, char **);
  virtual void compute(int, int);
  double memory_usage();
  void setup_grid();

 protected:
  // For future kspace_hybrid
  int nstyles;      // # of sub-styles For future kspace_hybrid
  int **setflag;    // 0/1 = whether each i,j has been set
  int *density_flags;

  int nfactors;
  int *factors;

  double *fkx, *fky, *fkz;
  double *fkx2, *fky2, *fkz2;

  FFT_SCALAR ***vg;
  FFT_SCALAR ***vg_hat;

  FFT_SCALAR *work1, *work2;

  FFT_SCALAR **rho1d, **rho_coeff, **drho_coeff;
  FFT_SCALAR ***grad_potent, ***grad_potent_hat, **potent, **potent_hat;
  FFT_SCALAR ****gradWtypex, ****gradWtypey, ****gradWtypez;

  FFT_SCALAR ****density_brick_types;
  FFT_SCALAR ****avg_density_brick_types;
  FFT_SCALAR *kappa_density;
  FFT_SCALAR **density_fft_types;
  FFT_SCALAR **density_hat_fft_types;
  FFT_SCALAR *ktmp;
  FFT_SCALAR *ktmpi;
  FFT_SCALAR *ktmpj;
  FFT_SCALAR *ktmp2;
  FFT_SCALAR *ktmp2i;
  FFT_SCALAR *ktmp2j;
  FFT_SCALAR *tmp;

  int **potent_map;
  double volume;
  int nmax;

  void precompute_density_hat_fft();
  void complex_multiply(FFT_SCALAR *, FFT_SCALAR *, FFT_SCALAR *, int);
  double **potent_param;
  int npot, *pot_map;
  int ***potent_type_map;
  double rho0, set_rho0, modified_rho0;

  int factorable(int);
  double **chi;
  double **a2;
  double **rp;
  double **xi;
  double **np_rho;
  double gridsize;
  virtual int modify_param(int, char **);
  int num_potent;
  double calculate_rho0();
  int get_style(const int, const int);

  // group-group interactions

  int group_allocate_flag;

  double rms(int, double, bigint, double);
  virtual void allocate();
  void deallocate();
  void init_potential(FFT_SCALAR *, const int, const double *);
  void init_potential_ft(FFT_SCALAR *, const int, const double *);
  void calc_work(FFT_SCALAR *, const int, const int);
  void calc_cross_work(const Interaction &);
  void init_cross_potentials();
  double get_k_alias(int, double *);
  void get_k_alias(FFT_SCALAR *, FFT_SCALAR **);
  void manually_flip_density_flags();

  int me, nprocs;
  double cutoff;
  double kappa;
  double delxinv, delyinv, delzinv, delvolinv;
  double h_x, h_y, h_z;
  double shift, shiftone;
  int peratom_allocate_flag;

  int nxlo_in, nylo_in, nzlo_in, nxhi_in, nyhi_in, nzhi_in;
  int nxlo_out, nylo_out, nzlo_out, nxhi_out, nyhi_out, nzhi_out;
  int nxlo_ghost, nxhi_ghost, nylo_ghost, nyhi_ghost, nzlo_ghost, nzhi_ghost;
  int nxlo_fft, nylo_fft, nzlo_fft, nxhi_fft, nyhi_fft, nzhi_fft;
  int nlower, nupper;
  int ngrid, nfft, nfft_both;
  int subtract_rho0, normalize_by_rho0, mix_flag, sub_flag, norm_flag;

  std::vector<Interaction> cross_iter;

  int write_grid_flag, grid_data_output_freq;
  int ave_grid_flag, nvalid_last, nvalid, nevery, irepeat, nrepeat, peratom_freq;
  char *grid_data_filename = new char[50];
  char *ave_grid_filename = new char[50];
  FILE *otp;

  // group-group interactions

  class FFT3d *fft1, *fft2;
  class Remap *remap;
  class GridComm *gc;

  FFT_SCALAR *gc_buf1, *gc_buf2;
  int ngc_buf1, ngc_buf2, npergrid;

  int **part2grid;    // storage for particle -> grid mapping

  double *boxlo;

  virtual void set_grid_global();
  void set_grid_local();

  virtual void particle_map();
  virtual void make_rho();
  void brick2fft();
  void ev_calculation(int, int, int);
  ;

  virtual void fieldforce_param();

  void procs2grid2d(int, int, int, int *, int *);
  void compute_rho1d(const FFT_SCALAR &, const FFT_SCALAR &, const FFT_SCALAR &, int, FFT_SCALAR **,
                     FFT_SCALAR **);
  void compute_rho_coeff();

  // grid communication

  virtual void pack_forward_grid(int, void *, int, int *);
  virtual void unpack_forward_grid(int, void *, int, int *);
  virtual void pack_reverse_grid(int, void *, int, int *);
  virtual void unpack_reverse_grid(int, void *, int, int *);

  // triclinic

  int triclinic;    // domain settings, orthog or triclinic

  // group-group interactions

  void accumulate_gradient();
  void vir_func_init();
  void write_grid_data(char *, const int);
  void pack_grid_data(double **);
  void pack_avg_grid_data(double **buf);
  void sum_grid_data();
  void multiply_ave_grid_data(const double);
  bigint nextvalid();
  void ave_grid();
};

enum cross_type { GAUSSIAN, ERFC, GAUSSIAN_ERFC, NONE, DELETE };

class Interaction {
 public:
  int i = -1, j = -1;
  int loc = -1;
  cross_type type = NONE;
  std::vector<double> parameters;
};

}    // namespace LAMMPS_NS

#endif
#endif
