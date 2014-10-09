/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator

   Original Version:
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   See the README file in the top-level LAMMPS directory.

   -----------------------------------------------------------------------

   USER-CUDA Package and associated modifications:
   https://sourceforge.net/projects/lammpscuda/

   Christian Trott, christian.trott@tu-ilmenau.de
   Lars Winterfeld, lars.winterfeld@tu-ilmenau.de
   Theoretical Physics II, University of Technology Ilmenau, Germany

   See the README file in the USER-CUDA directory.

   This software is distributed under the GNU General Public License.
------------------------------------------------------------------------- */

#ifdef KSPACE_CLASS

KSpaceStyle(pppm/cuda,PPPMCuda)

#else

#ifndef LMP_PPPM_CUDA_H
#define LMP_PPPM_CUDA_H

#include "pppm_old.h"
#include "cuda_data.h"
#include "cuda_precision.h"

namespace LAMMPS_NS {

class PPPMCuda : public PPPMOld {
 public:
  PPPMCuda(class LAMMPS *, int, char **);
  ~PPPMCuda();
  void init();
  void setup();
  void compute(int, int);
  int timing_1d(int, double &);
  int timing_3d(int, double &);

  double poissontime;

 protected:
  class Cuda *cuda;
  class FFT3dCuda *fft1c,*fft2c;
  double* work3;

  cCudaData<double     , FFT_CFLOAT      , x >* cu_work1;
  cCudaData<double     , FFT_CFLOAT      , x >* cu_work2;
  cCudaData<double     , FFT_CFLOAT      , x >* cu_work3;
  cCudaData<double     , PPPM_CFLOAT     , x >* cu_greensfn;
  cCudaData<double     , PPPM_CFLOAT     , x >* cu_gf_b;
  cCudaData<double     , PPPM_CFLOAT     , x >* cu_fkx;
  cCudaData<double     , PPPM_CFLOAT     , x >* cu_fky;
  cCudaData<double     , PPPM_CFLOAT     , x >* cu_fkz;
  cCudaData<double     , PPPM_CFLOAT     , xy>* cu_vg;
  cCudaData<double     , PPPM_CFLOAT     , x >* cu_density_brick;
  cCudaData<int        , int                     , x >* cu_density_brick_int;
  cCudaData<double     , PPPM_CFLOAT     , x >* cu_vdx_brick;
  cCudaData<double     , PPPM_CFLOAT     , x >* cu_vdy_brick;
  cCudaData<double     , PPPM_CFLOAT     , x >* cu_vdz_brick;
  cCudaData<double     , PPPM_CFLOAT     , x >* cu_density_fft;
  cCudaData<double     , ENERGY_CFLOAT   , x >* cu_energy;
  cCudaData<double     , ENERGY_CFLOAT   , x >* cu_virial;
  cCudaData<double     , X_CFLOAT                   , yx>* cu_x;
  cCudaData<double     , V_CFLOAT                   , yx>* cu_v;
  cCudaData<double     , F_CFLOAT                   , yx>* cu_f;
  cCudaData<double     , F_CFLOAT                   , yx>* cu_q;
  cCudaData<int        , int                           , yx>* cu_part2grid;
  cCudaData<double           , PPPM_CFLOAT                , x >* cu_rho_coeff;
  cCudaData<PPPM_CFLOAT , PPPM_CFLOAT                , x >* cu_debugdata;
  cCudaData<int        , int                           , x >* cu_flag;
  cCudaData<int        , int                           , x >* cu_pppm_grid_n;
  cCudaData<int        , int                           , x >* cu_pppm_grid_ids;

  ENERGY_CFLOAT* slabbuf;
  cCudaData<ENERGY_CFLOAT, ENERGY_CFLOAT, x >* cu_slabbuf;

  int*** density_brick_int;
  PPPM_CFLOAT density_intScale;
  int pppm_grid_nmax;
  int* pppm2partgrid;
  int* pppm_grid;
  PPPM_CFLOAT* debugdata;
  bool firstpass;

  void set_grid();
  void allocate();
  void deallocate();

  virtual void particle_map();
  virtual void make_rho();
  void poisson(int, int);
  virtual void fieldforce();
  virtual void slabcorr(int);
  double*** vdx_brick_tmp;
  int old_nmax;
  int global_flag;
  dev_array* adev_data_array;
};

}

#endif
#endif
