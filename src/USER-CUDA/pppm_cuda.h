/* ----------------------------------------------------------------------
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

#include "pppm.h"
#include "cuda_data.h"
#include "cuda_precision.h"

namespace LAMMPS_NS {

class PPPMCuda : public PPPM {
 public:
  PPPMCuda(class LAMMPS *, int, char **);
  ~PPPMCuda();
  void init();
  void setup();
  void compute(int, int);
  void timing(int, double &, double &);

  double poissontime;

 protected:
  class Cuda *cuda;
  class FFT3dCuda *fft1c,*fft2c;
  double* work3;
 
  cCudaData<double     , FFT_FLOAT      , x >* cu_work1;
  cCudaData<double     , FFT_FLOAT      , x >* cu_work2;
  cCudaData<double     , FFT_FLOAT      , x >* cu_work3;
  cCudaData<double     , PPPM_FLOAT     , x >* cu_greensfn;
  cCudaData<double     , PPPM_FLOAT     , x >* cu_gf_b;
  cCudaData<double     , PPPM_FLOAT     , x >* cu_fkx;
  cCudaData<double     , PPPM_FLOAT     , x >* cu_fky;
  cCudaData<double     , PPPM_FLOAT     , x >* cu_fkz;
  cCudaData<double     , PPPM_FLOAT     , xy>* cu_vg;
  cCudaData<double     , PPPM_FLOAT     , x >* cu_density_brick;
  cCudaData<int        , int     		, x >* cu_density_brick_int;
  cCudaData<double     , PPPM_FLOAT     , x >* cu_vdx_brick;
  cCudaData<double     , PPPM_FLOAT     , x >* cu_vdy_brick;
  cCudaData<double     , PPPM_FLOAT     , x >* cu_vdz_brick;
  cCudaData<double     , PPPM_FLOAT     , x >* cu_density_fft;
  cCudaData<double     , ENERGY_FLOAT   , x >* cu_energy;
  cCudaData<double     , ENERGY_FLOAT   , x >* cu_virial;
  cCudaData<double     , X_FLOAT   		, yx>* cu_x;
  cCudaData<double     , V_FLOAT   		, yx>* cu_v;
  cCudaData<double     , F_FLOAT   		, yx>* cu_f;	
  cCudaData<double     , F_FLOAT   		, yx>* cu_q;	
  cCudaData<int        , int   			, yx>* cu_part2grid;	
  cCudaData<double	   , PPPM_FLOAT		, x >* cu_rho_coeff;
  cCudaData<PPPM_FLOAT , PPPM_FLOAT		, x >* cu_debugdata;
  cCudaData<int        , int   			, x >* cu_flag;	
  cCudaData<int        , int   			, x >* cu_pppm_grid_n;	
  cCudaData<int        , int   			, x >* cu_pppm_grid_ids;	
  
  ENERGY_FLOAT* slabbuf;
  cCudaData<ENERGY_FLOAT, ENERGY_FLOAT, x >* cu_slabbuf;
  
  int*** density_brick_int;
  PPPM_FLOAT density_intScale;
  int pppm_grid_nmax;
  int* pppm2partgrid;
  int* pppm_grid; 
  PPPM_FLOAT* debugdata;
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
  double qqrd2e;
};

}

#endif
#endif
