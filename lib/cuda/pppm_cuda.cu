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

#include "cuda_precision.h"
//#define FFT_CUFFT
#define MY_PREFIX pppm
#include "cuda_shared.h"
#include "cuda_common.h"
#include "pppm_cuda_cu.h"
#include "cuda_runtime.h"
#include <stdio.h>

//#include "crm_cuda_utils.cu"
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

__device__ __constant__ FFT_FLOAT* work1;
__device__ __constant__ FFT_FLOAT* work2;
__device__ __constant__ FFT_FLOAT* work3;
__device__ __constant__ PPPM_FLOAT* greensfn;
__device__ __constant__ PPPM_FLOAT* gf_b;
__device__ __constant__ PPPM_FLOAT* fkx;
__device__ __constant__ PPPM_FLOAT* fky;
__device__ __constant__ PPPM_FLOAT* fkz;
__device__ __constant__ PPPM_FLOAT* vg;
__device__ __constant__ int* part2grid;
__device__ __constant__ PPPM_FLOAT* density_brick;
__device__ __constant__ int* density_brick_int;
__device__ __constant__ PPPM_FLOAT density_intScale;
__device__ __constant__ PPPM_FLOAT* vdx_brick;
__device__ __constant__ PPPM_FLOAT* vdy_brick;
__device__ __constant__ PPPM_FLOAT* vdz_brick;
__device__ __constant__ PPPM_FLOAT* density_fft;
__device__ __constant__ ENERGY_FLOAT* energy;
__device__ __constant__ ENERGY_FLOAT* virial;
__device__ __constant__ int nxlo_in;
__device__ __constant__ int nxhi_in;
__device__ __constant__ int nxlo_out;
__device__ __constant__ int nxhi_out;
__device__ __constant__ int nylo_in;
__device__ __constant__ int nyhi_in;
__device__ __constant__ int nylo_out;
__device__ __constant__ int nyhi_out;
__device__ __constant__ int nzlo_in;
__device__ __constant__ int nzhi_in;
__device__ __constant__ int nzlo_out;
__device__ __constant__ int nzhi_out;
__device__ __constant__ int nxlo_fft;
__device__ __constant__ int nxhi_fft;
__device__ __constant__ int nylo_fft;
__device__ __constant__ int nyhi_fft;
__device__ __constant__ int nzlo_fft;
__device__ __constant__ int nzhi_fft;
__device__ __constant__ int nx_pppm;
__device__ __constant__ int ny_pppm;
__device__ __constant__ int nz_pppm;
__device__ __constant__ int slabflag;
__device__ __constant__ PPPM_FLOAT qqrd2e;
__device__ __constant__ int order;
//__device__ __constant__ float3 sublo;
__device__ __constant__ PPPM_FLOAT* rho_coeff;
__device__ __constant__ int nmax;
__device__ __constant__ int nlocal;
__device__ __constant__ PPPM_FLOAT* debugdata;
__device__ __constant__ PPPM_FLOAT delxinv;
__device__ __constant__ PPPM_FLOAT delyinv;
__device__ __constant__ PPPM_FLOAT delzinv;
__device__ __constant__ int nlower;
__device__ __constant__ int nupper;
__device__ __constant__ PPPM_FLOAT shiftone;


#include "pppm_cuda_kernel.cu"
#include "stdio.h"
void pppm_device_init(void* cu_density_brick, void* cu_vdx_brick, void* cu_vdy_brick, void* cu_vdz_brick, void* cu_density_fft, void* cu_energy, void* cu_virial
                      , void* cu_work1, void* cu_work2, void* cu_work3, void* cu_greensfn, void* cu_fkx, void* cu_fky, void* cu_fkz, void* cu_vg
                      , int cu_nxlo_in, int cu_nxhi_in, int cu_nylo_in, int cu_nyhi_in, int cu_nzlo_in, int cu_nzhi_in, int cu_nxlo_out, int cu_nxhi_out, int cu_nylo_out, int cu_nyhi_out, int cu_nzlo_out, int cu_nzhi_out, int cu_nx_pppm, int cu_ny_pppm, int cu_nz_pppm
                      , int cu_nxlo_fft, int cu_nxhi_fft, int cu_nylo_fft, int cu_nyhi_fft, int cu_nzlo_fft, int cu_nzhi_fft, void* cu_gf_b
                      , double cu_qqrd2e, int cu_order, void* cu_rho_coeff, void* cu_debugdata, void* cu_density_brick_int, int cu_slabflag
                     )
{
  CUT_CHECK_ERROR("ERROR-CUDA poisson_init Start");
  cudaMemcpyToSymbol(density_brick, &cu_density_brick, sizeof(PPPM_FLOAT*));
  cudaMemcpyToSymbol(density_brick_int, &cu_density_brick_int, sizeof(PPPM_FLOAT*));
  cudaMemcpyToSymbol(vdx_brick, &cu_vdx_brick, sizeof(PPPM_FLOAT*));
  cudaMemcpyToSymbol(vdy_brick, &cu_vdy_brick, sizeof(PPPM_FLOAT*));
  cudaMemcpyToSymbol(vdz_brick, &cu_vdz_brick, sizeof(PPPM_FLOAT*));
  cudaMemcpyToSymbol(density_fft, &cu_density_fft, sizeof(PPPM_FLOAT*));
  cudaMemcpyToSymbol(energy, &cu_energy, sizeof(ENERGY_FLOAT*));
  cudaMemcpyToSymbol(virial, &cu_virial, sizeof(ENERGY_FLOAT*));
  cudaMemcpyToSymbol(nxlo_in, &cu_nxlo_in, sizeof(int));
  cudaMemcpyToSymbol(nxhi_in, &cu_nxhi_in, sizeof(int));
  cudaMemcpyToSymbol(nxlo_out, &cu_nxlo_out, sizeof(int));
  cudaMemcpyToSymbol(nxhi_out, &cu_nxhi_out, sizeof(int));
  cudaMemcpyToSymbol(nylo_in, &cu_nylo_in, sizeof(int));
  cudaMemcpyToSymbol(nyhi_in, &cu_nyhi_in, sizeof(int));
  cudaMemcpyToSymbol(nylo_out, &cu_nylo_out, sizeof(int));
  cudaMemcpyToSymbol(nyhi_out, &cu_nyhi_out, sizeof(int));
  cudaMemcpyToSymbol(nzlo_in, &cu_nzlo_in, sizeof(int));
  cudaMemcpyToSymbol(nzhi_in, &cu_nzhi_in, sizeof(int));
  cudaMemcpyToSymbol(nzlo_out, &cu_nzlo_out, sizeof(int));
  cudaMemcpyToSymbol(nzhi_out, &cu_nzhi_out, sizeof(int));
  cudaMemcpyToSymbol(nxlo_fft, &cu_nxlo_fft, sizeof(int));
  cudaMemcpyToSymbol(nxhi_fft, &cu_nxhi_fft, sizeof(int));
  cudaMemcpyToSymbol(nylo_fft, &cu_nylo_fft, sizeof(int));
  cudaMemcpyToSymbol(nyhi_fft, &cu_nyhi_fft, sizeof(int));
  cudaMemcpyToSymbol(nzlo_fft, &cu_nzlo_fft, sizeof(int));
  cudaMemcpyToSymbol(nzhi_fft, &cu_nzhi_fft, sizeof(int));
  cudaMemcpyToSymbol(slabflag, &cu_slabflag, sizeof(int));
  cudaMemcpyToSymbol(nx_pppm, &cu_nx_pppm, sizeof(int));
  cudaMemcpyToSymbol(ny_pppm, &cu_ny_pppm, sizeof(int));
  cudaMemcpyToSymbol(nz_pppm, &cu_nz_pppm, sizeof(int));
  cudaMemcpyToSymbol(work1, &cu_work1, sizeof(FFT_FLOAT*));
  cudaMemcpyToSymbol(work2, &cu_work2, sizeof(FFT_FLOAT*));
  cudaMemcpyToSymbol(work3, &cu_work3, sizeof(FFT_FLOAT*));
  cudaMemcpyToSymbol(greensfn, &cu_greensfn, sizeof(PPPM_FLOAT*));
  cudaMemcpyToSymbol(gf_b, &cu_gf_b, sizeof(PPPM_FLOAT*));
  cudaMemcpyToSymbol(fkx, &cu_fkx, sizeof(PPPM_FLOAT*));
  cudaMemcpyToSymbol(fky, &cu_fky, sizeof(PPPM_FLOAT*));
  cudaMemcpyToSymbol(fkz, &cu_fkz, sizeof(PPPM_FLOAT*));
  cudaMemcpyToSymbol(vg, &cu_vg, sizeof(PPPM_FLOAT*));

  PPPM_FLOAT cu_qqrd2e_a = cu_qqrd2e;
  cudaMemcpyToSymbol(qqrd2e, &cu_qqrd2e_a, sizeof(PPPM_FLOAT));
  cudaMemcpyToSymbol(order, &cu_order, sizeof(int));
  cudaMemcpyToSymbol(rho_coeff, &cu_rho_coeff, sizeof(PPPM_FLOAT*));
  cudaMemcpyToSymbol(debugdata, &cu_debugdata, sizeof(PPPM_FLOAT*));

  CUT_CHECK_ERROR("ERROR-CUDA poisson_init");

  /*if(sizeof(CUDA_FLOAT)==sizeof(float)) printf("PPPMCuda Kernel: Using single precision\n");

  #ifdef PPPM_PRECISION
  if(sizeof(PPPM_FLOAT)==sizeof(float)) printf("PPPMCuda Kernel: Using single precision for pppm core\n");
  if(sizeof(PPPM_FLOAT)==sizeof(double)) printf("PPPMCuda Kernel: Using double precision for pppm core\n");
  #endif
  #ifdef ENERGY_PRECISION
  if(sizeof(ENERGY_FLOAT)==sizeof(float)) printf("PPPMCuda Kernel: Using single precision for energy\n");
  if(sizeof(ENERGY_FLOAT)==sizeof(double)) printf("PPPMCuda Kernel: Using double precision for energy\n");
  #endif
  #ifdef ENERGY_PRECISION
  if(sizeof(FFT_FLOAT)==sizeof(float)) printf("PPPMCuda Kernel: Using single precision for fft\n");
  if(sizeof(FFT_FLOAT)==sizeof(double)) printf("PPPMCuda Kernel: Using double precision for fft\n");
  #endif
  #ifdef X_PRECISION
  if(sizeof(X_FLOAT)==sizeof(float)) printf("PPPMCuda Kernel: Using single precision for positions\n");
  if(sizeof(X_FLOAT)==sizeof(double)) printf("PPPMCuda Kernel: Using double precision for positions\n");
  #endif
  #ifdef F_PRECISION
  if(sizeof(F_FLOAT)==sizeof(float)) printf("PPPMCuda Kernel: Using single precision for forces\n");
  if(sizeof(F_FLOAT)==sizeof(double)) printf("PPPMCuda Kernel: Using double precision for forces\n");
  #endif*/
}

void pppm_device_init_setup(cuda_shared_data* sdata, PPPM_FLOAT cu_shiftone, PPPM_FLOAT cu_delxinv, PPPM_FLOAT cu_delyinv, PPPM_FLOAT cu_delzinv, int cu_nlower, int cu_nupper)
{
  cudaMemcpyToSymbol(delxinv, &cu_delxinv, sizeof(PPPM_FLOAT));
  cudaMemcpyToSymbol(delyinv, &cu_delyinv, sizeof(PPPM_FLOAT));
  cudaMemcpyToSymbol(delzinv, &cu_delzinv, sizeof(PPPM_FLOAT));
  cudaMemcpyToSymbol(shiftone, &cu_shiftone, sizeof(PPPM_FLOAT));
  cudaMemcpyToSymbol(nlower, &cu_nlower, sizeof(int));
  cudaMemcpyToSymbol(nupper, &cu_nupper, sizeof(int));
  cudaMemcpyToSymbol(MY_AP(sublo)   , sdata->domain.sublo, 3 * sizeof(X_FLOAT));
  cudaMemcpyToSymbol(MY_AP(subhi)   , sdata->domain.subhi, 3 * sizeof(X_FLOAT));
  cudaMemcpyToSymbol(MY_AP(boxlo)   , sdata->domain.boxlo, 3 * sizeof(X_FLOAT));
  CUT_CHECK_ERROR("ERROR-CUDA pppm_init_setup");
}

void pppm_device_update(cuda_shared_data* sdata, void* cu_part2grid, int nlocala, int nmaxa)
{
  cudaMemcpyToSymbol(part2grid, &cu_part2grid, sizeof(int*));
  cudaMemcpyToSymbol(MY_AP(x)   , & sdata->atom.x   .dev_data, sizeof(X_FLOAT*));
  cudaMemcpyToSymbol(MY_AP(f)   , & sdata->atom.f   .dev_data, sizeof(F_FLOAT*));
  cudaMemcpyToSymbol(MY_AP(q)   , & sdata->atom.q   .dev_data, sizeof(F_FLOAT*));
  cudaMemcpyToSymbol(MY_AP(tag)   , & sdata->atom.tag   .dev_data, sizeof(int*));
  //cudaMemcpyToSymbol(MY_AP(nlocal)   , & sdata->atom.nlocal   .dev_data, sizeof(int));
  cudaMemcpyToSymbol(nlocal   , &nlocala, sizeof(int));
  cudaMemcpyToSymbol(nmax   , &nmaxa, sizeof(int));
  CUT_CHECK_ERROR("ERROR-CUDA pppm_device_update");

}

void pppm_update_nlocal(int nlocala)
{
  cudaMemcpyToSymbol(nlocal   , &nlocala, sizeof(int));
  CUT_CHECK_ERROR("ERROR-CUDA update_nlocal b");
}


void Cuda_PPPM_Setup_fkxyz_vg(int nx_pppma, int ny_pppma, int nz_pppma, PPPM_FLOAT unitkx, PPPM_FLOAT unitky, PPPM_FLOAT unitkz, PPPM_FLOAT g_ewald)
{
  dim3 grid;
  dim3 threads;
  grid.x = nz_pppma;
  grid.y = ny_pppma;
  grid.z = 1;
  threads.x = nx_pppma;
  threads.y = 1;
  threads.z = 1;
  setup_fkxyz_vg <<< grid, threads, 0>>>(unitkx, unitky, unitkz, g_ewald);
  cudaThreadSynchronize();

  CUT_CHECK_ERROR("ERROR-CUDA Cuda_PPPM_Setup_fkxyz_vg ");
}

void Cuda_PPPM_setup_greensfn(int nx_pppma, int ny_pppma, int nz_pppma, PPPM_FLOAT unitkx, PPPM_FLOAT unitky, PPPM_FLOAT unitkz, PPPM_FLOAT g_ewald,
                              int nbx, int nby, int nbz, PPPM_FLOAT xprd, PPPM_FLOAT yprd, PPPM_FLOAT zprd_slab)
{
  dim3 grid;
  dim3 threads;
  grid.x = nz_pppma;
  grid.y = ny_pppma;
  grid.z = 1;
  threads.x = nx_pppma;
  threads.y = 1;
  threads.z = 1;
  setup_greensfn <<< grid, threads, 0>>>(unitkx, unitky, unitkz, g_ewald, nbx, nby, nbz, xprd, yprd, zprd_slab);
  cudaThreadSynchronize();
  CUT_CHECK_ERROR("ERROR-CUDA Cuda_PPPM_Setup_greensfn ");
}

void poisson_scale(int nx_pppma, int ny_pppma, int nz_pppma)
{
  dim3 grid;
  dim3 threads;
  grid.x = nz_pppma;
  grid.y = ny_pppma;
  grid.z = 1;
  threads.x = nx_pppma;
  threads.y = 1;
  threads.z = 1;
  poisson_scale_kernel <<< grid, threads, 0>>>();
  CUT_CHECK_ERROR("ERROR-CUDA poisson_scale ");

}

void poisson_xgrad(int nx_pppma, int ny_pppma, int nz_pppma)
{
  dim3 grid;
  dim3 threads;
  grid.x = nz_pppma;
  grid.y = ny_pppma;
  grid.z = 1;
  threads.x = nx_pppma;
  threads.y = 1;
  threads.z = 1;
  poisson_xgrad_kernel <<< grid, threads, 0>>>();
  CUT_CHECK_ERROR("ERROR-CUDA poisson_xgrad ");
}

void poisson_ygrad(int nx_pppma, int ny_pppma, int nz_pppma)
{
  dim3 grid;
  dim3 threads;
  grid.x = nz_pppma;
  grid.y = ny_pppma;
  grid.z = 1;
  threads.x = nx_pppma;
  threads.y = 1;
  threads.z = 1;
  poisson_ygrad_kernel <<< grid, threads, 0>>>();
  CUT_CHECK_ERROR("ERROR-CUDA poisson_ygrad ");
}

void poisson_zgrad(int nx_pppma, int ny_pppma, int nz_pppma)
{
  dim3 grid;
  dim3 threads;
  grid.x = nz_pppma;
  grid.y = ny_pppma;
  grid.z = 1;
  threads.x = nx_pppma;
  threads.y = 1;
  threads.z = 1;
  poisson_zgrad_kernel <<< grid, threads, 0>>>();
  CUT_CHECK_ERROR("ERROR-CUDA poisson_zgrad ");
}

void poisson_vdx_brick(int ihi, int ilo, int jhi, int jlo, int khi, int klo, int nx_pppma, int ny_pppma, int nz_pppma)
{

  dim3 grid;
  dim3 threads;
  grid.x = khi - klo + 1;
  grid.y = jhi - jlo + 1;
  grid.z = 1;
  threads.x = ihi - ilo + 1;
  threads.y = 1;
  threads.z = 1;
  //printf("VDX_BRICK CUDA: %i %i %i\n",grid.x,grid.y,threads.x);
  poisson_vdx_brick_kernel <<< grid, threads, 0>>>(ilo, jlo, klo);
  CUT_CHECK_ERROR("ERROR-CUDA poisson_vdxbrick ");
  cudaThreadSynchronize();
}

void poisson_vdy_brick(int ihi, int ilo, int jhi, int jlo, int khi, int klo, int nx_pppm, int ny_pppm, int nz_pppm)
{
  dim3 grid;
  dim3 threads;
  grid.x = khi - klo + 1;
  grid.y = jhi - jlo + 1;
  grid.z = 1;
  threads.x = ihi - ilo + 1;
  threads.y = 1;
  threads.z = 1;
  poisson_vdy_brick_kernel <<< grid, threads, 0>>>(ilo, jlo, klo);
  CUT_CHECK_ERROR("ERROR-CUDA poisson_vdybrick ");
  cudaThreadSynchronize();
}

void poisson_vdz_brick(int ihi, int ilo, int jhi, int jlo, int khi, int klo, int nx_pppm, int ny_pppm, int nz_pppm)
{
  dim3 grid;
  dim3 threads;
  grid.x = khi - klo + 1;
  grid.y = jhi - jlo + 1;
  grid.z = 1;
  threads.x = ihi - ilo + 1;
  threads.y = 1;
  threads.z = 1;
  poisson_vdz_brick_kernel <<< grid, threads, 0>>>(ilo, jlo, klo);
  CUT_CHECK_ERROR("ERROR-CUDA poisson_vdzbrick ");
  cudaThreadSynchronize();
}


void poisson_energy(int nxlo_fft, int nxhi_fft, int nylo_fft, int nyhi_fft, int nzlo_fft, int nzhi_fft, int vflag)
{
  //printf("VFLAG_GPU: %i\n",vflag);
  CUT_CHECK_ERROR("ERROR-CUDA poisson_energy start ");
  dim3 grid;
  dim3 threads;
  grid.x = nzhi_fft - nzlo_fft + 1;
  grid.y = nyhi_fft - nylo_fft + 1;
  grid.z = 1;
  threads.x = nxhi_fft - nxlo_fft + 1;
  threads.y = 1;
  threads.z = 1;
  poisson_energy_kernel <<< grid, threads, threads.x* sizeof(ENERGY_FLOAT)>>>(nxlo_fft, nylo_fft, nzlo_fft, vflag);

  cudaThreadSynchronize();
  CUT_CHECK_ERROR("ERROR-CUDA poisson_energy end ");
}

ENERGY_FLOAT sum_energy(void* cu_virial, void* cu_energy, int nx_pppma, int ny_pppma, int nz_pppma, int vflag, ENERGY_FLOAT* cpu_virial)
{
  ENERGY_FLOAT host_energy = 0;
  dim3 grid;
  dim3 threads;

  grid.x = nz_pppma;
  grid.y = 1;
  grid.z = 1;
  threads.x = ny_pppma;
  threads.y = 1;
  threads.z = 1;
  sum_energy_kernel1 <<< grid, threads, ny_pppma* sizeof(ENERGY_FLOAT)>>>(vflag);
  cudaThreadSynchronize();
  CUT_CHECK_ERROR("ERROR-CUDA sumenergy_kernel1 ");

  grid.x = 1;
  grid.y = 1;
  grid.z = 1;
  threads.x = nz_pppma;
  threads.y = 1;
  threads.z = 1;
  sum_energy_kernel2 <<< grid, threads, nz_pppma* sizeof(ENERGY_FLOAT)>>>(vflag);
  cudaThreadSynchronize();
  CUT_CHECK_ERROR("ERROR-CUDA sumenergy_kernel2 ");

  cudaMemcpy((void*)(&host_energy), cu_energy, sizeof(ENERGY_FLOAT), cudaMemcpyDeviceToHost);

  if(vflag)
    cudaMemcpy((void*) cpu_virial, (void*) cu_virial, 6 * sizeof(ENERGY_FLOAT), cudaMemcpyDeviceToHost);
  CUT_CHECK_ERROR("ERROR-CUDA sumenergy_memcopy");

  return host_energy;
}

void cuda_make_rho(cuda_shared_data* sdata, void* flag, PPPM_FLOAT* cu_density_intScale, int ihi, int ilo, int jhi, int jlo, int khi, int klo, void* cu_density_brick, void* cu_density_brick_int)
{
  CUT_CHECK_ERROR("cuda_make_rho begin");
  dim3 grid, threads;
  int cpu_flag[3];
  grid.x = (sdata->atom.nlocal + 31) / 32;
  grid.y = 1;
  grid.z = 1;
  threads.x = 32;
  threads.y = 1;
  threads.z = 1;
  int sharedmemsize = (32 + 32 * (sdata->pppm.nupper - sdata->pppm.nlower + 1) + sdata->pppm.order * (sdata->pppm.order / 2 - (1 - sdata->pppm.order) / 2 + 1)) * sizeof(PPPM_FLOAT);

  do {
    cpu_flag[0] = 0;
    cpu_flag[1] = 0;
    cpu_flag[2] = 0;
    cudaMemcpyToSymbol(density_intScale, cu_density_intScale, sizeof(PPPM_FLOAT*));
    CUT_CHECK_ERROR("ERROR-CUDA make_rho pre Z");
    cudaMemset(flag, 0, 3 * sizeof(int));
    CUT_CHECK_ERROR("ERROR-CUDA make_rho pre A");
    cudaMemset(cu_density_brick, 0, (khi - klo + 1) * (jhi - jlo + 1) * (ihi - ilo + 1)*sizeof(PPPM_FLOAT));
    CUT_CHECK_ERROR("ERROR-CUDA make_rho pre B");
    cudaMemset(cu_density_brick_int, 0, (khi - klo + 1) * (jhi - jlo + 1) * (ihi - ilo + 1)*sizeof(int));
    CUT_CHECK_ERROR("ERROR-CUDA make_rho pre C");
    make_rho_kernel <<< grid, threads, sharedmemsize>>>((int*) flag, 32 / (sdata->pppm.nupper - sdata->pppm.nlower + 1));
    cudaThreadSynchronize();
    CUT_CHECK_ERROR("ERROR-CUDA make_rho A");
    cudaMemcpy((void*) &cpu_flag, flag, 3 * sizeof(int), cudaMemcpyDeviceToHost);

    if(cpu_flag[0] != 0) {
      (*cu_density_intScale) /= 2;
      MYDBG(printf("PPPM_Cuda::cuda_make_rho: Decrease cu_density_intScale to: %e\n", *cu_density_intScale);)
    }
    if((cpu_flag[0] == 0) && (cpu_flag[1] == 0)) {
      (*cu_density_intScale) *= 2;
      MYDBG(printf("PPPM_Cuda::cuda_make_rho: Increase cu_density_intScale to: %e\n", *cu_density_intScale);)
    }
    /* if((*cu_density_intScale)>0xe0000000)
     {
     	printf("Error Scaling\n");
         cpu_flag[0]=0;
         cpu_flag[1]=1;
     }*/
    CUT_CHECK_ERROR("ERROR-CUDA make_rho B");
  } while((cpu_flag[0] != 0) || (cpu_flag[1] == 0));


  grid.x = khi - klo + 1;
  grid.y = jhi - jlo + 1;
  threads.x = ihi - ilo + 1;
  scale_rho_kernel <<< grid, threads, 0>>>();
  cudaThreadSynchronize();
  CUT_CHECK_ERROR("ERROR-CUDA make_rho_scale");
}


int cuda_particle_map(cuda_shared_data* sdata, void* flag)
{
  dim3 grid, threads;
  int cpu_flag;
  grid.x = (sdata->atom.nlocal + 31) / 32;
  grid.y = 1;
  grid.z = 1;
  threads.x = 32;
  threads.y = 1;
  threads.z = 1;
  CUT_CHECK_ERROR("ERROR-CUDA particla_map ..pre");
  particle_map_kernel <<< grid, threads, 0>>>((int*) flag);
  cudaThreadSynchronize();
  CUT_CHECK_ERROR("ERROR-CUDA particla_map a");
  cudaMemcpy((void*) &cpu_flag, flag, sizeof(int), cudaMemcpyDeviceToHost);
  CUT_CHECK_ERROR("ERROR-CUDA particla_map b");
  return cpu_flag;
}


void cuda_fieldforce(cuda_shared_data* sdata, void* flag)
{
  dim3 grid, threads;
  grid.x = (sdata->atom.nlocal + 31) / 32;
  grid.y = 1;
  grid.z = 1;
  threads.x = 32;
  threads.y = 1;
  threads.z = 1;
  int sharedmemsize = (32 + 3 * 32 * (sdata->pppm.nupper - sdata->pppm.nlower + 1) + sdata->pppm.order * (sdata->pppm.order / 2 - (1 - sdata->pppm.order) / 2 + 1)) * sizeof(PPPM_FLOAT);
  fieldforce_kernel <<< grid, threads, sharedmemsize>>>
  (sdata->pppm.nupper - sdata->pppm.nlower + 1, 32 / (sdata->pppm.nupper - sdata->pppm.nlower + 1), (int*) flag);
  cudaThreadSynchronize();
  CUT_CHECK_ERROR("ERROR-CUDA fieldforce");
}

double cuda_slabcorr_energy(cuda_shared_data* sdata, ENERGY_FLOAT* buf, ENERGY_FLOAT* dev_buf)
{
  dim3 grid, threads;
  grid.x = (sdata->atom.nlocal + 31) / 32;
  grid.y = 1;
  grid.z = 1;
  threads.x = 32;
  threads.y = 1;
  threads.z = 1;
  slabcorr_energy_kernel <<< grid, threads, 32* sizeof(ENERGY_FLOAT)>>>(dev_buf);
  cudaThreadSynchronize();
  cudaMemcpy((void*) buf, dev_buf, grid.x* sizeof(ENERGY_FLOAT), cudaMemcpyDeviceToHost);

  double dipole_all = 0.0;

  for(int i = 0; i < grid.x; i++)
    dipole_all += buf[i];

  return dipole_all;
}

void cuda_slabcorr_force(cuda_shared_data* sdata, F_FLOAT ffact)
{
  dim3 grid, threads;
  grid.x = (sdata->atom.nlocal + 31) / 32;
  grid.y = 1;
  grid.z = 1;
  threads.x = 32;
  threads.y = 1;
  threads.z = 1;
  slabcorr_force_kernel <<< grid, threads>>>(ffact);
  cudaThreadSynchronize();
}

void sum_virial(double* host_virial)
{
}

void pppm_initfftdata(cuda_shared_data* sdata, PPPM_FLOAT* in, FFT_FLOAT* out)
{
  int nslow = sdata->pppm.nzhi_in - sdata->pppm.nzlo_in;
  int nmid = sdata->pppm.nyhi_in - sdata->pppm.nylo_in;
  int nfast = sdata->pppm.nxhi_in - sdata->pppm.nxlo_in;
  int nrimz = MAX(sdata->pppm.nzlo_in - sdata->pppm.nzlo_out, sdata->pppm.nzhi_out - sdata->pppm.nzhi_in);
  int nrimy = MAX(sdata->pppm.nylo_in - sdata->pppm.nylo_out, sdata->pppm.nyhi_out - sdata->pppm.nyhi_in);
  int nrimx = MAX(sdata->pppm.nxlo_in - sdata->pppm.nxlo_out, sdata->pppm.nxhi_out - sdata->pppm.nxhi_in);
  dim3 grid;
  grid.x = nslow + 1;
  grid.y = nmid + 1;
  grid.z = 1;
  dim3 threads;
  threads.x = nfast + 1;
  threads.y = 1;
  threads.z = 1;
  cudaThreadSynchronize();
  initfftdata_core_kernel <<< grid, threads, 0>>>(in, out);
  cudaThreadSynchronize();
  grid.x = nrimz;
  grid.y = nmid + 1;
  threads.x = nfast + 1;
  initfftdata_z_kernel <<< grid, threads, 0>>>(in, out);
  cudaThreadSynchronize();
  grid.x = nslow + 1;
  grid.y = nrimy;
  threads.x = nfast + 1;
  initfftdata_y_kernel <<< grid, threads, 0>>>(in, out);
  cudaThreadSynchronize();
  grid.x = nslow + 1;
  grid.y = nmid + 1;
  threads.x = nrimx;
  initfftdata_x_kernel <<< grid, threads, 0>>>(in, out);
  cudaThreadSynchronize();
  grid.x = nrimz;
  grid.y = nrimy;
  threads.x = nfast + 1;
  initfftdata_yz_kernel <<< grid, threads, 0>>>(in, out);
  cudaThreadSynchronize();
  grid.x = nrimz;
  grid.y = nmid + 1;
  threads.x = nrimx;
  initfftdata_xz_kernel <<< grid, threads, 0>>>(in, out);
  cudaThreadSynchronize();
  grid.x = nslow + 1;
  grid.y = nrimy;
  threads.x = nrimx;
  initfftdata_xy_kernel <<< grid, threads, 0>>>(in, out);
  cudaThreadSynchronize();
  grid.x = nrimz;
  grid.y = nrimy;
  threads.x = nrimx;
  initfftdata_xyz_kernel <<< grid, threads, 0>>>(in, out);
  cudaThreadSynchronize();
  CUT_CHECK_ERROR("ERROR-CUDA initfftdata_kernel");
}


