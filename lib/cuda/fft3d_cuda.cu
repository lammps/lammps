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

//#define CUDA_PRECISION 1
#include "cuda_precision.h"
#include "cuda_common.h"
struct  FFT_DATA {
  FFT_FLOAT re;
  FFT_FLOAT im;
};

#include "fft3d_cuda_cu.h"
#include "fft3d_cuda_kernel.cu"
#include <stdio.h>

void initfftdata(double* in, FFT_FLOAT* out, int nfast, int nmid, int nslow)
{

  dim3 grid;
  grid.x = nslow;
  grid.y = nmid;
  grid.z = 1;
  dim3 threads;
  threads.x = nfast;
  threads.y = 1;
  threads.z = 1;
  cudaThreadSynchronize();
  initfftdata_kernel <<< grid, threads, 0>>>(in, out);
  cudaThreadSynchronize();
  MYDBG(printf("ERROR-CUDA initfftdata_kernel: %s\n", cudaGetErrorString(cudaGetLastError())));
}


void permute(FFT_DATA* in, FFT_DATA* out, int nfast, int nmid, int nslow)
{

  dim3 grid;
  grid.x = nslow;
  grid.y = nmid;
  grid.z = 1;
  dim3 threads;
  threads.x = nfast * 2;
  threads.y = 1;
  threads.z = 1;
  permute_kernel <<< grid, threads, 0>>>((FFT_FLOAT*)in, (FFT_FLOAT*)out);
  cudaThreadSynchronize();
  MYDBG(printf("ERROR-CUDA permute_kernel: %s\n", cudaGetErrorString(cudaGetLastError())));
}

void permute_scale(FFT_DATA* in, FFT_DATA* out, int nfast, int nmid, int nslow)
{

  dim3 grid;
  grid.x = nslow;
  grid.y = nmid;
  grid.z = 1;
  dim3 threads;
  threads.x = nfast * 2;
  threads.y = 1;
  threads.z = 1;
  permute_kernel <<< grid, threads, 0>>>((FFT_FLOAT*)in, (FFT_FLOAT*)out);
  cudaThreadSynchronize();
}
void permute_part(FFT_DATA* in, FFT_DATA* out, int nfast, int nmid, int nslow, int ihi, int ilo, int jhi, int jlo, int khi, int klo)
{

  dim3 grid;
  grid.x = (ihi - ilo + 1);
  grid.y = (jhi - jlo + 1);
  grid.z = 1;
  dim3 threads;
  threads.x = (khi - klo + 1) * 2;
  threads.y = 1;
  threads.z = 1;
  permute_part_kernel <<< grid, threads, 0>>>((FFT_FLOAT*)in, (FFT_FLOAT*)out, nfast, nmid, nslow, ihi, ilo, jhi, jlo, khi, klo);
  cudaThreadSynchronize();
}

void FFTsyncthreads()
{
  cudaThreadSynchronize();
}

