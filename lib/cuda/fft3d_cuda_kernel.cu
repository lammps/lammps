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

__global__ void initfftdata_kernel(double* in, FFT_FLOAT* out)
{
  out[2 * (((blockIdx.x * gridDim.y + blockIdx.y)*blockDim.x) + threadIdx.x)] = in[((blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x) + threadIdx.x];
  out[2 * (((blockIdx.x * gridDim.y + blockIdx.y)*blockDim.x) + threadIdx.x) + 1] = 0;
}


__global__ void permute_kernel(FFT_FLOAT* in, FFT_FLOAT* out)
{
  out[2 * (((threadIdx.x / 2)*gridDim.x + blockIdx.x)*gridDim.y + blockIdx.y) + threadIdx.x - 2 * (threadIdx.x / 2)] = in[((blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x) + threadIdx.x];
}

__global__ void permute_scale_kernel(FFT_FLOAT* in, FFT_FLOAT* out)
{
  out[2 * (((threadIdx.x / 2)*gridDim.x + blockIdx.x)*gridDim.y + blockIdx.y) + threadIdx.x - 2 * (threadIdx.x / 2)] = in[((blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x) + threadIdx.x] * gridDim.x * gridDim.y * blockDim.x * 0.5;
}

__global__ void permute_part_kernel(FFT_FLOAT* in, FFT_FLOAT* out, int nfast, int nmid, int nslow, int ihi, int ilo, int jhi, int jlo, int khi, int klo)
{
  {
    out[2 * ((threadIdx.x / 2) * (ihi - ilo + 1) * (jhi - jlo + 1) + (blockIdx.x) * (jhi - jlo + 1) + blockIdx.y - jlo) + threadIdx.x - 2 * (threadIdx.x / 2)] = in[2 * (blockIdx.x + ilo) * nmid * nslow + 2 * (blockIdx.y + jlo) * nmid + threadIdx.x + 2 * klo];
  }
}
