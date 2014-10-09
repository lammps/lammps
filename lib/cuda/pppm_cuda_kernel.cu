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

#define OFFSET 4096
__device__ int negativCUDA(float f)
{
  return ((unsigned int)1 << 31 & (__float_as_int(f))) >> 31;
}

__device__ void reduceBlock(float* data)
{
  int p2 = 1;

  while(p2 * 2 < blockDim.x) p2 *= 2;

  if(threadIdx.x < blockDim.x - p2)
    data[threadIdx.x] += data[threadIdx.x + p2];

  __syncthreads();

  for(int i = 2; i <= p2; i *= 2) {
    if(threadIdx.x < p2 / i)
      data[threadIdx.x] += data[threadIdx.x + p2 / i];

    __syncthreads();
  }
}

__device__ void reduceBlock(double* data)
{
  int p2 = 1;

  while(p2 * 2 < blockDim.x) p2 *= 2;

  if(threadIdx.x < blockDim.x - p2)
    data[threadIdx.x] += data[threadIdx.x + p2];

  __syncthreads();

  for(int i = 2; i <= p2; i *= 2) {
    if(threadIdx.x < p2 / i)
      data[threadIdx.x] += data[threadIdx.x + p2 / i];

    __syncthreads();
  }
}

extern __shared__ PPPM_CFLOAT sharedmem[];

__global__ void setup_fkxyz_vg(PPPM_CFLOAT unitkx, PPPM_CFLOAT unitky, PPPM_CFLOAT unitkz, PPPM_CFLOAT g_ewald)
{
  PPPM_CFLOAT my_fkx = unitkx * (int(threadIdx.x) - nx_pppm * (2 * int(threadIdx.x) / nx_pppm));
  PPPM_CFLOAT my_fky = unitky * (int(blockIdx.y) - ny_pppm * (2 * int(blockIdx.y) / ny_pppm));
  PPPM_CFLOAT my_fkz = unitkz * (int(blockIdx.x) - nz_pppm * (2 * int(blockIdx.x) / nz_pppm));

  if((blockIdx.x == 0) && (blockIdx.y == 0)) fkx[threadIdx.x] = my_fkx;

  if((blockIdx.x == 0) && (threadIdx.x == 0)) fky[blockIdx.y] = my_fky;

  if((threadIdx.x == 0) && (blockIdx.y == 0)) fkz[blockIdx.x] = my_fkz;

  __syncthreads();

  if((blockIdx.x >= nzlo_fft) && (blockIdx.x <= nzhi_fft) &&
      (blockIdx.y >= nylo_fft) && (blockIdx.y <= nyhi_fft) &&
      (threadIdx.x >= nxlo_fft) && (threadIdx.x <= nxhi_fft)) {
    int n = ((int(blockIdx.x) - nzlo_fft) * (nyhi_fft - nylo_fft + 1) + int(blockIdx.y) - nylo_fft) * (nxhi_fft - nxlo_fft + 1) + int(threadIdx.x) - nxlo_fft;
    PPPM_CFLOAT sqk = my_fkx * my_fkx + my_fky * my_fky + my_fkz * my_fkz;
    PPPM_CFLOAT vterm = (sqk == PPPM_F(0.0)) ? PPPM_F(0.0) : PPPM_F(-2.0) * (PPPM_F(1.0) / sqk + PPPM_F(0.25) / (g_ewald * g_ewald));
    vg[6 * n + 0] = (sqk == PPPM_F(0.0)) ? PPPM_F(0.0) : PPPM_F(1.0) + vterm * my_fkx * my_fkx;
    vg[6 * n + 1] = (sqk == PPPM_F(0.0)) ? PPPM_F(0.0) : PPPM_F(1.0) + vterm * my_fky * my_fky;
    vg[6 * n + 2] = (sqk == PPPM_F(0.0)) ? PPPM_F(0.0) : PPPM_F(1.0) + vterm * my_fkz * my_fkz;
    vg[6 * n + 3] = (sqk == PPPM_F(0.0)) ? PPPM_F(0.0) : vterm * my_fkx * my_fky;
    vg[6 * n + 4] = (sqk == PPPM_F(0.0)) ? PPPM_F(0.0) : vterm * my_fkx * my_fkz;
    vg[6 * n + 5] = (sqk == PPPM_F(0.0)) ? PPPM_F(0.0) : vterm * my_fky * my_fkz;

  }
}

__device__ PPPM_CFLOAT gf_denom(PPPM_CFLOAT x, PPPM_CFLOAT y, PPPM_CFLOAT z)
{
  PPPM_CFLOAT sx, sy, sz;
  sz = sy = sx = PPPM_F(0.0);

  for(int l = order - 1; l >= 0; l--) {
    sx = gf_b[l] + sx * x;
    sy = gf_b[l] + sy * y;
    sz = gf_b[l] + sz * z;
  }

  PPPM_CFLOAT s = sx * sy * sz;
  return s * s;
}

__global__ void setup_greensfn(PPPM_CFLOAT unitkx, PPPM_CFLOAT unitky, PPPM_CFLOAT unitkz, PPPM_CFLOAT g_ewald,
                               int nbx, int nby, int nbz,
                               PPPM_CFLOAT xprd, PPPM_CFLOAT yprd, PPPM_CFLOAT zprd_slab)
{
  PPPM_CFLOAT sqk;
  int nx, ny, nz, kper, lper, mper, k, l, m;
  PPPM_CFLOAT snx, sny, snz, snx2, sny2, snz2;
  PPPM_CFLOAT argx, argy, argz, wx, wy, wz, sx, sy, sz, qx, qy, qz;
  PPPM_CFLOAT sum1, dot1, dot2;
  PPPM_CFLOAT numerator, denominator;

  PPPM_CFLOAT form = PPPM_F(1.0);
  int n = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;
  m = blockIdx.x;
  l = blockIdx.y;
  k = threadIdx.x;

  mper = m - nz_pppm * (2 * m / nz_pppm);
  snz = sin(PPPM_F(0.5) * unitkz * mper * zprd_slab / nz_pppm);
  snz2 = snz * snz;


  lper = l - ny_pppm * (2 * l / ny_pppm);
  sny = sin(PPPM_F(0.5) * unitky * lper * yprd / ny_pppm);
  sny2 = sny * sny;

  kper = k - nx_pppm * (2 * k / nx_pppm);
  snx = sin(PPPM_F(0.5) * unitkx * kper * xprd / nx_pppm);
  snx2 = snx * snx;

  sqk = pow(unitkx * kper, PPPM_F(2.0)) + pow(unitky * lper, PPPM_F(2.0)) +
        pow(unitkz * mper, PPPM_F(2.0));

  if(sqk != PPPM_F(0.0)) {
    numerator = form * PPPM_F(12.5663706) / sqk;
    denominator = gf_denom(snx2, sny2, snz2);
    sum1 = PPPM_F(0.0);

    for(nx = -nbx; nx <= nbx; nx++) {
      qx = unitkx * (kper + nx_pppm * nx);
      sx = exp(PPPM_F(-.25) * pow(qx / g_ewald, PPPM_F(2.0)));
      wx = PPPM_F(1.0);
      argx = PPPM_F(0.5) * qx * xprd / nx_pppm;

      if(argx != PPPM_F(0.0)) wx = pow(sin(argx) / argx, order);

      for(ny = -nby; ny <= nby; ny++) {
        qy = unitky * (lper + ny_pppm * ny);
        sy = exp(PPPM_F(-.25) * pow(qy / g_ewald, PPPM_F(2.0)));
        wy = PPPM_F(1.0);
        argy = PPPM_F(0.5) * qy * yprd / ny_pppm;

        if(argy != PPPM_F(0.0)) wy = pow(sin(argy) / argy, order);

        for(nz = -nbz; nz <= nbz; nz++) {
          qz = unitkz * (mper + nz_pppm * nz);
          sz = exp(PPPM_F(-.25) * pow(qz / g_ewald, PPPM_F(2.0)));
          wz = PPPM_F(1.0);
          argz = PPPM_F(0.5) * qz * zprd_slab / nz_pppm;

          if(argz != PPPM_F(0.0)) wz = pow(sin(argz) / argz, order);

          dot1 = unitkx * kper * qx + unitky * lper * qy + unitkz * mper * qz;
          dot2 = qx * qx + qy * qy + qz * qz;
          sum1 += (dot1 / dot2) * sx * sy * sz * pow(wx * wy * wz, PPPM_F(2.0));
        }
      }
    }

    greensfn[n] = numerator * sum1 / denominator;
  } else greensfn[n] = PPPM_F(0.0);
}

__global__ void poisson_scale_kernel()
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;
  FFT_CFLOAT scaleinv = FFT_F(1.0) / (gridDim.x * gridDim.y * blockDim.x);
  work1[2 * i] *= scaleinv * greensfn[i];
  work1[2 * i + 1] *= scaleinv * greensfn[i];
}

__global__ void poisson_xgrad_kernel()
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;
  work2[2 * i] = fkx[threadIdx.x] * work1[2 * i + 1];
  work2[2 * i + 1] = -fkx[threadIdx.x] * work1[2 * i];
}

__global__ void poisson_ygrad_kernel()
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;
  work2[2 * i] = fky[blockIdx.y] * work1[2 * i + 1];
  work2[2 * i + 1] = -fky[blockIdx.y] * work1[2 * i];
}

__global__ void poisson_zgrad_kernel()
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;
  work2[2 * i] = fkz[blockIdx.x] * work1[2 * i + 1];
  work2[2 * i + 1] = -fkz[blockIdx.x] * work1[2 * i];
}

__global__ void poisson_vdx_brick_kernel(int ilo, int jlo, int klo)
{
  int k = blockIdx.x + klo;
  k += nz_pppm * negativCUDA(CUDA_F(1.0) * k) - nz_pppm * negativCUDA(CUDA_F(1.0) * (nz_pppm - k - 1));
  int j = blockIdx.y + jlo;
  j += ny_pppm * negativCUDA(CUDA_F(1.0) * j) - ny_pppm * negativCUDA(CUDA_F(1.0) * (ny_pppm - j - 1));
  int i = threadIdx.x + ilo;
  i += nx_pppm * negativCUDA(CUDA_F(1.0) * i) - nx_pppm * negativCUDA(CUDA_F(1.0) * (nx_pppm - i - 1));
  vdx_brick[((blockIdx.x) * (nyhi_out - nylo_out + 1) + blockIdx.y) * (nxhi_out - nxlo_out + 1) + threadIdx.x] = work3[2 * (((k) * ny_pppm + (j)) * nx_pppm + i)];
}

__global__ void poisson_vdy_brick_kernel(int ilo, int jlo, int klo)
{
  int k = blockIdx.x + klo;
  k += nz_pppm * negativCUDA(CUDA_F(1.0) * k) - nz_pppm * negativCUDA(CUDA_F(1.0) * (nz_pppm - k - 1));
  int j = blockIdx.y + jlo;
  j += ny_pppm * negativCUDA(CUDA_F(1.0) * j) - ny_pppm * negativCUDA(CUDA_F(1.0) * (ny_pppm - j - 1));
  int i = threadIdx.x + ilo;
  i += nx_pppm * negativCUDA(CUDA_F(1.0) * i) - nx_pppm * negativCUDA(CUDA_F(1.0) * (nx_pppm - i - 1));
  vdy_brick[((blockIdx.x) * (nyhi_out - nylo_out + 1) + blockIdx.y) * (nxhi_out - nxlo_out + 1) + threadIdx.x] = work3[2 * (((k) * ny_pppm + (j)) * nx_pppm + i)];
}

__global__ void poisson_vdz_brick_kernel(int ilo, int jlo, int klo)
{
  int k = blockIdx.x + klo;
  k += nz_pppm * negativCUDA(CUDA_F(1.0) * k) - nz_pppm * negativCUDA(CUDA_F(1.0) * (nz_pppm - k - 1));
  int j = blockIdx.y + jlo;
  j += ny_pppm * negativCUDA(CUDA_F(1.0) * j) - ny_pppm * negativCUDA(CUDA_F(1.0) * (ny_pppm - j - 1));
  int i = threadIdx.x + ilo;
  i += nx_pppm * negativCUDA(CUDA_F(1.0) * i) - nx_pppm * negativCUDA(CUDA_F(1.0) * (nx_pppm - i - 1));
  vdz_brick[((blockIdx.x) * (nyhi_out - nylo_out + 1) + blockIdx.y) * (nxhi_out - nxlo_out + 1) + threadIdx.x] = work3[2 * (((k) * ny_pppm + (j)) * nx_pppm + i)];
}

__global__ void poisson_energy_kernel(int nxlo_fft, int nylo_fft, int nzlo_fft, int vflag)
{
  ENERGY_CFLOAT scaleinv = FFT_F(1.0) / (nx_pppm * ny_pppm * nz_pppm);
  int i = (blockIdx.x + nzlo_fft) * ny_pppm * nx_pppm + (blockIdx.y + nylo_fft) * nx_pppm + threadIdx.x + nxlo_fft;
  ENERGY_CFLOAT* s_energy = (ENERGY_CFLOAT*) sharedmem;
  ENERGY_CFLOAT myenergy = scaleinv * scaleinv * greensfn[i] * (work1[2 * i] * work1[2 * i] + work1[2 * i + 1] * work1[2 * i + 1]);
  s_energy[threadIdx.x] = myenergy;

  __syncthreads();
  reduceBlock(s_energy);

  if(threadIdx.x == 0)
    energy[blockIdx.x * ny_pppm + blockIdx.y] = s_energy[0];

  if(vflag) {
    __syncthreads();

    for(int j = 0; j < 6; j++) {
      s_energy[threadIdx.x] = myenergy * vg[((blockIdx.x * gridDim.y + blockIdx.y) * (blockDim.x) + threadIdx.x) * 6 + j];
      __syncthreads();
      reduceBlock(s_energy);

      if(threadIdx.x == 0)
        virial[blockIdx.x * ny_pppm + blockIdx.y + j * nz_pppm * ny_pppm] = s_energy[0];
    }
  }
}


__global__ void sum_energy_kernel1(int vflag)
{
  ENERGY_CFLOAT myenergy = energy[(blockIdx.x * ny_pppm + threadIdx.x)];
  ENERGY_CFLOAT* s_energy = (ENERGY_CFLOAT*) sharedmem;
  s_energy[threadIdx.x] = myenergy;
  __syncthreads();
  reduceBlock(s_energy);

  if(threadIdx.x == 0)
    energy[blockIdx.x * ny_pppm] = s_energy[0];

  if(vflag) {
    __syncthreads();

    for(int j = 0; j < 6; j++) {
      myenergy = virial[blockIdx.x * ny_pppm + threadIdx.x + j * ny_pppm * nz_pppm];
      s_energy[threadIdx.x] = myenergy;
      __syncthreads();
      reduceBlock(s_energy);

      if(threadIdx.x == 0)
        virial[blockIdx.x * ny_pppm + j * ny_pppm * nz_pppm] = s_energy[0];
    }
  }

}

__global__ void sum_energy_kernel2(int vflag)
{
  ENERGY_CFLOAT myenergy = energy[threadIdx.x * ny_pppm];
  ENERGY_CFLOAT* s_energy = (ENERGY_CFLOAT*) sharedmem;
  s_energy[threadIdx.x] = myenergy;
  __syncthreads();
  reduceBlock(s_energy);

  if(threadIdx.x == 0)
    energy[0] = s_energy[0];

  if(vflag) {
    __syncthreads();

    for(int j = 0; j < 6; j++) {
      myenergy = virial[threadIdx.x * ny_pppm + j * ny_pppm * nz_pppm];
      s_energy[threadIdx.x] = myenergy;
      __syncthreads();
      reduceBlock(s_energy);

      if(threadIdx.x == 0)
        virial[j] = s_energy[0];
    }
  }
}

__device__ PPPM_CFLOAT rho1d(int k, PPPM_CFLOAT d, PPPM_CFLOAT* srho_coeff)
{
  PPPM_CFLOAT rho1d_tmp = PPPM_F(0.0);

  for(int l = order - 1; l >= 0; l--)
    rho1d_tmp = srho_coeff[l * order + k - (1 - order) / 2] + rho1d_tmp * d;

  return rho1d_tmp;
}

__global__ void particle_map_kernel(int* flag)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;

  if(i < nlocal) {
    int nx, ny, nz;
    PPPM_CFLOAT shift = PPPM_F(0.5) - shiftone; //+OFFSET;
    nx = (int)((_x[i] - _boxlo[0]) * delxinv + shift); // - OFFSET;
    ny = (int)((_x[i + nmax] - _boxlo[1]) * delyinv + shift); // - OFFSET;
    nz = (int)((_x[i + 2 * nmax] - _boxlo[2]) * delzinv + shift); // - OFFSET;

    part2grid[i] = nx;
    part2grid[i + nmax] = ny;
    part2grid[i + 2 * nmax] = nz;

    // check that entire stencil around nx,ny,nz will fit in my 3d brick
    if(nx + nlower < nxlo_out || nx + nupper > nxhi_out ||
        ny + nlower < nylo_out || ny + nupper > nyhi_out ||
        nz + nlower < nzlo_out || nz + nupper > nzhi_out) {
      flag[0]++;
      debugdata[0] = i;
      debugdata[1] = _boxlo[0];
      debugdata[2] = _boxlo[1];
      debugdata[3] = _boxlo[2];
      debugdata[4] = nx;
      debugdata[5] = ny;
      debugdata[6] = nz;
      debugdata[7] = _x[i];
      debugdata[8] = _x[i + _nmax];
      debugdata[9] = _x[i + 2 * _nmax];
      debugdata[10] = nlocal;

    }
  }
}

__global__ void make_rho_kernelA()
{
  int i, l, m, n, nx, ny, nz, mx, my, mz;

  // clear 3d density array


  // loop over my charges, add their contribution to nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt

  i = blockIdx.x * blockDim.x + threadIdx.x;

  if(i < nlocal) {

    PPPM_CFLOAT dx, dy, dz, x0, y0, z0;
    nx = part2grid[i];
    ny = part2grid[i + nmax];
    nz = part2grid[i + 2 * nmax];
    dx = nx + shiftone - (_x[i] - _boxlo[0]) * delxinv;
    dy = ny + shiftone - (_x[i + nmax] - _boxlo[1]) * delyinv;
    dz = nz + shiftone - (_x[i + 2 * nmax] - _boxlo[2]) * delzinv;

    z0 = delxinv * delyinv * delzinv * _q[i];

    for(n = nlower; n <= nupper; n++) {
      mz = n + nz;
      y0 = z0 * rho1d(n, dz, rho_coeff);

      for(m = nlower; m <= nupper; m++) {
        my = m + ny;
        x0 = y0 * rho1d(m, dy, rho_coeff);

        for(l = nlower; l <= nupper; l++) {
          mx = l + nx;
          int mzyx = ((mz - nzlo_out) * (nyhi_out - nylo_out + 1) + my - nylo_out) * (nxhi_out - nxlo_out + 1) + mx - nxlo_out;

          while(atomicAdd(&density_brick_int[mzyx], 1) != 0) atomicAdd(&density_brick_int[mzyx], -1);

          density_brick[mzyx] += x0 * rho1d(l, dx, rho_coeff);
          __threadfence();
          atomicAdd(&density_brick_int[mzyx], -1);
          __syncthreads();

        }
      }
    }
  }
}

__global__ void make_rho_kernel(int* flag, int read_threads_at_same_time)
{
  int i, l, m, n, nx, ny, nz, mx, my, mz, a, b;

  // clear 3d density array


  // loop over my charges, add their contribution to nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // int nzxy=blockIdx.x*gridDim.y+blockIdx.y;

  int nelements = nupper - nlower + 1;
  int* idx = (int*) sharedmem;
  int* sdensity_brick_int = &idx[blockDim.x];
  PPPM_CFLOAT* srho_coeff = (PPPM_CFLOAT*) &sdensity_brick_int[nelements * blockDim.x];

  if(threadIdx.x < order * (order / 2 - (1 - order) / 2 + 1))
    srho_coeff[threadIdx.x] = rho_coeff[threadIdx.x];

  __syncthreads();

  i = blockIdx.x * blockDim.x + threadIdx.x;

  if(false) {
    if(i < nlocal) {

      PPPM_CFLOAT dx, dy, dz, x0, y0, z0;
      nx = part2grid[i];
      ny = part2grid[i + nmax];
      nz = part2grid[i + 2 * nmax];
      dx = nx + shiftone - (_x[i] - _boxlo[0]) * delxinv;
      dy = ny + shiftone - (_x[i + nmax] - _boxlo[1]) * delyinv;
      dz = nz + shiftone - (_x[i + 2 * nmax] - _boxlo[2]) * delzinv;

      z0 = delxinv * delyinv * delzinv * _q[i];

      for(n = nlower; n <= nupper; n++) {
        mz = n + nz;
        y0 = z0 * rho1d(n, dz, srho_coeff);

        for(m = nlower; m <= nupper; m++) {
          my = m + ny;
          x0 = y0 * rho1d(m, dy, srho_coeff);

          for(l = nlower; l <= nupper; l++) {
            mx = l + nx;
            int mzyx = ((mz - nzlo_out) * (nyhi_out - nylo_out + 1) + my - nylo_out) * (nxhi_out - nxlo_out + 1) + mx - nxlo_out;

            a = int(x0 * rho1d(l, dx, srho_coeff) * density_intScale);
            b = (atomicAdd(&density_brick_int[mzyx], a) | a);

            if(((b) & (0x7c000000)) && (not((b) & (0x80000000)))) {
              flag[1]++;

              if((b) & (0x60000000)) flag[0]++;
            }

            __syncthreads();
          }
        }
      }
    }

    return;
  }

  i = blockIdx.x * blockDim.x + threadIdx.x;
  {

    PPPM_CFLOAT dx, dy, dz, x0, y0, z0, qtmp;

    if(i < nlocal) {
      qtmp = _q[i];
      nx = part2grid[i];
      ny = part2grid[i + nmax];
      nz = part2grid[i + 2 * nmax];
      dx = nx + shiftone - (_x[i] - _boxlo[0]) * delxinv;
      dy = ny + shiftone - (_x[i + nmax] - _boxlo[1]) * delyinv;
      dz = nz + shiftone - (_x[i + 2 * nmax] - _boxlo[2]) * delzinv;
      z0 = delxinv * delyinv * delzinv * qtmp;
    } else {
      nx = ny = nz = 1;
      dx = dy = dz = PPPM_F(0.1);
    }

    __syncthreads();

    for(n = nlower; n <= nupper; n++) {
      mz = n + nz;
      y0 = z0 * rho1d(n, dz, srho_coeff);

      for(m = nlower; m <= nupper; m++) {
        my = m + ny;
        x0 = y0 * rho1d(m, dy, srho_coeff);

        if(i < nlocal) {
          idx[threadIdx.x] = ((mz - nzlo_out) * (nyhi_out - nylo_out + 1) + my - nylo_out) * (nxhi_out - nxlo_out + 1) + nx + nlower - nxlo_out;

          for(l = nlower; l <= nupper; l++) {
            sdensity_brick_int[threadIdx.x * nelements + l - nlower] = int(x0 * rho1d(l, dx, srho_coeff) * density_intScale);
          }
        } else idx[threadIdx.x] = -1;

        __syncthreads();

        for(int ii = 0; ii < blockDim.x; ii += read_threads_at_same_time) {
          int kk = threadIdx.x / nelements;

          if((threadIdx.x < nelements * read_threads_at_same_time) && (kk + ii < blockDim.x) && (idx[ii + kk] > -1)) {
            a = sdensity_brick_int[ii * nelements + threadIdx.x];
            //if(a*a>1e-100)
            b = (atomicAdd(&density_brick_int[idx[ii + kk] + threadIdx.x - kk * nelements], a) | a);

            //else
            //b=(density_brick_int[idx[ii+kk]+threadIdx.x-kk*nelements]|a);
            if(((b) & (0x7c000000)) && (not((b) & (0x80000000)))) {
              flag[1]++;

              if((b) & (0x60000000)) flag[0]++;
            }
          }
        }

        __syncthreads();	   //*/
      }
    }

  }
}

__global__ void scale_rho_kernel()
{
  int i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x;
  density_brick[i] = (1.0 / density_intScale) * density_brick_int[i];
}

__global__ void fieldforce_kernel(int elements_per_thread, int read_threads_at_same_time, int* flag) //20*x64 0.36
{
  int i;

  // loop over my charges, interpolate electric field from nearby grid points
  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt
  // ek = 3 components of E-field on particle
  i = blockIdx.x * blockDim.x + threadIdx.x;
  int* idx = (int*) sharedmem;
  PPPM_CFLOAT* tmp_brick = (PPPM_CFLOAT*) &idx[blockDim.x];
  PPPM_CFLOAT* srho_coeff = (PPPM_CFLOAT*) &tmp_brick[3 * blockDim.x * elements_per_thread];

  if(threadIdx.x < order * (order / 2 - (1 - order) / 2 + 1))
    srho_coeff[threadIdx.x] = rho_coeff[threadIdx.x];

  __syncthreads();
  {
    int l, m, n, nx, ny, nz, my, mz;
    PPPM_CFLOAT dx, dy, dz, x0, y0, z0;
    PPPM_CFLOAT ek[3];

    if(i < nlocal) {
      nx = part2grid[i];
      ny = part2grid[i + nmax];
      nz = part2grid[i + 2 * nmax];
      dx = nx + shiftone - (_x[i] - _boxlo[0]) * delxinv;
      dy = ny + shiftone - (_x[i + nmax] - _boxlo[1]) * delyinv;
      dz = nz + shiftone - (_x[i + 2 * nmax] - _boxlo[2]) * delzinv;

      ek[0] = ek[1] = ek[2] = PPPM_F(0.0);
    } else {
      nx = ny = nz = 1;
      dx = dy = dz = PPPM_F(0.1);
    }

    __syncthreads();

    for(n = nlower; n <= nupper; n++) {
      mz = n + nz;
      z0 = rho1d(n, dz, srho_coeff);

      for(m = nlower; m <= nupper; m++) {
        my = m + ny;
        y0 = z0 * rho1d(m, dy, srho_coeff);


        if(i < nlocal)
          idx[threadIdx.x] = ((mz - nzlo_out) * (nyhi_out - nylo_out + 1) + my - nylo_out) * (nxhi_out - nxlo_out + 1) + nx + nlower - nxlo_out;
        else idx[threadIdx.x] = -1;

        __syncthreads();

        for(int ii = 0; ii < blockDim.x; ii += read_threads_at_same_time) {
          int kk = threadIdx.x / elements_per_thread;

          if((threadIdx.x < elements_per_thread * read_threads_at_same_time) && (kk + ii < blockDim.x) && (idx[ii + kk] > -1)) {
            tmp_brick[ii * elements_per_thread + threadIdx.x] = vdx_brick[idx[ii + kk] + threadIdx.x - kk * elements_per_thread];
            tmp_brick[(ii + blockDim.x)*elements_per_thread + threadIdx.x] = vdy_brick[idx[ii + kk] + threadIdx.x - kk * elements_per_thread];
            tmp_brick[(ii + 2 * blockDim.x)*elements_per_thread + threadIdx.x] = vdz_brick[idx[ii + kk] + threadIdx.x - kk * elements_per_thread];
          }
        }

        __syncthreads();

        if(i < nlocal)
          for(l = nlower; l <= nupper; l++) {
            x0 = y0 * rho1d(l, dx, srho_coeff);

            ek[0] -= x0 * tmp_brick[threadIdx.x * elements_per_thread + l - nlower];
            ek[1] -= x0 * tmp_brick[threadIdx.x * elements_per_thread + l - nlower + blockDim.x * elements_per_thread];
            ek[2] -= x0 * tmp_brick[threadIdx.x * elements_per_thread + l - nlower + 2 * blockDim.x * elements_per_thread];
          }

        __syncthreads();
      }
    }

    // convert E-field to force


    _f[i] += qqrd2e * _q[i] * ek[0];
    _f[i + nmax] += qqrd2e * _q[i] * ek[1];
    _f[i + 2 * nmax] += qqrd2e * _q[i] * ek[2];
  }
}

__global__ void slabcorr_energy_kernel(ENERGY_CFLOAT* buf)
{
  ENERGY_CFLOAT* dipole = (ENERGY_CFLOAT*) sharedmem;
  int i = blockIdx.x * blockDim.x + threadIdx.x;

  if(i < nlocal)
    dipole[threadIdx.x] = _q[i] * _x[i + 2 * nmax];
  else
    dipole[threadIdx.x] = ENERGY_F(0.0);

  __syncthreads();
  reduceBlock(dipole);

  if(threadIdx.x == 0) buf[blockIdx.x] = dipole[0];
}

__global__ void slabcorr_force_kernel(F_CFLOAT ffact)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;

  if(i < nlocal)
    _f[i + 2 * nmax] += qqrd2e * _q[i] * ffact;
}


__global__ void initfftdata_core_kernel(PPPM_CFLOAT* in, FFT_CFLOAT* out)
{
  out[2 * (((blockIdx.x * (nyhi_in - nylo_in + 1) + blockIdx.y) * (nxhi_in - nxlo_in + 1)) + threadIdx.x)] = in[(((blockIdx.x + nzlo_in - nzlo_out) * (nyhi_out - nylo_out + 1) + blockIdx.y + nylo_in - nylo_out) * (nxhi_out - nxlo_out + 1)) + threadIdx.x + nxlo_in - nxlo_out];
  out[2 * (((blockIdx.x * (nyhi_in - nylo_in + 1) + blockIdx.y) * (nxhi_in - nxlo_in + 1)) + threadIdx.x) + 1] = 0;
}

__global__ void initfftdata_z_kernel(PPPM_CFLOAT* in, FFT_CFLOAT* out)
{
  if(slabflag) {
    if(blockIdx.x < nzlo_in - nzlo_out)
      out[2 * (((nzhi_in - nzlo_in + 2 - nupper - slabflag + blockIdx.x) * (nyhi_in - nylo_in + 1) + blockIdx.y) * (nxhi_in - nxlo_in + 1) + threadIdx.x)] += in[(((blockIdx.x) * (nyhi_out - nylo_out + 1) + blockIdx.y + nylo_in - nylo_out) * (nxhi_out - nxlo_out + 1)) + threadIdx.x + nxlo_in - nxlo_out];
  } else {
    if(blockIdx.x < nzlo_in - nzlo_out)
      out[2 * (((blockIdx.x + 2 * (nzhi_in + 1) - nzlo_in - nzhi_out) * (nyhi_in - nylo_in + 1) + blockIdx.y) * (nxhi_in - nxlo_in + 1) + threadIdx.x)] += in[(((blockIdx.x) * (nyhi_out - nylo_out + 1) + blockIdx.y + nylo_in - nylo_out) * (nxhi_out - nxlo_out + 1)) + threadIdx.x + nxlo_in - nxlo_out];
  }

  if(blockIdx.x < nzhi_out - nzhi_in)
    out[2 * ((((blockIdx.x) * (nyhi_in - nylo_in + 1) + blockIdx.y) * (nxhi_in - nxlo_in + 1)) + threadIdx.x)] += in[(((blockIdx.x + (nzhi_out - nzlo_in)) * (nyhi_out - nylo_out + 1) + blockIdx.y + nylo_in - nylo_out) * (nxhi_out - nxlo_out + 1)) + threadIdx.x + nxlo_in - nxlo_out];
}

__global__ void initfftdata_y_kernel(PPPM_CFLOAT* in, FFT_CFLOAT* out)
{
  if(blockIdx.y < nylo_in - nylo_out)
    out[2 * (((blockIdx.x * (nyhi_in - nylo_in + 1) + (2 * (nyhi_in + 1) - nylo_in - nyhi_out) + blockIdx.y) * (nxhi_in - nxlo_in + 1)) + threadIdx.x)] += in[(((blockIdx.x + nzlo_in - nzlo_out) * (nyhi_out - nylo_out + 1) + blockIdx.y) * (nxhi_out - nxlo_out + 1)) + threadIdx.x + nxlo_in - nxlo_out];

  if(blockIdx.y < nyhi_out - nyhi_in)
    out[2 * (((blockIdx.x * (nyhi_in - nylo_in + 1) + blockIdx.y) * (nxhi_in - nxlo_in + 1)) + threadIdx.x)] += in[(((blockIdx.x + nzlo_in - nzlo_out) * (nyhi_out - nylo_out + 1) + blockIdx.y + (nyhi_out - nylo_in)) * (nxhi_out - nxlo_out + 1)) + threadIdx.x + nxlo_in - nxlo_out];
}

__global__ void initfftdata_x_kernel(PPPM_CFLOAT* in, FFT_CFLOAT* out)
{
  if(threadIdx.x < nxlo_in - nxlo_out)
    out[2 * (((blockIdx.x * (nyhi_in - nylo_in + 1) + blockIdx.y) * (nxhi_in - nxlo_in + 1)) + threadIdx.x + 2 * (nxhi_in + 1) - nxlo_in - nxhi_out)] += in[(((blockIdx.x + nzlo_in - nzlo_out) * (nyhi_out - nylo_out + 1) + blockIdx.y + nylo_in - nylo_out) * (nxhi_out - nxlo_out + 1)) + threadIdx.x];

  if(threadIdx.x < nxhi_out - nxhi_in)
    out[2 * (((blockIdx.x * (nyhi_in - nylo_in + 1) + blockIdx.y) * (nxhi_in - nxlo_in + 1)) + threadIdx.x)] += in[(((blockIdx.x + nzlo_in - nzlo_out) * (nyhi_out - nylo_out + 1) + blockIdx.y + nylo_in - nylo_out) * (nxhi_out - nxlo_out + 1)) + threadIdx.x + nxhi_in - nxlo_out + 1];
}

__global__ void initfftdata_yz_kernel(PPPM_CFLOAT* in, FFT_CFLOAT* out)
{
  if(slabflag) {
    if(blockIdx.x < nzlo_in - nzlo_out)
      if(blockIdx.y < nyhi_out - nyhi_in)
        out[2 * ((((nzhi_in - nzlo_in + 2 - nupper - slabflag + blockIdx.x) * (nyhi_in - nylo_in + 1) + blockIdx.y) * (nxhi_in - nxlo_in + 1)) + threadIdx.x)] += in[(((blockIdx.x) * (nyhi_out - nylo_out + 1) + blockIdx.y + nyhi_in - nylo_out + 1) * (nxhi_out - nxlo_out + 1)) + threadIdx.x + nxlo_in - nxlo_out];

    if(blockIdx.x < nzlo_in - nzlo_out)
      if(blockIdx.y < nylo_in - nylo_out)
        out[2 * ((((nzhi_in - nzlo_in + 2 - nupper - slabflag + blockIdx.x) * (nyhi_in - nylo_in + 1) + blockIdx.y + 2 * (nyhi_in + 1) - nylo_in - nyhi_out) * (nxhi_in - nxlo_in + 1)) + threadIdx.x)] += in[(((blockIdx.x) * (nyhi_out - nylo_out + 1) + blockIdx.y) * (nxhi_out - nxlo_out + 1)) + threadIdx.x + nxlo_in - nxlo_out];
  } else {
    if(blockIdx.x < nzlo_in - nzlo_out)
      if(blockIdx.y < nyhi_out - nyhi_in)
        out[2 * ((((blockIdx.x + 2 * (nzhi_in + 1) - nzlo_in - nzhi_out) * (nyhi_in - nylo_in + 1) + blockIdx.y) * (nxhi_in - nxlo_in + 1)) + threadIdx.x)] += in[(((blockIdx.x) * (nyhi_out - nylo_out + 1) + blockIdx.y + nyhi_in - nylo_out + 1) * (nxhi_out - nxlo_out + 1)) + threadIdx.x + nxlo_in - nxlo_out];

    if(blockIdx.x < nzlo_in - nzlo_out)
      if(blockIdx.y < nylo_in - nylo_out)
        out[2 * ((((blockIdx.x + 2 * (nzhi_in + 1) - nzlo_in - nzhi_out) * (nyhi_in - nylo_in + 1) + blockIdx.y + 2 * (nyhi_in + 1) - nylo_in - nyhi_out) * (nxhi_in - nxlo_in + 1)) + threadIdx.x)] += in[(((blockIdx.x) * (nyhi_out - nylo_out + 1) + blockIdx.y) * (nxhi_out - nxlo_out + 1)) + threadIdx.x + nxlo_in - nxlo_out];
  }

  if(blockIdx.x < nzhi_out - nzhi_in)
    if(blockIdx.y < nyhi_out - nyhi_in)
      out[2 * (((blockIdx.x * (nyhi_in - nylo_in + 1) + blockIdx.y) * (nxhi_in - nxlo_in + 1)) + threadIdx.x)] += in[(((blockIdx.x + nzhi_in - nzlo_out + 1) * (nyhi_out - nylo_out + 1) + blockIdx.y + nyhi_in - nylo_out + 1) * (nxhi_out - nxlo_out + 1)) + threadIdx.x + nxlo_in - nxlo_out];

  if(blockIdx.x < nzhi_out - nzhi_in)
    if(blockIdx.y < nylo_in - nylo_out)
      out[2 * (((blockIdx.x * (nyhi_in - nylo_in + 1) + blockIdx.y + 2 * (nyhi_in + 1) - nylo_in - nyhi_out) * (nxhi_in - nxlo_in + 1)) + threadIdx.x)] += in[(((blockIdx.x + nzhi_in - nzlo_out + 1) * (nyhi_out - nylo_out + 1) + blockIdx.y) * (nxhi_out - nxlo_out + 1)) + threadIdx.x + nxlo_in - nxlo_out];
}

__global__ void initfftdata_xz_kernel(PPPM_CFLOAT* in, FFT_CFLOAT* out)
{
  if(blockIdx.x < nzhi_out - nzhi_in)
    if(threadIdx.x < nxlo_in - nxlo_out)
      out[2 * (((blockIdx.x * (nyhi_in - nylo_in + 1) + blockIdx.y) * (nxhi_in - nxlo_in + 1)) + threadIdx.x + 2 * (nxhi_in + 1) - nxlo_in - nxhi_out)] += in[(((blockIdx.x + nzhi_in - nzlo_out + 1) * (nyhi_out - nylo_out + 1) + blockIdx.y + nylo_in - nylo_out) * (nxhi_out - nxlo_out + 1)) + threadIdx.x];

  if(blockIdx.x < nzhi_out - nzhi_in)
    if(threadIdx.x < nxhi_out - nxhi_in)
      out[2 * (((blockIdx.x * (nyhi_in - nylo_in + 1) + blockIdx.y) * (nxhi_in - nxlo_in + 1)) + threadIdx.x)] += in[(((blockIdx.x + nzhi_in - nzlo_out + 1) * (nyhi_out - nylo_out + 1) + blockIdx.y + nylo_in - nylo_out) * (nxhi_out - nxlo_out + 1)) + threadIdx.x + nxhi_in - nxlo_out + 1];

  if(slabflag) {
    if(blockIdx.x < nzlo_in - nzlo_out)
      if(threadIdx.x < nxlo_in - nxlo_out)
        out[2 * ((((nzhi_in - nzlo_in + 2 - nupper - slabflag + blockIdx.x) * (nyhi_in - nylo_in + 1) + blockIdx.y) * (nxhi_in - nxlo_in + 1)) + threadIdx.x + 2 * (nxhi_in + 1) - nxlo_in - nxhi_out)] += in[(((blockIdx.x) * (nyhi_out - nylo_out + 1) + blockIdx.y + nylo_in - nylo_out) * (nxhi_out - nxlo_out + 1)) + threadIdx.x];

    if(blockIdx.x < nzlo_in - nzlo_out)
      if(threadIdx.x < nxhi_out - nxhi_in)
        out[2 * ((((nzhi_in - nzlo_in + 2 - nupper - slabflag + blockIdx.x) * (nyhi_in - nylo_in + 1) + blockIdx.y) * (nxhi_in - nxlo_in + 1)) + threadIdx.x)] += in[(((blockIdx.x) * (nyhi_out - nylo_out + 1) + blockIdx.y + nylo_in - nylo_out) * (nxhi_out - nxlo_out + 1)) + threadIdx.x + nxhi_in - nxlo_out + 1];
  } else {
    if(blockIdx.x < nzlo_in - nzlo_out)
      if(threadIdx.x < nxlo_in - nxlo_out)
        out[2 * ((((blockIdx.x + 2 * (nzhi_in + 1) - nzlo_in - nzhi_out) * (nyhi_in - nylo_in + 1) + blockIdx.y) * (nxhi_in - nxlo_in + 1)) + threadIdx.x + 2 * (nxhi_in + 1) - nxlo_in - nxhi_out)] += in[(((blockIdx.x) * (nyhi_out - nylo_out + 1) + blockIdx.y + nylo_in - nylo_out) * (nxhi_out - nxlo_out + 1)) + threadIdx.x];

    if(blockIdx.x < nzlo_in - nzlo_out)
      if(threadIdx.x < nxhi_out - nxhi_in)
        out[2 * ((((blockIdx.x + 2 * (nzhi_in + 1) - nzlo_in - nzhi_out) * (nyhi_in - nylo_in + 1) + blockIdx.y) * (nxhi_in - nxlo_in + 1)) + threadIdx.x)] += in[(((blockIdx.x) * (nyhi_out - nylo_out + 1) + blockIdx.y + nylo_in - nylo_out) * (nxhi_out - nxlo_out + 1)) + threadIdx.x + nxhi_in - nxlo_out + 1];
  }
}

__global__ void initfftdata_xy_kernel(PPPM_CFLOAT* in, FFT_CFLOAT* out)
{
  if(blockIdx.y < nyhi_out - nyhi_in)
    if(threadIdx.x < nxlo_in - nxlo_out)
      out[2 * (((blockIdx.x * (nyhi_in - nylo_in + 1) + blockIdx.y) * (nxhi_in - nxlo_in + 1)) + threadIdx.x + 2 * (nxhi_in + 1) - nxlo_in - nxhi_out)] += in[(((blockIdx.x + nzlo_in - nzlo_out) * (nyhi_out - nylo_out + 1) + blockIdx.y + nyhi_in - nylo_out + 1) * (nxhi_out - nxlo_out + 1)) + threadIdx.x];

  if(blockIdx.y < nyhi_out - nyhi_in)
    if(threadIdx.x < nxhi_out - nxhi_in)
      out[2 * (((blockIdx.x * (nyhi_in - nylo_in + 1) + blockIdx.y) * (nxhi_in - nxlo_in + 1)) + threadIdx.x)] += in[(((blockIdx.x + nzlo_in - nzlo_out) * (nyhi_out - nylo_out + 1) + blockIdx.y + nyhi_in - nylo_out + 1) * (nxhi_out - nxlo_out + 1)) + threadIdx.x + nxhi_in - nxlo_out + 1];

  if(blockIdx.y < nylo_in - nylo_out)
    if(threadIdx.x < nxlo_in - nxlo_out)
      out[2 * (((blockIdx.x * (nyhi_in - nylo_in + 1) + blockIdx.y + 2 * (nyhi_in + 1) - nylo_in - nyhi_out) * (nxhi_in - nxlo_in + 1)) + threadIdx.x + 2 * (nxhi_in + 1) - nxlo_in - nxhi_out)] += in[(((blockIdx.x + nzlo_in - nzlo_out) * (nyhi_out - nylo_out + 1) + blockIdx.y) * (nxhi_out - nxlo_out + 1)) + threadIdx.x];

  if(blockIdx.y < nylo_in - nylo_out)
    if(threadIdx.x < nxhi_out - nxhi_in)
      out[2 * (((blockIdx.x * (nyhi_in - nylo_in + 1) + blockIdx.y + 2 * (nyhi_in + 1) - nylo_in - nyhi_out) * (nxhi_in - nxlo_in + 1)) + threadIdx.x)] += in[(((blockIdx.x + nzlo_in - nzlo_out) * (nyhi_out - nylo_out + 1) + blockIdx.y) * (nxhi_out - nxlo_out + 1)) + threadIdx.x + nxhi_in - nxlo_out + 1];
}

__global__ void initfftdata_xyz_kernel(PPPM_CFLOAT* in, FFT_CFLOAT* out)
{
  if(blockIdx.x < nzhi_out - nzhi_in)
    if(blockIdx.y < nyhi_out - nyhi_in)
      if(threadIdx.x < nxlo_in - nxlo_out)
        out[2 * (((blockIdx.x * (nyhi_in - nylo_in + 1) + blockIdx.y) * (nxhi_in - nxlo_in + 1)) + threadIdx.x + 2 * (nxhi_in + 1) - nxlo_in - nxhi_out)] += in[(((blockIdx.x + nzhi_in - nzlo_out + 1) * (nyhi_out - nylo_out + 1) + blockIdx.y + nyhi_in - nylo_out + 1) * (nxhi_out - nxlo_out + 1)) + threadIdx.x];

  if(blockIdx.x < nzhi_out - nzhi_in)
    if(blockIdx.y < nyhi_out - nyhi_in)
      if(threadIdx.x < nxhi_out - nxhi_in)
        out[2 * (((blockIdx.x * (nyhi_in - nylo_in + 1) + blockIdx.y) * (nxhi_in - nxlo_in + 1)) + threadIdx.x)] += in[(((blockIdx.x + nzhi_in - nzlo_out + 1) * (nyhi_out - nylo_out + 1) + blockIdx.y + nyhi_in - nylo_out + 1) * (nxhi_out - nxlo_out + 1)) + threadIdx.x + nxhi_in - nxlo_out + 1];

  if(blockIdx.x < nzhi_out - nzhi_in)
    if(blockIdx.y < nylo_in - nylo_out)
      if(threadIdx.x < nxlo_in - nxlo_out)
        out[2 * (((blockIdx.x * (nyhi_in - nylo_in + 1) + blockIdx.y + 2 * (nyhi_in + 1) - nylo_in - nyhi_out) * (nxhi_in - nxlo_in + 1)) + threadIdx.x + 2 * (nxhi_in + 1) - nxlo_in - nxhi_out)] += in[(((blockIdx.x + nzhi_in - nzlo_out + 1) * (nyhi_out - nylo_out + 1) + blockIdx.y) * (nxhi_out - nxlo_out + 1)) + threadIdx.x];

  if(blockIdx.x < nzhi_out - nzhi_in)
    if(blockIdx.y < nylo_in - nylo_out)
      if(threadIdx.x < nxhi_out - nxhi_in)
        out[2 * (((blockIdx.x * (nyhi_in - nylo_in + 1) + blockIdx.y + 2 * (nyhi_in + 1) - nylo_in - nyhi_out) * (nxhi_in - nxlo_in + 1)) + threadIdx.x)] += in[(((blockIdx.x + nzhi_in - nzlo_out + 1) * (nyhi_out - nylo_out + 1) + blockIdx.y) * (nxhi_out - nxlo_out + 1)) + threadIdx.x + nxhi_in - nxlo_out + 1];

  if(slabflag) {
    if(blockIdx.x < nzlo_in - nzlo_out)
      if(blockIdx.y < nyhi_out - nyhi_in)
        if(threadIdx.x < nxlo_in - nxlo_out)
          out[2 * ((((nzhi_in - nzlo_in + 2 - nupper - slabflag + blockIdx.x) * (nyhi_in - nylo_in + 1) + blockIdx.y) * (nxhi_in - nxlo_in + 1)) + threadIdx.x + 2 * (nxhi_in + 1) - nxlo_in - nxhi_out)] += in[(((blockIdx.x) * (nyhi_out - nylo_out + 1) + blockIdx.y + nyhi_in - nylo_out + 1) * (nxhi_out - nxlo_out + 1)) + threadIdx.x];

    if(blockIdx.x < nzlo_in - nzlo_out)
      if(blockIdx.y < nyhi_out - nyhi_in)
        if(threadIdx.x < nxhi_out - nxhi_in)
          out[2 * ((((nzhi_in - nzlo_in + 2 - nupper - slabflag + blockIdx.x) * (nyhi_in - nylo_in + 1) + blockIdx.y) * (nxhi_in - nxlo_in + 1)) + threadIdx.x)] += in[(((blockIdx.x) * (nyhi_out - nylo_out + 1) + blockIdx.y + nyhi_in - nylo_out + 1) * (nxhi_out - nxlo_out + 1)) + threadIdx.x + nxhi_in - nxlo_out + 1];

    if(blockIdx.x < nzlo_in - nzlo_out)
      if(blockIdx.y < nylo_in - nylo_out)
        if(threadIdx.x < nxlo_in - nxlo_out)
          out[2 * ((((nzhi_in - nzlo_in + 2 - nupper - slabflag + blockIdx.x) * (nyhi_in - nylo_in + 1) + blockIdx.y + 2 * (nyhi_in + 1) - nylo_in - nyhi_out) * (nxhi_in - nxlo_in + 1)) + threadIdx.x + 2 * (nxhi_in + 1) - nxlo_in - nxhi_out)] += in[(((blockIdx.x) * (nyhi_out - nylo_out + 1) + blockIdx.y) * (nxhi_out - nxlo_out + 1)) + threadIdx.x];

    if(blockIdx.x < nzlo_in - nzlo_out)
      if(blockIdx.y < nylo_in - nylo_out)
        if(threadIdx.x < nxhi_out - nxhi_in)
          out[2 * ((((nzhi_in - nzlo_in + 2 - nupper - slabflag + blockIdx.x) * (nyhi_in - nylo_in + 1) + blockIdx.y + 2 * (nyhi_in + 1) - nylo_in - nyhi_out) * (nxhi_in - nxlo_in + 1)) + threadIdx.x)] += in[(((blockIdx.x) * (nyhi_out - nylo_out + 1) + blockIdx.y) * (nxhi_out - nxlo_out + 1)) + threadIdx.x + nxhi_in - nxlo_out + 1];
  } else {
    if(blockIdx.x < nzlo_in - nzlo_out)
      if(blockIdx.y < nyhi_out - nyhi_in)
        if(threadIdx.x < nxlo_in - nxlo_out)
          out[2 * ((((blockIdx.x + 2 * (nzhi_in + 1) - nzlo_in - nzhi_out) * (nyhi_in - nylo_in + 1) + blockIdx.y) * (nxhi_in - nxlo_in + 1)) + threadIdx.x + 2 * (nxhi_in + 1) - nxlo_in - nxhi_out)] += in[(((blockIdx.x) * (nyhi_out - nylo_out + 1) + blockIdx.y + nyhi_in - nylo_out + 1) * (nxhi_out - nxlo_out + 1)) + threadIdx.x];

    if(blockIdx.x < nzlo_in - nzlo_out)
      if(blockIdx.y < nyhi_out - nyhi_in)
        if(threadIdx.x < nxhi_out - nxhi_in)
          out[2 * ((((blockIdx.x + 2 * (nzhi_in + 1) - nzlo_in - nzhi_out) * (nyhi_in - nylo_in + 1) + blockIdx.y) * (nxhi_in - nxlo_in + 1)) + threadIdx.x)] += in[(((blockIdx.x) * (nyhi_out - nylo_out + 1) + blockIdx.y + nyhi_in - nylo_out + 1) * (nxhi_out - nxlo_out + 1)) + threadIdx.x + nxhi_in - nxlo_out + 1];

    if(blockIdx.x < nzlo_in - nzlo_out)
      if(blockIdx.y < nylo_in - nylo_out)
        if(threadIdx.x < nxlo_in - nxlo_out)
          out[2 * ((((blockIdx.x + 2 * (nzhi_in + 1) - nzlo_in - nzhi_out) * (nyhi_in - nylo_in + 1) + blockIdx.y + 2 * (nyhi_in + 1) - nylo_in - nyhi_out) * (nxhi_in - nxlo_in + 1)) + threadIdx.x + 2 * (nxhi_in + 1) - nxlo_in - nxhi_out)] += in[(((blockIdx.x) * (nyhi_out - nylo_out + 1) + blockIdx.y) * (nxhi_out - nxlo_out + 1)) + threadIdx.x];

    if(blockIdx.x < nzlo_in - nzlo_out)
      if(blockIdx.y < nylo_in - nylo_out)
        if(threadIdx.x < nxhi_out - nxhi_in)
          out[2 * ((((blockIdx.x + 2 * (nzhi_in + 1) - nzlo_in - nzhi_out) * (nyhi_in - nylo_in + 1) + blockIdx.y + 2 * (nyhi_in + 1) - nylo_in - nyhi_out) * (nxhi_in - nxlo_in + 1)) + threadIdx.x)] += in[(((blockIdx.x) * (nyhi_out - nylo_out + 1) + blockIdx.y) * (nxhi_out - nxlo_out + 1)) + threadIdx.x + nxhi_in - nxlo_out + 1];
  }
}
