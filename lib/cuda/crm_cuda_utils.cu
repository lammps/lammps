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

#ifndef CRM_CUDA_UTILS
#define CRM_CUDA_UTILS

//split n threads into 2 dimensional grid + threads, return values are grid.x grid.y and threads.x
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

inline int3 getgrid(int n, int shared_per_thread = 0, int threadsmax = 256, bool p2 = false)
{
  int3 gridparams;
  int sharedsize = 16000;

  if(shared_per_thread > 0) threadsmax = sharedsize / shared_per_thread < threadsmax ? sharedsize / shared_per_thread : threadsmax;

  if((n < 60 * 32) || (threadsmax < 64))
    gridparams.z = 32;
  else if((n < 60 * 64) || (threadsmax < 128))
    gridparams.z = 64;
  else if((n < 60 * 128) || (threadsmax < 256))
    gridparams.z = 128;
  else if((n < 60 * 256) || (threadsmax < 512))
    gridparams.z = 256;
  else gridparams.z = 512;

  if(p2) {
    gridparams.z = 16;

    while(gridparams.z * 2 <= threadsmax) gridparams.z *= 2;
  }


  int blocks = (n + gridparams.z - 1) / gridparams.z;

  if(blocks > 10000)
    gridparams.x = gridparams.y = int(sqrt(blocks));
  else {
    gridparams.x = blocks;
    gridparams.y = 1;
  }

  while(gridparams.x * gridparams.y * gridparams.z < n) gridparams.x++;

  if(gridparams.x == 0) gridparams.x = 1;

  return gridparams;
}

//return value: 1 if f<0; else: 0
//take care if working with values as "blockId.x-n" for f: it might be interpreted as a unsigned int
static inline __device__ int negativCUDA(float f)
{
  return ((unsigned int)1 << 31 & (__float_as_int(f))) >> 31;
}

//return value: -1 if f<0; else +1
static inline __device__ float fsignCUDA(float f)
{
  return f < 0.0f ? -1.0f : 1.0f;
}

//functions to copy data between global and shared memory (indeed you can copy data between two arbitrary memory regims on device - as long as you have read respectively write rights)
//blockDim.y and blockDim.z are assumed to be 1
static inline __device__ void copySharedToGlob(int* shared, int* glob, const int &n)
{
  int i, k;
  k = n - blockDim.x;

  for(i = 0; i < k; i += blockDim.x) {
    glob[i + threadIdx.x] = shared[i + threadIdx.x];
  }

  if(threadIdx.x < n - i) {
    glob[i + threadIdx.x] = shared[i + threadIdx.x];
  }

  __syncthreads();
}

static inline __device__ void copySharedToGlob(float* shared, float* glob, const int &n)
{
  int i, k;
  k = n - blockDim.x;

  for(i = 0; i < k; i += blockDim.x) {
    glob[i + threadIdx.x] = shared[i + threadIdx.x];
  }

  if(threadIdx.x < n - i) {
    glob[i + threadIdx.x] = shared[i + threadIdx.x];
  }

  __syncthreads();
}

static inline __device__ void copySharedToGlob(double* shared, double* glob, const int &n)
{
  int i, k;
  k = n - blockDim.x;

  for(i = 0; i < k; i += blockDim.x) {
    glob[i + threadIdx.x] = shared[i + threadIdx.x];
  }

  if(threadIdx.x < n - i) {
    glob[i + threadIdx.x] = shared[i + threadIdx.x];
  }

  __syncthreads();
}

static inline __device__ void copyGlobToShared(int* glob, int* shared, const int &n)
{
  int i, k;
  k = n - blockDim.x;

  for(i = 0; i < k; i += blockDim.x) {
    shared[i + threadIdx.x] = glob[i + threadIdx.x];
  }

  if(threadIdx.x < n - i) {
    shared[i + threadIdx.x] = glob[i + threadIdx.x];
  }

  __syncthreads();
}

static __device__ inline void copyGlobToShared(float* glob, float* shared, const int &n)
{
  int i, k;
  k = n - blockDim.x;

  for(i = 0; i < k; i += blockDim.x) {
    shared[i + threadIdx.x] = glob[i + threadIdx.x];
  }

  if(threadIdx.x < n - i) {
    shared[i + threadIdx.x] = glob[i + threadIdx.x];
  }

  __syncthreads();
}

static __device__ inline void copyGlobToShared(double* glob, double* shared, const int &n)
{
  int i;

  for(i = 0; i < n - blockDim.x; i += blockDim.x) {
    shared[i + threadIdx.x] = glob[i + threadIdx.x];
  }

  if(threadIdx.x < n - i) {
    shared[i + threadIdx.x] = glob[i + threadIdx.x];
  }

  __syncthreads();
}

//copy data between two memory areas on device, 3d BlockDims are allowed
static __device__ inline void copyData(double* source, double* target, const int &n)
{
  int i;
  int offset = threadIdx.x * blockDim.y * blockDim.z + threadIdx.y * blockDim.z + threadIdx.z;

  for(i = 0; i < n - blockDim.x * blockDim.y * blockDim.z; i += blockDim.x * blockDim.y * blockDim.z) {
    target[i + offset] = source[i + offset];
  }

  if(offset < n - i) {
    target[i + offset] = source[i + offset];
  }

  __syncthreads();
}

static __device__ inline void copyData(float* source, float* target, const int &n)
{
  int i;
  int offset = threadIdx.x * blockDim.y * blockDim.z + threadIdx.y * blockDim.z + threadIdx.z;

  for(i = 0; i < n - blockDim.x * blockDim.y * blockDim.z; i += blockDim.x * blockDim.y * blockDim.z) {
    target[i + offset] = source[i + offset];
  }

  if(offset < n - i) {
    target[i + offset] = source[i + offset];
  }

  __syncthreads();
}

static __device__ inline void copyData(int* source, int* target, const int &n)
{
  int i;
  int offset = threadIdx.x * blockDim.y * blockDim.z + threadIdx.y * blockDim.z + threadIdx.z;

  for(i = 0; i < n - blockDim.x * blockDim.y * blockDim.z; i += blockDim.x * blockDim.y * blockDim.z) {
    target[i + offset] = source[i + offset];
  }

  if(offset < n - i) {
    target[i + offset] = source[i + offset];
  }

  __syncthreads();
}

static __device__ inline void copyData(unsigned int* source, unsigned int* target, const int &n)
{
  int i;
  int offset = threadIdx.x * blockDim.y * blockDim.z + threadIdx.y * blockDim.z + threadIdx.z;

  for(i = 0; i < n - blockDim.x * blockDim.y * blockDim.z; i += blockDim.x * blockDim.y * blockDim.z) {
    target[i + offset] = source[i + offset];
  }

  if(offset < n - i) {
    target[i + offset] = source[i + offset];
  }

  __syncthreads();
}

//functions in order to sum over values of one block. P2 means blockdim MUST be a power of 2 otherwise the behaviour is not well defined
//in the end in data[0]=sum_i=0^blockDim.x data[i]
//for reduceBlockP2 and reduceBlock blockDim.y=1 and blockDim.z=1
static __device__ inline void reduceBlockP2(int* data)
{
  __syncthreads();

  for(int i = 2; i <= blockDim.x; i *= 2) {
    if(threadIdx.x < blockDim.x / i)
      data[threadIdx.x] += data[threadIdx.x + blockDim.x / i];

    __syncthreads();
  }
}

static __device__ inline void reduceBlockP2(unsigned int* data)
{
  __syncthreads();

  for(int i = 2; i <= blockDim.x; i *= 2) {
    if(threadIdx.x < blockDim.x / i)
      data[threadIdx.x] += data[threadIdx.x + blockDim.x / i];

    __syncthreads();
  }
}

static __device__ inline void reduceBlockP2(float* data)
{
  __syncthreads();

  for(int i = 2; i <= blockDim.x; i *= 2) {
    if(threadIdx.x < blockDim.x / i)
      data[threadIdx.x] += data[threadIdx.x + blockDim.x / i];

    __syncthreads();
  }
}

static __device__ inline void reduceBlockP2(double* data)
{
  __syncthreads();

  for(int i = 2; i <= blockDim.x; i *= 2) {
    if(threadIdx.x < blockDim.x / i)
      data[threadIdx.x] += data[threadIdx.x + blockDim.x / i];

    __syncthreads();
  }
}

static __device__ inline void reduceBlock(float* data)
{
  __syncthreads();
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

static __device__ inline void reduceBlock(int* data)
{
  __syncthreads();
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

static __device__ inline void reduceBlock(unsigned int* data)
{
  __syncthreads();
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

static __device__ inline void reduceBlock(double* data)
{
  __syncthreads();
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

static __device__ inline void cudaFillBlockData_int(int* data, const int &n, const int &value)
{
  int i;

  for(i = 0; i < n - blockDim.x; i += blockDim.x) {
    data[i + threadIdx.x] = value;
  }

  if(threadIdx.x < n - i) data[i + threadIdx.x] = value;
}

static __device__ inline void cudaFillBlockData_float(float* data, const int &n, const float &value)
{
  int i;

  for(i = 0; i < n - blockDim.x; i += blockDim.x) {
    data[i + threadIdx.x] = value;
  }

  if(threadIdx.x < n - i) data[i + threadIdx.x] = value;
}

static __device__ inline void reduce(float* data, int n) //cautious not sure if working
{
  __syncthreads();
  int p2 = 1;

  while(p2 * 2 < n) p2 *= 2;

  int j = 0;

  while((threadIdx.x + blockDim.x * j) * 2 < n - p2) {
    data[threadIdx.x + blockDim.x * j] += data[(threadIdx.x + blockDim.x * j) + p2];
    j++;
  }

  __syncthreads();

  for(int i = 2; i <= p2; i *= 2) {
    while((threadIdx.x + blockDim.x * j) < p2 / i) {
      data[threadIdx.x + blockDim.x * j] += data[(threadIdx.x + blockDim.x * j) + p2 / i];
      j++;
    }

    __syncthreads();
  }
}

static __device__ inline void reduce(double* data, int n) //cautious not sure if working
{
  __syncthreads();
  int p2 = 1;

  while(p2 * 2 < n) p2 *= 2;

  int j = 0;

  while((threadIdx.x + blockDim.x * j) * 2 < n - p2) {
    data[threadIdx.x + blockDim.x * j] += data[(threadIdx.x + blockDim.x * j) + p2];
    j++;
  }

  __syncthreads();

  for(int i = 2; i <= p2; i *= 2) {
    while((threadIdx.x + blockDim.x * j) < p2 / i) {
      data[threadIdx.x + blockDim.x * j] += data[(threadIdx.x + blockDim.x * j) + p2 / i];
      j++;
    }

    __syncthreads();
  }
}

static __device__ inline void minOfBlock(float* data)
{
  __syncthreads();
  int p2 = 1;

  while(p2 * 2 < blockDim.x) p2 *= 2;

  if(threadIdx.x < blockDim.x - p2)
    data[threadIdx.x] = MIN(data[threadIdx.x + p2], data[threadIdx.x]);

  __syncthreads();

  for(int i = 2; i <= p2; i *= 2) {
    if(threadIdx.x < p2 / i)
      data[threadIdx.x] = MIN(data[threadIdx.x + p2 / i], data[threadIdx.x]);

    __syncthreads();
  }
}

static __device__ inline void maxOfBlock(float* data)
{
  __syncthreads();
  int p2 = 1;

  while(p2 * 2 < blockDim.x) p2 *= 2;

  if(threadIdx.x < blockDim.x - p2)
    data[threadIdx.x] = MAX(data[threadIdx.x + p2], data[threadIdx.x]);

  __syncthreads();

  for(int i = 2; i <= p2; i *= 2) {
    if(threadIdx.x < p2 / i)
      data[threadIdx.x] = MAX(data[threadIdx.x + p2 / i], data[threadIdx.x]);

    __syncthreads();
  }
}

static __device__ inline void minOfBlock(double* data)
{
  __syncthreads();
  int p2 = 1;

  while(p2 * 2 < blockDim.x) p2 *= 2;

  if(threadIdx.x < blockDim.x - p2)
    data[threadIdx.x] = MIN(data[threadIdx.x + p2], data[threadIdx.x]);

  __syncthreads();

  for(int i = 2; i <= p2; i *= 2) {
    if(threadIdx.x < p2 / i)
      data[threadIdx.x] = MIN(data[threadIdx.x + p2 / i], data[threadIdx.x]);

    __syncthreads();
  }
}

static __device__ inline void maxOfBlock(double* data)
{
  __syncthreads();
  int p2 = 1;

  while(p2 * 2 < blockDim.x) p2 *= 2;

  if(threadIdx.x < blockDim.x - p2)
    data[threadIdx.x] = MAX(data[threadIdx.x + p2], data[threadIdx.x]);

  __syncthreads();

  for(int i = 2; i <= p2; i *= 2) {
    if(threadIdx.x < p2 / i)
      data[threadIdx.x] = MAX(data[threadIdx.x + p2 / i], data[threadIdx.x]);

    __syncthreads();
  }
}


static __device__ inline void minOfData(double* data, int n) //cautious not sure if working
{
  __syncthreads();
  int p2 = 1;

  while(p2 * 2 < n) p2 *= 2;

  int j = 0;

  while((threadIdx.x + blockDim.x * j) < n - p2) {
    data[threadIdx.x + blockDim.x * j] = MIN(data[threadIdx.x + blockDim.x * j], data[(threadIdx.x + blockDim.x * j) + p2]);
    j++;
  }

  __syncthreads();

  for(int i = 2; i <= p2; i *= 2) {
    while((threadIdx.x + blockDim.x * j) < p2 / i) {
      data[threadIdx.x + blockDim.x * j] = MIN(data[threadIdx.x + blockDim.x * j], data[(threadIdx.x + blockDim.x * j) + p2 / i]);
      j++;
    }

    __syncthreads();
  }
}

static __device__ inline void maxOfData(double* data, int n) //cautious not sure if working
{
  __syncthreads();
  int p2 = 1;

  while(p2 * 2 < n) p2 *= 2;

  int j = 0;

  while((threadIdx.x + blockDim.x * j) < n - p2) {
    data[threadIdx.x + blockDim.x * j] = MAX(data[threadIdx.x + blockDim.x * j], data[(threadIdx.x + blockDim.x * j) + p2]);
    j++;
  }

  __syncthreads();

  for(int i = 2; i <= p2; i *= 2) {
    while((threadIdx.x + blockDim.x * j) < p2 / i) {
      data[threadIdx.x + blockDim.x * j] = MAX(data[threadIdx.x + blockDim.x * j], data[(threadIdx.x + blockDim.x * j) + p2 / i]);
      j++;
    }

    __syncthreads();
  }
}

static __device__ inline void minOfData(float* data, int n) //cautious not sure if working
{
  __syncthreads();
  int p2 = 1;

  while(p2 * 2 < n) p2 *= 2;

  int j = 0;

  while((threadIdx.x + blockDim.x * j) < n - p2) {
    data[threadIdx.x + blockDim.x * j] = MIN(data[threadIdx.x + blockDim.x * j], data[(threadIdx.x + blockDim.x * j) + p2]);
    j++;
  }

  __syncthreads();

  for(int i = 2; i <= p2; i *= 2) {
    while((threadIdx.x + blockDim.x * j) < p2 / i) {
      data[threadIdx.x + blockDim.x * j] = MIN(data[threadIdx.x + blockDim.x * j], data[(threadIdx.x + blockDim.x * j) + p2 / i]);
      j++;
    }

    __syncthreads();
  }
}

static __device__ inline void maxOfData(float* data, int n) //cautious not sure if working
{
  __syncthreads();
  int p2 = 1;

  while(p2 * 2 < n) p2 *= 2;

  int j = 0;

  while((threadIdx.x + blockDim.x * j) < n - p2) {
    data[threadIdx.x + blockDim.x * j] = MAX(data[threadIdx.x + blockDim.x * j], data[(threadIdx.x + blockDim.x * j) + p2]);
    j++;
  }

  __syncthreads();

  for(int i = 2; i <= p2; i *= 2) {
    while((threadIdx.x + blockDim.x * j) < p2 / i) {
      data[threadIdx.x + blockDim.x * j] = MAX(data[threadIdx.x + blockDim.x * j], data[(threadIdx.x + blockDim.x * j) + p2 / i]);
      j++;
    }

    __syncthreads();
  }
}

#if X_PRECISION == 2
static __device__ inline double tex1Dfetch_double(texture<int2, 1> t, int i)
{
  int2 v = tex1Dfetch(t, i);
  return __hiloint2double(v.y, v.x);
}

static __device__ inline X_CFLOAT4 tex1Dfetch_double(texture<int4, 1> t, int i)
{
  int4 v = tex1Dfetch(t, 2 * i);
  int4 u = tex1Dfetch(t, 2 * i + 1);
  X_CFLOAT4 w;

  w.x = __hiloint2double(v.y, v.x);
  w.y = __hiloint2double(v.w, v.z);
  w.z = __hiloint2double(u.y, u.x);
  w.w = __hiloint2double(u.w, u.z);
  return w;
}
#endif

inline void BindXTypeTexture(cuda_shared_data* sdata)
{
#ifdef CUDA_USE_TEXTURE
  _x_type_tex.normalized = false;                      // access with normalized texture coordinates
  _x_type_tex.filterMode = cudaFilterModePoint;        // Point mode, so no
  _x_type_tex.addressMode[0] = cudaAddressModeWrap;    // wrap texture coordinates
  const textureReference* x_type_texture_ptr = &MY_AP(x_type_tex);

#if X_PRECISION == 1
  cudaChannelFormatDesc channelDescXType = cudaCreateChannelDesc<float4>();
  cudaBindTexture(0, x_type_texture_ptr, sdata->atom.x_type.dev_data, &channelDescXType, sdata->atom.nmax * sizeof(X_CFLOAT4));
#else
  cudaChannelFormatDesc channelDescXType = cudaCreateChannelDesc<int4>();
  cudaBindTexture(0, x_type_texture_ptr, sdata->atom.x_type.dev_data, &channelDescXType, sdata->atom.nmax * 2 * sizeof(int4));
#endif
#endif
}

static __device__ inline X_CFLOAT4 fetchXType(int i)
{
#ifdef CUDA_USE_TEXTURE
#if X_PRECISION == 1
  return tex1Dfetch(_x_type_tex, i);
#else
  return tex1Dfetch_double(_x_type_tex, i);
#endif
#else
  return _x_type[i];
#endif
}

#if V_PRECISION == 2
static __device__ inline double tex1Dfetch_double_v(texture<int2, 1> t, int i)
{
  int2 v = tex1Dfetch(t, i);
  return __hiloint2double(v.y, v.x);
}

static __device__ inline V_CFLOAT4 tex1Dfetch_double_v(texture<int4, 1> t, int i)
{
  int4 v = tex1Dfetch(t, 2 * i);
  int4 u = tex1Dfetch(t, 2 * i + 1);
  V_CFLOAT4 w;

  w.x = __hiloint2double(v.y, v.x);
  w.y = __hiloint2double(v.w, v.z);
  w.z = __hiloint2double(u.y, u.x);
  w.w = __hiloint2double(u.w, u.z);
  return w;
}
#endif

inline void BindVRadiusTexture(cuda_shared_data* sdata)
{
#ifdef CUDA_USE_TEXTURE
  _v_radius_tex.normalized = false;                      // access with normalized texture coordinates
  _v_radius_tex.filterMode = cudaFilterModePoint;        // Point mode, so no
  _v_radius_tex.addressMode[0] = cudaAddressModeWrap;    // wrap texture coordinates
  const textureReference* v_radius_texture_ptr = &MY_AP(v_radius_tex);

#if V_PRECISION == 1
  cudaChannelFormatDesc channelDescVRadius = cudaCreateChannelDesc<float4>();
  cudaBindTexture(0, v_radius_texture_ptr, sdata->atom.v_radius.dev_data, &channelDescVRadius, sdata->atom.nmax * sizeof(X_CFLOAT4));
#else
  cudaChannelFormatDesc channelDescVRadius = cudaCreateChannelDesc<int4>();
  cudaBindTexture(0, v_radius_texture_ptr, sdata->atom.v_radius.dev_data, &channelDescVRadius, sdata->atom.nmax * 2 * sizeof(int4));
#endif
#endif
}

static __device__ inline V_CFLOAT4 fetchVRadius(int i)
{
#ifdef CUDA_USE_TEXTURE
#if V_PRECISION == 1
  return tex1Dfetch(_v_radius_tex, i);
#else
  return tex1Dfetch_double_v(_v_radius_tex, i);
#endif
#else
  return _v_radius[i];
#endif
}

inline void BindOmegaRmassTexture(cuda_shared_data* sdata)
{
#ifdef CUDA_USE_TEXTURE
  _omega_rmass_tex.normalized = false;                      // access with normalized texture coordinates
  _omega_rmass_tex.filterMode = cudaFilterModePoint;        // Point mode, so no
  _omega_rmass_tex.addressMode[0] = cudaAddressModeWrap;    // wrap texture coordinates
  const textureReference* omega_rmass_texture_ptr = &MY_AP(omega_rmass_tex);

#if V_PRECISION == 1
  cudaChannelFormatDesc channelDescOmegaRmass = cudaCreateChannelDesc<float4>();
  cudaBindTexture(0, omega_rmass_texture_ptr, sdata->atom.omega_rmass.dev_data, &channelDescOmegaRmass, sdata->atom.nmax * sizeof(X_CFLOAT4));
#else
  cudaChannelFormatDesc channelDescOmegaRmass = cudaCreateChannelDesc<int4>();
  cudaBindTexture(0, omega_rmass_texture_ptr, sdata->atom.omega_rmass.dev_data, &channelDescOmegaRmass, sdata->atom.nmax * 2 * sizeof(int4));
#endif
#endif
}

static __device__ inline V_CFLOAT4 fetchOmegaRmass(int i)
{
#ifdef CUDA_USE_TEXTURE
#if V_PRECISION == 1
  return tex1Dfetch(_omega_rmass_tex, i);
#else
  return tex1Dfetch_double_v(_omega_rmass_tex, i);
#endif
#else
  return _omega_rmass[i];
#endif
}

#if F_PRECISION == 2
static __device__ inline double tex1Dfetch_double_f(texture<int2, 1> t, int i)
{
  int2 v = tex1Dfetch(t, i);
  return __hiloint2double(v.y, v.x);
}

static __device__ inline F_CFLOAT4 tex1Dfetch_double_f(texture<int4, 1> t, int i)
{
  int4 v = tex1Dfetch(t, 2 * i);
  int4 u = tex1Dfetch(t, 2 * i + 1);
  F_CFLOAT4 w;

  w.x = __hiloint2double(v.y, v.x);
  w.y = __hiloint2double(v.w, v.z);
  w.z = __hiloint2double(u.y, u.x);
  w.w = __hiloint2double(u.w, u.z);
  return w;
}
#endif

inline void BindQTexture(cuda_shared_data* sdata)
{
#ifdef CUDA_USE_TEXTURE
  _q_tex.normalized = false;                      // access with normalized texture coordinates
  _q_tex.filterMode = cudaFilterModePoint;        // Point mode, so no
  _q_tex.addressMode[0] = cudaAddressModeWrap;    // wrap texture coordinates
  const textureReference* q_texture_ptr = &MY_AP(q_tex);

#if F_PRECISION == 1
  cudaChannelFormatDesc channelDescQ = cudaCreateChannelDesc<float>();
  cudaBindTexture(0, q_texture_ptr, sdata->atom.q.dev_data, &channelDescQ, sdata->atom.nmax * sizeof(F_CFLOAT));
#else
  cudaChannelFormatDesc channelDescQ = cudaCreateChannelDesc<int2>();
  cudaBindTexture(0, q_texture_ptr, sdata->atom.q.dev_data, &channelDescQ, sdata->atom.nmax * sizeof(int2));
#endif
#endif
}

static __device__ inline F_CFLOAT fetchQ(int i)
{
#ifdef CUDA_USE_TEXTURE
#if F_PRECISION == 1
  return tex1Dfetch(_q_tex, i);
#else
  return tex1Dfetch_double_f(_q_tex, i);
#endif
#else
  return _q[i];
#endif
}

#endif

/*

inline void BindPairCoeffTypeTexture(cuda_shared_data* sdata,coeff_tex)
{
	#ifdef CUDA_USE_TEXTURE
		_coeff_tex.normalized = false;                      // access with normalized texture coordinates
		_coeff_tex.filterMode = cudaFilterModePoint;        // Point mode, so no
		_coeff_tex.addressMode[0] = cudaAddressModeWrap;    // wrap texture coordinates
		const textureReference* coeff_texture_ptr;
		cudaGetTextureReference(&coeff_texture_ptr, &MY_AP(coeff_tex));

		#if F_PRECISION == 1
		cudaChannelFormatDesc channelDescXType = cudaCreateChannelDesc<float4>();
		cudaBindTexture(0,x_type_texture_ptr, sdata->atom.x_type.dev_data, &channelDescXType, sdata->atom.nmax*sizeof(X_CFLOAT4));
		#else
		cudaChannelFormatDesc channelDescXType = cudaCreateChannelDesc<int4>();
		cudaBindTexture(0,x_type_texture_ptr, sdata->atom.x_type.dev_data, &channelDescXType, sdata->atom.nmax*2*sizeof(int4));
		#endif
	#endif
}

static __device__ inline X_CFLOAT4 fetchXType(int i)
{
		#ifdef CUDA_USE_TEXTURE
		  #if X_PRECISION == 1
		     return tex1Dfetch(_x_type_tex,i);
		  #else
		     return tex1Dfetch_double(_x_type_tex,i);
		  #endif
		#else
		  return _x_type[i];
		#endif
}
*/
#define SBBITS 30

static inline __device__ int sbmask(int j)
{
  return j >> SBBITS & 3;
}

static inline __device__ void minimum_image(X_CFLOAT4 &delta)
{
  if(_triclinic == 0) {
    if(_periodicity[0]) {
      delta.x += delta.x < -X_F(0.5) * _prd[0] ? _prd[0] :
                 (delta.x >  X_F(0.5) * _prd[0] ? -_prd[0] : X_F(0.0));
    }

    if(_periodicity[1]) {
      delta.y += delta.y < -X_F(0.5) * _prd[1] ? _prd[1] :
                 (delta.y >  X_F(0.5) * _prd[1] ? -_prd[1] : X_F(0.0));
    }

    if(_periodicity[2]) {
      delta.z += delta.z < -X_F(0.5) * _prd[2] ? _prd[2] :
                 (delta.z >  X_F(0.5) * _prd[2] ? -_prd[2] : X_F(0.0));
    }

  } else {
    if(_periodicity[1]) {
      delta.z += delta.z < -X_F(0.5) * _prd[2] ? _prd[2] :
                 (delta.z >  X_F(0.5) * _prd[2] ? -_prd[2] : X_F(0.0));
      delta.y += delta.z < -X_F(0.5) * _prd[2] ? _h[3] :
                 (delta.z >  X_F(0.5) * _prd[2] ? -_h[3] : X_F(0.0));
      delta.x += delta.z < -X_F(0.5) * _prd[2] ? _h[4] :
                 (delta.z >  X_F(0.5) * _prd[2] ? -_h[4] : X_F(0.0));

    }

    if(_periodicity[1]) {
      delta.y += delta.y < -X_F(0.5) * _prd[1] ? _prd[1] :
                 (delta.y >  X_F(0.5) * _prd[1] ? -_prd[1] : X_F(0.0));
      delta.x += delta.y < -X_F(0.5) * _prd[1] ? _h[5] :
                 (delta.y >  X_F(0.5) * _prd[1] ? -_h[5] : X_F(0.0));

    }

    if(_periodicity[0]) {
      delta.x += delta.x < -X_F(0.5) * _prd[0] ? _prd[0] :
                 (delta.x >  X_F(0.5) * _prd[0] ? -_prd[0] : X_F(0.0));
    }
  }
}

static inline __device__ void closest_image(X_CFLOAT4 &x1, X_CFLOAT4 &x2, X_CFLOAT4 &ci)
{
  ci.x = x2.x - x1.x;
  ci.y = x2.y - x1.y;
  ci.z = x2.z - x1.z;
  minimum_image(ci);
  ci.x += x1.x;
  ci.y += x1.y;
  ci.z += x1.z;
}
