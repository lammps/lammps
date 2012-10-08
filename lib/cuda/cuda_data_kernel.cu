__global__ void CudaData_Upload_Kernel_DoubleFloat(double* buffer, float* dev_data,
    unsigned nx, unsigned ny, unsigned nz, copy_mode mode)
{
  if(mode == x) mode = xx;

  unsigned length = nx;

  if(ny > 0) length *= ny;

  if(nz > 0) length *= nz;

  unsigned i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x, j, k, l;


  if(i >= length) return;

  switch(mode) {
    case xx: {
      dev_data[i] = buffer[i];
    }

    case xy: {
      dev_data[i] = buffer[i];
    }

    case yx: {
      j = i / ny;
      k = i % ny;
      dev_data[k * nx + j] = buffer[j * ny + k];
    }

    case xyz: {
      dev_data[i] = buffer[i];
    }

    case xzy: {
      j = i / (ny * nz);
      k = (i % (ny * nz)) / nz;
      l = i % nz;
      dev_data[j * ny * nz + l * ny + k] = buffer[j * ny * nz + k * nz + l];
    }
  }
}

__global__ void CudaData_Upload_Kernel_DoubleDouble(double* buffer, double* dev_data,
    unsigned nx, unsigned ny, unsigned nz, copy_mode mode)
{
  if(mode == x) mode = xx;

  unsigned length = nx;

  if(ny > 0) length *= ny;

  if(nz > 0) length *= nz;

  unsigned i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x, j, k, l;

  if(i >= length) return;

  switch(mode) {
    case xx:
      dev_data[i] = buffer[i];

    case xy:
      dev_data[i] = buffer[i];

    case yx:
      j = i / ny;
      k = i % ny;
      dev_data[k * nx + j] = buffer[j * ny + k];

    case xyz:
      dev_data[i] = buffer[i];

    case xzy:
      j = i / (ny * nz);
      k = (i % (ny * nz)) / nz;
      l = i % nz;
      dev_data[j * ny * nz + l * ny + k] = buffer[j * ny * nz + k * nz + l];
  }
}

__global__ void CudaData_Upload_Kernel_FloatDouble(float* buffer, double* dev_data,
    unsigned nx, unsigned ny, unsigned nz, copy_mode mode)
{
  if(mode == x) mode = xx;

  unsigned length = nx;

  if(ny > 0) length *= ny;

  if(nz > 0) length *= nz;

  unsigned i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x, j, k, l;

  if(i >= length) return;

  switch(mode) {
    case xx:
      dev_data[i] = buffer[i];

    case xy:
      dev_data[i] = buffer[i];

    case yx:
      j = i / ny;
      k = i % ny;
      dev_data[k * nx + j] = buffer[j * ny + k];

    case xyz:
      dev_data[i] = buffer[i];

    case xzy:
      j = i / (ny * nz);
      k = (i % (ny * nz)) / nz;
      l = i % nz;
      dev_data[j * ny * nz + l * ny + k] = buffer[j * ny * nz + k * nz + l];
  }
}

__global__ void CudaData_Upload_Kernel_FloatFloat(float* buffer, float* dev_data,
    unsigned nx, unsigned ny, unsigned nz, copy_mode mode)
{
  if(mode == x) mode = xx;

  unsigned length = nx;

  if(ny > 0) length *= ny;

  if(nz > 0) length *= nz;

  unsigned i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x, j, k, l;

  if(i >= length) return;

  switch(mode) {
    case xx:
      dev_data[i] = buffer[i];

    case xy:
      dev_data[i] = buffer[i];

    case yx:
      j = i / ny;
      k = i % ny;
      dev_data[k * nx + j] = buffer[j * ny + k];

    case xyz:
      dev_data[i] = buffer[i];

    case xzy:
      j = i / (ny * nz);
      k = (i % (ny * nz)) / nz;
      l = i % nz;
      dev_data[j * ny * nz + l * ny + k] = buffer[j * ny * nz + k * nz + l];
  }
}

__global__ void CudaData_Upload_Kernel_IntInt(int* buffer, int* dev_data,
    unsigned nx, unsigned ny, unsigned nz, copy_mode mode)
{
  if(mode == x) mode = xx;

  unsigned length = nx;

  if(ny > 0) length *= ny;

  if(nz > 0) length *= nz;

  unsigned i = (blockIdx.x * gridDim.y + blockIdx.y) * blockDim.x + threadIdx.x, j, k, l;

  if(i >= length) return;

  switch(mode) {
    case xx:
      dev_data[i] = buffer[i];

    case xy:
      dev_data[i] = buffer[i];

    case yx:
      j = i / ny;
      k = i % ny;
      dev_data[k * nx + j] = buffer[j * ny + k];

    case xyz:
      dev_data[i] = buffer[i];

    case xzy:
      j = i / (ny * nz);
      k = (i % (ny * nz)) / nz;
      l = i % nz;
      dev_data[j * ny * nz + l * ny + k] = buffer[j * ny * nz + k * nz + l];
  }
}
