enum copy_mode {x, xx, xy, yx, xyz, xzy}; // yxz, yzx, zxy, zyx not yet implemented since they were not needed yet

#include "cuda_data_cu.h"
#include "cuda_wrapper_cu.h"
#include "cuda_data_kernel.cu"
#include <cstdio>

void CudaData_Upload_DoubleFloat(void* host_data, void* dev_data, unsigned* n, copy_mode mode, void* buffer)
{
  int size = n[0];

  if(n[1] > 0) size *= n[1];

  if(n[2] > 0) size *= n[2];

  dim3 threads;
  threads.x = 1;
  threads.y = 1;
  threads.z = 1;
  dim3 grid;
  grid.x = 1;
  grid.y = 1;
  grid.z = 1;

  if(size <= 128 * 30)
    threads.x = 32;
  else if(size <= 256 * 30)
    threads.x = 64;
  else if(size <= 512 * 30)
    threads.x = 128;
  else
    threads.x = 256;

  grid.x = ((size - 1) + threads.x) / threads.x;

  if(grid.x > 32000)
    grid.x = 32000;

  while(grid.x * grid.y * threads.x < size) grid.y++;

  float debugdata[size];
  //int* cu_debug=(int*) CudaWrapper_AllocCudaData(size*sizeof(FLOAT));
  size *= sizeof(double);
  printf("size: %i (%i %i %i) (%i %i %i) %p\n", size, grid.x, grid.y, threads.x, n[0], n[1], n[2], buffer);
  CudaWrapper_UploadCudaData(host_data, buffer, size);
  CudaData_Upload_Kernel_DoubleFloat <<< grid, threads>>>((double*)buffer, (float*)dev_data, n[0], n[1], n[2], mode);
  cudaThreadSynchronize();
  CudaWrapper_DownloadCudaData(debugdata, dev_data, size / 2);
  double sum = 0;
  printf("debugdata: ");

  for(int i = 0; i < size / sizeof(double); i++) sum += (debugdata[i] - ((double*) host_data)[i]) * (debugdata[i] - ((double*) host_data)[i]);

  printf("%lf \n", sum);

}

void CudaData_Upload_DoubleDouble(void* host_data, void* dev_data, unsigned* n, copy_mode mode, void* buffer)
{
  int size = n[0];

  if(n[1] > 0) size *= n[1];

  if(n[2] > 0) size *= n[2];

  dim3 threads;
  threads.x = 1;
  threads.y = 1;
  threads.z = 1;
  dim3 grid;
  grid.x = 1;
  grid.y = 1;
  grid.z = 1;

  if(size <= 128 * 30)
    threads.x = 32;
  else if(size <= 256 * 30)
    threads.x = 64;
  else if(size <= 512 * 30)
    threads.x = 128;
  else
    threads.x = 256;

  grid.x = ((size - 1) + threads.x) / threads.x;

  if(grid.x > 32000)
    grid.x = 32000;

  while(grid.x * grid.y * threads.x < size) grid.y++;

  size *= sizeof(double);

  CudaWrapper_UploadCudaData(host_data, buffer, size);
  CudaData_Upload_Kernel_DoubleDouble <<< grid, threads>>>((double*)buffer, (double*)dev_data, n[0], n[1], n[2], mode);
  cudaThreadSynchronize();
}

void CudaData_Upload_FloatDouble(void* host_data, void* dev_data, unsigned* n, copy_mode mode, void* buffer)
{
  int size = n[0];

  if(n[1] > 0) size *= n[1];

  if(n[2] > 0) size *= n[2];

  dim3 threads;
  threads.x = 1;
  threads.y = 1;
  threads.z = 1;
  dim3 grid;
  grid.x = 1;
  grid.y = 1;
  grid.z = 1;

  if(size <= 128 * 30)
    threads.x = 32;
  else if(size <= 256 * 30)
    threads.x = 64;
  else if(size <= 512 * 30)
    threads.x = 128;
  else
    threads.x = 256;

  grid.x = ((size - 1) + threads.x) / threads.x;

  if(grid.x > 32000)
    grid.x = 32000;

  while(grid.x * grid.y * threads.x < size) grid.y++;

  size *= sizeof(float);

  CudaWrapper_UploadCudaData(host_data, buffer, size);
  CudaData_Upload_Kernel_FloatDouble <<< grid, threads>>>((float*)buffer, (double*)dev_data, n[0], n[1], n[2], mode);
  cudaThreadSynchronize();
}

void CudaData_Upload_FloatFloat(void* host_data, void* dev_data, unsigned* n, copy_mode mode, void* buffer)
{
  int size = n[0];

  if(n[1] > 0) size *= n[1];

  if(n[2] > 0) size *= n[2];

  dim3 threads;
  threads.x = 1;
  threads.y = 1;
  threads.z = 1;
  dim3 grid;
  grid.x = 1;
  grid.y = 1;
  grid.z = 1;

  if(size <= 128 * 30)
    threads.x = 32;
  else if(size <= 256 * 30)
    threads.x = 64;
  else if(size <= 512 * 30)
    threads.x = 128;
  else
    threads.x = 256;

  grid.x = ((size - 1) + threads.x) / threads.x;

  if(grid.x > 32000)
    grid.x = 32000;

  while(grid.x * grid.y * threads.x < size) grid.y++;

  size *= sizeof(float);

  CudaWrapper_UploadCudaData(host_data, buffer, size);
  CudaData_Upload_Kernel_FloatFloat <<< grid, threads>>>((float*)buffer, (float*)dev_data, n[0], n[1], n[2], mode);
  cudaThreadSynchronize();
}

void CudaData_Upload_IntInt(void* host_data, void* dev_data, unsigned* n, copy_mode mode, void* buffer)
{
  int size = n[0];

  if(n[1] > 0) size *= n[1];

  if(n[2] > 0) size *= n[2];

  dim3 threads;
  threads.x = 1;
  threads.y = 1;
  threads.z = 1;
  dim3 grid;
  grid.x = 1;
  grid.y = 1;
  grid.z = 1;

  if(size <= 128 * 30)
    threads.x = 32;
  else if(size <= 256 * 30)
    threads.x = 64;
  else if(size <= 512 * 30)
    threads.x = 128;
  else
    threads.x = 256;

  grid.x = ((size - 1) + threads.x) / threads.x;

  if(grid.x > 32000)
    grid.x = 32000;

  while(grid.x * grid.y * threads.x < size) grid.y++;

  size *= sizeof(int);

  CudaWrapper_UploadCudaData(host_data, buffer, size);
  CudaData_Upload_Kernel_IntInt <<< grid, threads>>>((int*)buffer, (int*)dev_data, n[0], n[1], n[2], mode);
  cudaThreadSynchronize();
}

void CudaData_Download(void* host_data, void* dev_data, int host_size, int dev_size, unsigned* n, copy_mode mode, void* buffer)
{
}
