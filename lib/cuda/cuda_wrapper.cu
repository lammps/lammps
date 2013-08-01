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

#include <stdio.h>
#include "cuda_shared.h"
#include "cuda_common.h"
#include "cuda_wrapper_cu.h"
#include "cuda_wrapper_kernel.cu"

static int CudaWrapper_total_gpu_mem = 0;
static double CudaWrapper_total_upload_time = 0;
static double CudaWrapper_total_download_time = 0;
static double CudaWrapper_cpubuffer_upload_time = 0;
static double CudaWrapper_cpubuffer_download_time = 0;
static cudaStream_t* streams;
static int nstreams = 0;

void CudaWrapper_Init(int argc, char** argv, int me, int ppn, int* devicelist)
{
  MYDBG(printf("# CUDA: debug mode on\n");)

#if __DEVICE_EMULATION__

  printf("# CUDA: emulation mode on\n");

#else

  // modified from cutil.h
  static int deviceCount = 0;
  static bool sharedmode = false;

  if(deviceCount && !sharedmode) return;

  if(deviceCount && sharedmode) cudaThreadExit();

  CUDA_SAFE_CALL_NO_SYNC(cudaGetDeviceCount(&deviceCount));

  if(deviceCount == 0) {
    fprintf(stderr, "cutil error: no devices supporting CUDA.\n");
    exit(EXIT_FAILURE);
  }

  MYDBG(printf("# CUDA There are %i devices supporting CUDA in this system.\n", deviceCount);)

  cudaDeviceProp deviceProp[deviceCount];

  for(int i = 0; i < deviceCount; i++)
    CUDA_SAFE_CALL_NO_SYNC(cudaGetDeviceProperties(&(deviceProp[i]), i));


  int dev_list[deviceCount];

  for(int i = 0; i < deviceCount; i++) dev_list[i] = i;

  for(int i = 0; i < deviceCount; i++) {
    for(int j = 0; j < deviceCount - 1 - i; j++)
      if(deviceProp[dev_list[j]].multiProcessorCount < deviceProp[dev_list[j + 1]].multiProcessorCount) {
        int k = dev_list[j];
        dev_list[j] = dev_list[j + 1];
        dev_list[j + 1] = k;
      }
  }

  for(int i = 0; i < deviceCount; i++) {
    if((deviceProp[dev_list[i]].computeMode == 0)) sharedmode = true;

    cudaSetDevice(i);
    cudaSetDeviceFlags(cudaDeviceMapHost);
  }

  if(sharedmode) {
    if(ppn && (me % ppn + 1) > deviceCount) {
      printf("Asking for more GPUs per node when there are. Reduce gpu/node setting.\n");
      exit(0);
    }

    int devicea = me % ppn;

    if(devicelist) devicea = devicelist[devicea];
    else
      devicea = dev_list[devicea];

    if(devicea >= deviceCount)  {
      printf("Asking for non existent GPU %i. Found only %i GPUs.\n", devicea, deviceCount);
      exit(0);
    }

    MYDBG(
      printf(" # CUDA  myid: %i take device: %i\n", me, devicea);
    )
    CUDA_SAFE_CALL(cudaSetDevice(devicea));
  } else {
    CUDA_SAFE_CALL(cudaSetValidDevices(dev_list, deviceCount));
  }

  cudaThreadSynchronize();

  int dev;
  CUDA_SAFE_CALL(cudaGetDevice(&dev));

  if(deviceProp[dev].major < 1) {
    fprintf(stderr, "CUDA error: device does not support CUDA.\n");
    exit(EXIT_FAILURE);
  } else if((deviceProp[dev].major == 1) && (deviceProp[dev].minor != 3)) {
    fprintf(stderr, "CUDA error: You need a device with compute capability 1.3 or higher (Device %i is a %s with CC %i.%i)\n", dev, deviceProp[dev].name, deviceProp[dev].major, deviceProp[dev].minor);
    exit(EXIT_FAILURE);
  }

  if((deviceProp[dev].major == 2) && (CUDA_ARCH < 20)) {
    fprintf(stderr, "CUDA warning: You are using a compute %i.%i or higher GPU while LAMMPScuda has been compiled for architecture 1.3\n", deviceProp[dev].major, deviceProp[dev].minor);
  }

  if((deviceProp[dev].major == 1) && (CUDA_ARCH >= 20)) {
    fprintf(stderr, "CUDA error: You are using a compute 1.3 GPU while LAMMPScuda has been compiled for architecture %i\n", CUDA_ARCH);
    exit(EXIT_FAILURE);
  }


  fprintf(stderr, "# Using device %d: %s\n", dev, deviceProp[dev].name);
  MYDBG(fprintf(stderr, "# Using device %d: %s\n", dev, deviceProp[dev].name);)

  MYDBG
  (
    printf("name = %s\n", deviceProp[dev].name);
    printf("totalGlobalMem = %u\n", deviceProp[dev].totalGlobalMem);
    printf("sharedMemPerBlock = %i\n", deviceProp[dev].sharedMemPerBlock);
    printf("regsPerBlock = %i\n", deviceProp[dev].regsPerBlock);
    printf("warpSize = %i\n", deviceProp[dev].warpSize);
    printf("memPitch = %i\n", deviceProp[dev].memPitch);
    printf("maxThreadsPerBlock = %i\n", deviceProp[dev].maxThreadsPerBlock);
    printf("maxThreadsDim = [%i, %i, %i]\n", deviceProp[dev].maxThreadsDim[0], deviceProp[dev].maxThreadsDim[1], deviceProp[dev].maxThreadsDim[2]);
    printf("maxGridSize = [%i, %i, %i]\n", deviceProp[dev].maxGridSize[0], deviceProp[dev].maxGridSize[1], deviceProp[dev].maxGridSize[2]);
    printf("totalConstMem = %i\n", deviceProp[dev].totalConstMem);
    printf("major . minor = %i . %i\n", deviceProp[dev].major, deviceProp[dev].minor);
    printf("clockRate = %i\n", deviceProp[dev].clockRate);
    printf("textureAlignment = %i\n", deviceProp[dev].textureAlignment);
    printf("deviceOverlap = %i\n", deviceProp[dev].deviceOverlap);
    printf("multiProcessorCount = %i\n", deviceProp[dev].multiProcessorCount);
    printf("computeMode = %i\n", deviceProp[dev].computeMode);
  )

#endif
}

void* CudaWrapper_AllocCudaData(unsigned nbytes)
{
  void* dev_data;
  CUDA_SAFE_CALL(cudaMalloc((void**)&dev_data, nbytes));
  MYDBG(printf("# CUDA: allocated %u bytes on device at dev%p\n", nbytes, dev_data);)
  CudaWrapper_total_gpu_mem += nbytes;
  return dev_data;
}

void CudaWrapper_UploadCudaData(void* host_data, void* dev_data, unsigned nbytes)
{
  MYDBG(printf("# CUDA: uploading %u bytes to device at dev%p from %p\n", nbytes, dev_data, host_data);)
  cudaThreadSynchronize();
  my_times time1, time2;
  my_gettime(CLOCK_REALTIME, &time1);
  CUDA_SAFE_CALL(cudaMemcpy(dev_data, host_data, nbytes, cudaMemcpyHostToDevice));
  my_gettime(CLOCK_REALTIME, &time2);
  CudaWrapper_total_upload_time +=
    time2.tv_sec - time1.tv_sec + 1.0 * (time2.tv_nsec - time1.tv_nsec) / 1000000000;
}

void CudaWrapper_UploadCudaDataAsync(void* host_data, void* dev_data, unsigned nbytes, int stream)
{
  MYDBG(printf("# CUDA: downloading %u bytes from device at dev%p\n", nbytes, dev_data);)
  cudaMemcpyAsync(dev_data, host_data, nbytes, cudaMemcpyHostToDevice, streams[stream]);
}

void CudaWrapper_DownloadCudaData(void* host_data, void* dev_data, unsigned nbytes)
{
  MYDBG(printf("# CUDA: downloading %u bytes from device at dev%p\n", nbytes, dev_data);)
  cudaThreadSynchronize();
  my_times time1, time2;
  my_gettime(CLOCK_REALTIME, &time1);
  CUDA_SAFE_CALL(cudaMemcpy(host_data, dev_data, nbytes, cudaMemcpyDeviceToHost));
  my_gettime(CLOCK_REALTIME, &time2);
  CudaWrapper_total_download_time +=
    time2.tv_sec - time1.tv_sec + 1.0 * (time2.tv_nsec - time1.tv_nsec) / 1000000000;
}

void CudaWrapper_DownloadCudaDataAsync(void* host_data, void* dev_data, unsigned nbytes, int stream)
{
  MYDBG(printf("# CUDA: downloading %u bytes from device at dev%p\n", nbytes, dev_data);)
  cudaMemcpyAsync(host_data, dev_data, nbytes, cudaMemcpyDeviceToHost, streams[stream]);
}

void CudaWrapper_FreeCudaData(void* dev_data, unsigned nbytes)
{
  MYDBG(printf("# CUDA: freeing memory at dev%p with %i bytes (last adress: %p)\n", dev_data, nbytes, (char*)dev_data + nbytes);)
  CUDA_SAFE_CALL(cudaFree(dev_data));
  CudaWrapper_total_gpu_mem -= nbytes;
}

void CudaWrapper_Memset(void* dev_data, int value, unsigned nbytes)
{
  MYDBG(printf("# CUDA: setting %u bytes to %i at dev%p\n", nbytes, value, dev_data);)
  CUDA_SAFE_CALL(cudaMemset(dev_data, value, nbytes));
}

void CudaWrapper_CopyData(void* dev_dest, void* dev_source, unsigned nbytes)
{
  MYDBG(printf("# CUDA: copy %u bytes from dev%p to dev%p\n", nbytes, dev_source, dev_dest);)
  CUDA_SAFE_CALL(cudaMemcpy(dev_dest, dev_source, nbytes, cudaMemcpyDeviceToDevice));
}

void* CudaWrapper_AllocPinnedHostData(unsigned nbytes, bool mapped, bool writeCombined)
{
  void* host_data;
  int flags = 0;

  if(mapped) flags = flags | cudaHostAllocMapped;

  if(writeCombined) flags = flags | cudaHostAllocWriteCombined;

  CUDA_SAFE_CALL(cudaHostAlloc((void**)&host_data, nbytes, flags));
  //	CUDA_SAFE_CALL( cudaMallocHost((void**)&host_data, nbytes) );
  MYDBG(printf("# CUDA: allocated %u bytes pinned memory on host at %p\n", nbytes, host_data);)
  return host_data;
}

void CudaWrapper_FreePinnedHostData(void* host_data)
{
  MYDBG(printf("# CUDA: freeing pinned host memory at %p \n", host_data);)

  if(host_data)
    CUDA_SAFE_CALL(cudaFreeHost(host_data));
}

void cuda_check_error(char* comment)
{
  printf("ERROR-CUDA %s %s\n", comment, cudaGetErrorString(cudaGetLastError()));
}

int CudaWrapper_CheckMemUsage()
{
  size_t free, total;
  cudaMemGetInfo(&free, &total);
  return total - free; //possible with cuda 3.0 ???
  //return CudaWrapper_total_gpu_mem;
}

double CudaWrapper_CheckUploadTime(bool reset)
{
  if(reset) CudaWrapper_total_upload_time = 0.0;

  return CudaWrapper_total_upload_time;
}

double CudaWrapper_CheckDownloadTime(bool reset)
{
  if(reset) CudaWrapper_total_download_time = 0.0;

  return CudaWrapper_total_download_time;
}

double CudaWrapper_CheckCPUBufUploadTime(bool reset)
{
  if(reset) CudaWrapper_cpubuffer_upload_time = 0.0;

  return CudaWrapper_cpubuffer_upload_time;
}

double CudaWrapper_CheckCPUBufDownloadTime(bool reset)
{
  if(reset) CudaWrapper_cpubuffer_download_time = 0.0;

  return CudaWrapper_cpubuffer_download_time;
}

void CudaWrapper_AddCPUBufUploadTime(double dt)
{
  CudaWrapper_cpubuffer_upload_time += dt;
}

void CudaWrapper_AddCPUBufDownloadTime(double dt)
{
  CudaWrapper_cpubuffer_download_time += dt;
}

void CudaWrapper_Sync()
{
  cudaThreadSynchronize();
}

void CudaWrapper_SyncStream(int stream)
{
  cudaStreamSynchronize(streams[stream]);
}

void CudaWrapper_AddStreams(int n)
{
  cudaStream_t* new_streams = new cudaStream_t[nstreams + n];

  for(int i = 0; i < nstreams; i++) new_streams[i] = streams[i];

  for(int i = nstreams; i < nstreams + n; i++) cudaStreamCreate(&new_streams[i]);

  if(nstreams > 0)
    delete [] streams;

  streams = new_streams;
  nstreams += n;
}

void* CudaWrapper_returnStreams()
{
  return (void*) streams;
}

int CudaWrapper_returnNStreams()
{
  return nstreams;
}

