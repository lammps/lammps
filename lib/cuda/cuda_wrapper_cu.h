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

#ifndef _CUDA_DATA_WRAPPER_H_
#define _CUDA_DATA_WRAPPER_H_

extern "C" void  CudaWrapper_Init(int argc, char** argv, int me = 0, int ppn = 2, int* devicelist = NULL);
extern "C" void* CudaWrapper_AllocCudaData(unsigned nbytes);
extern "C" void  CudaWrapper_UploadCudaData(void* host_data, void* dev_data, unsigned nbytes);
extern "C" void  CudaWrapper_UploadCudaDataAsync(void* host_data, void* dev_data, unsigned nbytes, int stream_id);
extern "C" void  CudaWrapper_DownloadCudaData(void* host_data, void* dev_data, unsigned nbytes);
extern "C" void  CudaWrapper_DownloadCudaDataAsync(void* host_data, void* dev_data, unsigned nbytes, int stream_id);
extern "C" void  CudaWrapper_FreeCudaData(void* dev_data, unsigned nbytes = 0);
extern "C" void  CudaWrapper_Memset(void* dev_data, int value, unsigned nbytes);
extern "C" void  CudaWrapper_CopyData(void* dev_dest, void* dev_source, unsigned nbytes);
extern "C" void* CudaWrapper_AllocPinnedHostData(unsigned nbytes, bool mapped = false, bool writeCombind = false);
extern "C" void  CudaWrapper_FreePinnedHostData(void* dev_data);
extern "C" void  cuda_check_error(char* comment);
extern "C" int   CudaWrapper_CheckMemUsage();
extern "C" double CudaWrapper_CheckUploadTime(bool reset = false);
extern "C" double CudaWrapper_CheckDownloadTime(bool reset = false);
extern "C" double CudaWrapper_CheckCPUBufUploadTime(bool reset = false);
extern "C" double CudaWrapper_CheckCPUBufDownloadTime(bool reset = false);
extern "C" void CudaWrapper_AddCPUBufUploadTime(double dt);
extern "C" void CudaWrapper_AddCPUBufDownloadTime(double dt);
extern "C" void CudaWrapper_Sync();
extern "C" void CudaWrapper_SyncStream(int n);
extern "C" void CudaWrapper_AddStreams(int n);
extern "C" void* CudaWrapper_returnStreams();
extern "C" int CudaWrapper_returnNStreams();

#endif // _CUDA_DATA_WRAPPER_H_
