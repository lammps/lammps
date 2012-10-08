#ifndef CUDA_DATA_CU_H_
#define CUDA_DATA_CU_H_

extern "C" void CudaData_Upload_DoubleFloat(void* host_data, void* dev_data, unsigned* n, copy_mode mode, void* buffer);
extern "C" void CudaData_Upload_DoubleDouble(void* host_data, void* dev_data, unsigned* n, copy_mode mode, void* buffer);
extern "C" void CudaData_Upload_FloatDouble(void* host_data, void* dev_data, unsigned* n, copy_mode mode, void* buffer);
extern "C" void CudaData_Upload_FloatFloat(void* host_data, void* dev_data, unsigned* n, copy_mode mode, void* buffer);
extern "C" void CudaData_Upload_IntInt(void* host_data, void* dev_data, unsigned* n, copy_mode mode, void* buffer);

extern "C" void CudaData_Download(void* host_data, void* dev_data, int host_size, int dev_size, unsigned* n, copy_mode mode, void* buffer);


#endif /*CUDA_DATA_CU_H_*/
