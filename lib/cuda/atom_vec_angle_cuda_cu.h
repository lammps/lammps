#ifndef ATOM_VEC_ANGLE_CUDA_CU_H_
#define ATOM_VEC_ANGLE_CUDA_CU_H_

extern "C" void Cuda_AtomVecAngleCuda_Init(cuda_shared_data* sdata);
extern "C" int Cuda_AtomVecAngleCuda_PackExchangeList(cuda_shared_data* sdata, int n, int dim, void* buf_send);
extern "C" int Cuda_AtomVecAngleCuda_PackExchange(cuda_shared_data* sdata, int nsend, void* buf_send, void* copylist);
extern "C" int Cuda_AtomVecAngleCuda_UnpackExchange(cuda_shared_data* sdata, int nsend, void* buf_send, void* copylist);
extern "C" int Cuda_AtomVecAngleCuda_PackBorder(cuda_shared_data* sdata, int nsend, int iswap, void* buf_send, int* pbc, int pbc_flag);
extern "C" int Cuda_AtomVecAngleCuda_PackBorderVel(cuda_shared_data* sdata, int n, int iswap, void* buf_send, int* pbc, int pbc_flag);
extern "C" int Cuda_AtomVecAngleCuda_PackBorder_Self(cuda_shared_data* sdata, int n, int iswap, int first, int* pbc, int pbc_flag);
extern "C" int Cuda_AtomVecAngleCuda_PackBorderVel_Self(cuda_shared_data* sdata, int n, int iswap, int first, int* pbc, int pbc_flag);
extern "C" int Cuda_AtomVecAngleCuda_UnpackBorder(cuda_shared_data* sdata, int n, int first, void* buf_recv);
extern "C" int Cuda_AtomVecAngleCuda_UnpackBorderVel(cuda_shared_data* sdata, int n, int first, void* buf_recv);

#endif /*ATOM_VEC_ANGLE2_CUDA_CU_H_*/
