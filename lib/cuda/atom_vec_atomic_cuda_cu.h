#ifndef ATOM_VEC_ATOMIC_CUDA_CU_H_
#define ATOM_VEC_ATOMIC_CUDA_CU_H_

extern "C" void Cuda_AtomVecAtomicCuda_Init(cuda_shared_data* sdata);
extern "C" int Cuda_AtomVecAtomicCuda_PackExchangeList(cuda_shared_data* sdata, int n, int dim, void* buf_send);
extern "C" int Cuda_AtomVecAtomicCuda_PackExchange(cuda_shared_data* sdata, int nsend, void* buf_send, void* copylist);
extern "C" int Cuda_AtomVecAtomicCuda_UnpackExchange(cuda_shared_data* sdata, int nsend, void* buf_send, void* copylist);
extern "C" int Cuda_AtomVecAtomicCuda_PackBorder(cuda_shared_data* sdata, int nsend, int iswap, void* buf_send, int* pbc, int pbc_flag);
extern "C" int Cuda_AtomVecAtomicCuda_PackBorderVel(cuda_shared_data* sdata, int n, int iswap, void* buf_send, int* pbc, int pbc_flag);
extern "C" int Cuda_AtomVecAtomicCuda_PackBorder_Self(cuda_shared_data* sdata, int n, int iswap, int first, int* pbc, int pbc_flag);
extern "C" int Cuda_AtomVecAtomicCuda_PackBorderVel_Self(cuda_shared_data* sdata, int n, int iswap, int first, int* pbc, int pbc_flag);
extern "C" int Cuda_AtomVecAtomicCuda_UnpackBorder(cuda_shared_data* sdata, int n, int first, void* buf_recv);
extern "C" int Cuda_AtomVecAtomicCuda_UnpackBorderVel(cuda_shared_data* sdata, int n, int first, void* buf_recv);

#endif /*ATOM_VEC_ATOMIC2_CUDA_CU_H_*/
