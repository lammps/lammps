#ifndef ATOM_VEC_CHARGE_CUDA_CU_H_
#define ATOM_VEC_CHARGE_CUDA_CU_H_

extern "C" void Cuda_AtomVecChargeCuda_Init(cuda_shared_data* sdata);
extern "C" int Cuda_AtomVecChargeCuda_PackExchangeList(cuda_shared_data* sdata, int n, int dim, void* buf_send);
extern "C" int Cuda_AtomVecChargeCuda_PackExchange(cuda_shared_data* sdata, int nsend, void* buf_send, void* copylist);
extern "C" int Cuda_AtomVecChargeCuda_UnpackExchange(cuda_shared_data* sdata, int nsend, void* buf_send, void* copylist);
extern "C" int Cuda_AtomVecChargeCuda_PackBorder(cuda_shared_data* sdata, int nsend, int iswap, void* buf_send, int* pbc, int pbc_flag);
extern "C" int Cuda_AtomVecChargeCuda_PackBorderVel(cuda_shared_data* sdata, int n, int iswap, void* buf_send, int* pbc, int pbc_flag);
extern "C" int Cuda_AtomVecChargeCuda_PackBorder_Self(cuda_shared_data* sdata, int n, int iswap, int first, int* pbc, int pbc_flag);
extern "C" int Cuda_AtomVecChargeCuda_PackBorderVel_Self(cuda_shared_data* sdata, int n, int iswap, int first, int* pbc, int pbc_flag);
extern "C" int Cuda_AtomVecChargeCuda_UnpackBorder(cuda_shared_data* sdata, int n, int first, void* buf_recv);
extern "C" int Cuda_AtomVecChargeCuda_UnpackBorderVel(cuda_shared_data* sdata, int n, int first, void* buf_recv);

#endif /*ATOM_VEC_CHARGE2_CUDA_CU_H_*/
