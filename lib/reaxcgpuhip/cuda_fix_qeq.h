#ifndef __CUDA_ALLOCATE_H_
#define __CUDA_ALLOCATE_H_

#include "reaxc_types.h"


#ifdef __cplusplus
extern "C"  {
#endif

void  CudaAllocateStorageForFixQeq(int nmax, int dual_enabled, fix_qeq_gpu *qeq_gpu);
void  CudaInitStorageForFixQeq(fix_qeq_gpu *qeq_gpu, double *Hdia_inv, double *b_s,double *b_t,double *b_prc,double *b_prm,double *s,double *t, int N);
void  Cuda_Calculate_H_Matrix(reax_list **gpu_lists,  reax_system *system,fix_qeq_gpu *qeq_gpu, control_params *control, int inum);
void  Cuda_Init_Taper(fix_qeq_gpu *qeq_gpu,double *Tap, int numTap);
void  Cuda_Allocate_Matrix( sparse_matrix *, int, int );
void Cuda_Deallocate_Matrix( sparse_matrix *H );
void Cuda_Init_Sparse_Matrix_Indices( reax_system *system, sparse_matrix *H);
void  Cuda_Init_Fix_Atoms(reax_system *system,fix_qeq_gpu *qeq_gpu);
void  Cuda_Init_Matvec_Fix(int nn, fix_qeq_gpu *qeq_gpu, reax_system *system)
void  Cuda_Copy_Pertype_Parameters_To_Device(double *chi,double *eta,double *gamma,int ntypes,fix_qeq_gpu *qeq_gpu);
void Cuda_Copy_From_Device_Comm_Fix(double *buf, double *x, int n, int offset);
void Cuda_Copy_To_Device_Comm_Fix(double *buf, double *x, int n, int offset);
void  Cuda_Sparse_Matvec_Compute(sparse_matrix *H,double *x, double *q, double *eta, reax_atom *d_fix_my_atoms, int nn, int NN);
void Cuda_Vector_Sum_Fix( real *, real, real *, real, real *, int );
void Cuda_CG_Preconditioner_Fix( real *, real *, real *, int );
void  Cuda_Copy_Vector_From_Device(real *host_vector, real *device_vector, int nn);
void Cuda_Calculate_Q(int nn,fix_qeq_gpu *qeq_gpu,int blocks);
void  Cuda_Parallel_Vector_Acc(int nn,double *x);
void  Cuda_UpdateQ_And_Copy_To_Device_Comm_Fix(double *buf,fix_qeq_gpu *qeq_gpu,int n);

#ifdef __cplusplus
}
#endif

#endif
