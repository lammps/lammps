#ifndef PODCOMMON_H
#define PODCOMMON_H

#include <string>

#define DDOT ddot_
#define DGEMV dgemv_
#define DGEMM dgemm_
#define DGETRF dgetrf_
#define DGETRI dgetri_
#define DSYEV dsyev_
#define DPOSV dposv_

#define PODMIN(a,b) ((a) < (b) ? (a) : (b))
#define PODMAX(a,b) ((a) > (b) ? (a) : (b))

#define CPUFREE(x)                    \
{                                     \
  if (x != NULL) {                    \
    free(x);                          \
    x = NULL;                         \
  }                                   \
}

extern "C" {
  double DNRM2(int*,double*,int*);
  double DDOT(int*,double*,int*,double*,int*);
  void DAXPY(int*,double*,double*,int*,double*,int*);
  void DGEMV(char*,int*,int*,double*,double*,int*,double*,int*,double*,double*,int*);
  void DGEMM(char*,char*,int*,int*,int*,double*,double*,int*,
       double*,int*,double*,double*,int*);
  void DGETRF(int*,int*,double*,int*,int*,int*);
  void DGETRI(int*,double*,int*,int*,double*,int*,int*);
  void DTRSM(char *, char*, char*, char *, int *, int *, double*, double*, int*,
       double*, int*);

  void DSYEV( char* jobz, char* uplo, int* n, double* a, int* lda,
    double* w, double* work, int* lwork, int* info );

  void DPOSV( char* uplo, int* n, int* nrhs, double* a, int* lda,
        double* b, int* ldb, int* info );
}
template <typename T> static void TemplateMalloc(T **data, int n, int backend)
{
  if (backend == 0)     // One thread CPU

    // allocate the memory on the CPU

    *data = (T *) malloc(n*sizeof(T));
  if (backend == 1)  // Open MP

    // allocate the memory on the CPU

    *data = (T *) malloc(n*sizeof(T));
#ifdef USE_CUDA
  if (backend == 2)  // CUDA C

    // allocate the memory on the GPU

    CUDA_CHECK( cudaMalloc( (void**)data, n * sizeof(T) ) );
#endif
}

template <typename T> static void TemplateCopytoDevice(T *d_data, T *h_data, int n, int backend)
{
  if (backend == 0)
    for (int i=0; i<n; i++) d_data[i] = h_data[i];
  if (backend == 1)
    for (int i=0; i<n; i++) d_data[i] = h_data[i];
#ifdef USE_CUDA
  if (backend == 2)  // CUDA
    cudaCopytoDevice(d_data, h_data, n);
#endif
#ifdef USE_HIP
  if (backend == 3)  // HIP
    hipCopytoDevice(d_data, h_data, n);
#endif
}

template <typename T> static void TemplateCopytoHost(T *h_data, T *d_data, int n, int backend)
{
  if (backend == 0)
    for (int i=0; i<n; i++) h_data[i] = d_data[i];
  if (backend == 1)
    for (int i=0; i<n; i++) h_data[i] = d_data[i];
#ifdef USE_CUDA
  if (backend == 2)  // CUDA
    cudaCopytoHost(h_data, d_data, n);
#endif
#ifdef USE_HIP
  if (backend == 3)  // HIP
    hipCopytoHost(h_data, d_data, n);
#endif
}

template <typename T> static void TemplateFree(T *data, int backend)
{
  if (backend == 0)
    CPUFREE(data);
  if (backend == 1)
    CPUFREE(data);
#ifdef USE_CUDA
  if (backend == 2)  // CUDA
    GPUFREE(data);
#endif
#ifdef USE_HIP
  if (backend == 3)  // HIP
    HIPFREE(data);
#endif
}

#endif
