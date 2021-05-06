#ifndef __CUDA_UTILS_H_
#define __CUDA_UTILS_H_

#include "reaxc_types.h"


#ifdef __cplusplus
extern "C"  {
#endif

void cuda_malloc( void **, size_t, int, const char * );

void cuda_free( void *, const char * );

void cuda_memset( void *, int , size_t , const char * );

void copy_host_device( void *, void *, size_t, enum hipMemcpyKind, const char * );

void copy_device( void *, void *, size_t, const char * );

void compute_blocks( int *, int *, int );

void compute_matvec_blocks( int *, int );

void compute_nearest_pow_2( int, int * );

void Cuda_Init_Block_Sizes( reax_system *, control_params * );

void Cuda_Print_Mem_Usage( );


#ifdef __cplusplus
#define cudaCheckError()    __cudaCheckError( __FILE__, __LINE__ )
static inline void __cudaCheckError( const char *file, const int line )
{
    hipError_t err;

    err = hipGetLastError();
    if ( hipSuccess != err )
    {
        fprintf( stderr, "[ERROR] runtime error encountered: %s:%d\n", file, line );
        fprintf( stderr, "    [INFO] CUDA API error code: %d\n", err );
        exit( RUNTIME_ERROR );
    }

#if defined(DEBUG_FOCUS)
    /* More careful checking. However, this will affect performance. */
    err = hipDeviceSynchronize( );
    if( hipSuccess != err )
    {
       exit( RUNTIME_ERROR );
    }
#endif

    return;
}
#endif

#ifdef __cplusplus
}
#endif


#endif
