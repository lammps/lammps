#include "hip/hip_runtime.h"

#include "cuda_reduction.h"

#include "cuda_shuffle.h"
#include "cuda_utils.h"

#include "vector.h"


#include <hipcub/hipcub.hpp>

//struct RvecSum
//{
//    template <typename T>
//    __device__ __forceinline__
//    T operator()(const T &a, const T &b) const
//    {
//        T c;
//        c[0] = a[0] + b[0];
//        c[1] = a[1] + b[1];
//        c[2] = a[2] + b[2];
//        return c;
//    }
//};


/* Perform a device-wide reduction (sum operation)
 *
 * d_array: device array to reduce
 * d_dest: device pointer to hold result of reduction */
void Cuda_Reduction_Sum( int *d_array, int *d_dest, size_t n )
{
    void *d_temp_storage = NULL;
    size_t temp_storage_bytes = 0;

    /* determine temporary device storage requirements */
    hipcub::DeviceReduce::Sum( d_temp_storage, temp_storage_bytes,
            d_array, d_dest, n );
    hipDeviceSynchronize( );
    cudaCheckError( );

    /* allocate temporary storage */
    cuda_malloc( &d_temp_storage, temp_storage_bytes, FALSE,
            "Cuda_Reduction_Sum::temp_storage" );

    /* run sum-reduction */
    hipcub::DeviceReduce::Sum( d_temp_storage, temp_storage_bytes,
            d_array, d_dest, n );
    hipDeviceSynchronize( );
    cudaCheckError( );

    /* deallocate temporary storage */
    cuda_free( d_temp_storage, "Cuda_Reduction_Sum::temp_storage" );
}


/* Perform a device-wide reduction (sum operation)
 *
 * d_array: device array to reduce
 * d_dest: device pointer to hold result of reduction */
void Cuda_Reduction_Sum( real *d_array, real *d_dest, size_t n )
{
    void *d_temp_storage = NULL;
    size_t temp_storage_bytes = 0;

    /* determine temporary device storage requirements */
    hipcub::DeviceReduce::Sum( d_temp_storage, temp_storage_bytes,
            d_array, d_dest, n );
    hipDeviceSynchronize( );
    cudaCheckError( );

    /* allocate temporary storage */
    cuda_malloc( &d_temp_storage, temp_storage_bytes, FALSE,
            "Cuda_Reduction_Sum::temp_storage" );

    /* run sum-reduction */
    hipcub::DeviceReduce::Sum( d_temp_storage, temp_storage_bytes,
            d_array, d_dest, n );
    hipDeviceSynchronize( );
    cudaCheckError( );

    /* deallocate temporary storage */
    cuda_free( d_temp_storage, "Cuda_Reduction_Sum::temp_storage" );
}


///* Perform a device-wide reduction (sum operation)
// *
// * d_array: device array to reduce
// * d_dest: device pointer to hold result of reduction */
//void Cuda_Reduction_Sum( rvec *d_array, rvec *d_dest, size_t n )
//{
//    void *d_temp_storage = NULL;
//    size_t temp_storage_bytes = 0;
//    RvecSum sum_op;
//    rvec init = {0.0, 0.0, 0.0};
//
//    /* determine temporary device storage requirements */
//    cub::DeviceReduce::Reduce( d_temp_storage, temp_storage_bytes,
//            d_array, d_dest, n, sum_op, init );
//    hipDeviceSynchronize( );
//    cudaCheckError( );
//
//    /* allocate temporary storage */
//    cuda_malloc( &d_temp_storage, temp_storage_bytes, FALSE,
//            "cub::reduce::temp_storage" );
//
//    /* run sum-reduction */
//    cub::DeviceReduce::Reduce( d_temp_storage, temp_storage_bytes,
//            d_array, d_dest, n, sum_op, init );
//    hipDeviceSynchronize( );
//    cudaCheckError( );
//
//    /* deallocate temporary storage */
//    cuda_free( d_temp_storage, "cub::reduce::temp_storage" );
//}


/* Perform a device-wide reduction (max operation)
 *
 * d_array: device array to reduce
 * d_dest: device pointer to hold result of reduction */
void Cuda_Reduction_Max( int *d_array, int *d_dest, size_t n )
{
    void *d_temp_storage = NULL;
    size_t temp_storage_bytes = 0;

    /* determine temporary device storage requirements */
    hipcub::DeviceReduce::Max( d_temp_storage, temp_storage_bytes,
            d_array, d_dest, n );
    hipDeviceSynchronize( );
    cudaCheckError( );

    /* allocate temporary storage */
    cuda_malloc( &d_temp_storage, temp_storage_bytes, FALSE,
            "Cuda_Reduction_Max::temp_storage" );

    /* run exclusive prefix sum */
    hipcub::DeviceReduce::Max( d_temp_storage, temp_storage_bytes,
            d_array, d_dest, n );
    hipDeviceSynchronize( );
    cudaCheckError( );

    /* deallocate temporary storage */
    cuda_free( d_temp_storage, "Cuda_Reduction_Max::temp_storage" );
}


/* Perform a device-wide scan (partial sum operation)
 *
 * d_src: device array to scan
 * d_dest: device array to hold result of scan */
void Cuda_Scan_Excl_Sum( int *d_src, int *d_dest, size_t n )
{
    void *d_temp_storage = NULL;
    size_t temp_storage_bytes = 0;

    /* determine temporary device storage requirements */
    hipcub::DeviceScan::ExclusiveSum( d_temp_storage, temp_storage_bytes,
            d_src, d_dest, n );
    hipDeviceSynchronize( );
    cudaCheckError( );

    /* allocate temporary storage */
    cuda_malloc( &d_temp_storage, temp_storage_bytes, FALSE,
            "Cuda_Scan_Excl_Sum::temp_storage" );

    /* run exclusive prefix sum */
    hipcub::DeviceScan::ExclusiveSum( d_temp_storage, temp_storage_bytes,
            d_src, d_dest, n );
    hipDeviceSynchronize( );
    cudaCheckError( );

    /* deallocate temporary storage */
    cuda_free( d_temp_storage, "Cuda_Scan_Excl_Sum::temp_storage" );
}


CUDA_GLOBAL void k_reduction( const real *input, real *per_block_results,
        const size_t n )
{
#if defined(__SM_35__)
    HIP_DYNAMIC_SHARED( real, my_results)
    real sdata;
    unsigned int i;
    int z, offset;
    real x;

    i = blockIdx.x * blockDim.x + threadIdx.x;

    if ( i < n )
    {
        x = input[i];
    }
    else
    {
        x = 0.0;
    }

    sdata = x;
    __syncthreads( );

    for ( z = 16; z >= 1; z /= 2 )
    {
        sdata += shfl( sdata, z );
    }

    if ( threadIdx.x % 32 == 0 )
    {
        my_results[threadIdx.x >> 5] = sdata;
    }
    __syncthreads( );

    for ( offset = blockDim.x >> 6; offset > 0; offset >>= 1 )
    {
        if ( threadIdx.x < offset )
        {
            my_results[threadIdx.x] += my_results[threadIdx.x + offset];
        }
        __syncthreads( );
    }

    if ( threadIdx.x == 0 )
    {
        per_block_results[blockIdx.x] = my_results[0];
    }

#else
    HIP_DYNAMIC_SHARED( real, sdata)
    unsigned int i;
    int offset;
    real x;

    i = blockIdx.x * blockDim.x + threadIdx.x;

    if ( i < n )
    {
        x = input[i];
    }
    else
    {
        x = 0.0;
    }
    sdata[threadIdx.x] = x;
    __syncthreads( );

    for ( offset = blockDim.x / 2; offset > 0; offset >>= 1 )
    {
        if ( threadIdx.x < offset )
        {
            sdata[threadIdx.x] += sdata[threadIdx.x + offset];
        }

        __syncthreads( );
    }

    if ( threadIdx.x == 0 )
    {
        per_block_results[blockIdx.x] = sdata[0];
    }
#endif
}


CUDA_GLOBAL void k_reduction_rvec( rvec *input, rvec *results, size_t n )
{
#if defined(__SM_35__)
    HIP_DYNAMIC_SHARED( rvec, my_rvec)
    rvec sdata;
    unsigned int i;
    int z, offset;

    i = blockIdx.x * blockDim.x + threadIdx.x, z, offset;

    rvec_MakeZero( sdata );

    if ( i < n )
    {
        rvec_Copy( sdata, input[i] );
    }

    __syncthreads( );

    for ( z = 16; z >= 1; z /= 2 )
    {
        sdata[0] += shfl( sdata[0], z );
        sdata[1] += shfl( sdata[1], z );
        sdata[2] += shfl( sdata[2], z );
    }

    if ( threadIdx.x % 32 == 0 )
    {
        rvec_Copy( my_rvec[threadIdx.x >> 5], sdata );
    }

    __syncthreads( );

    for ( offset = blockDim.x >> 6; offset > 0; offset >>= 1 )
    {
        if ( threadIdx.x < offset )
        {
            rvec_Add( my_rvec[threadIdx.x], my_rvec[threadIdx.x + offset] );
        }

        __syncthreads( );
    }

    if ( threadIdx.x == 0 )
    {
        rvec_Add( results[blockIdx.x], my_rvec[0] );
    }

#else
    HIP_DYNAMIC_SHARED( rvec, svec_data)
    unsigned int i;
    int offset;
    rvec x;

    i = blockIdx.x * blockDim.x + threadIdx.x;

    rvec_MakeZero( x );

    if ( i < n )
    {
        rvec_Copy( x, input[i] );
    }

    rvec_Copy( svec_data[threadIdx.x], x );
    __syncthreads( );

    for ( offset = blockDim.x / 2; offset > 0; offset >>= 1 )
    {
        if ( threadIdx.x < offset )
        {
            rvec_Add( svec_data[threadIdx.x], svec_data[threadIdx.x + offset] );
        }

        __syncthreads( );
    }

    if ( threadIdx.x == 0 )
    {
        //rvec_Copy( results[blockIdx.x], svec_data[0] );
        rvec_Add( results[blockIdx.x], svec_data[0] );
    }
#endif
}


CUDA_GLOBAL void k_reduction_rvec2( rvec2 *input, rvec2 *results, size_t n )
{
#if defined(__SM_35__)
    HIP_DYNAMIC_SHARED( rvec2, my_rvec2)
    rvec2 sdata;
    unsigned int i;
    int z, offset;

    i = blockIdx.x * blockDim.x + threadIdx.x;

    sdata[0] = 0.0;
    sdata[1] = 0.0;

    if ( i < n )
    {
        sdata[0] = input[i][0];
        sdata[1] = input[i][1];
    }

    __syncthreads( );

    for ( z = 16; z >= 1; z /= 2 )
    {
        sdata[0] += shfl( sdata[0], z );
        sdata[1] += shfl( sdata[1], z );
    }

    if ( threadIdx.x % 32 == 0 )
    {
        my_rvec2[threadIdx.x >> 5][0] = sdata[0];
        my_rvec2[threadIdx.x >> 5][1] = sdata[1];
    }

    __syncthreads( );

    for ( offset = blockDim.x >> 6; offset > 0; offset >>= 1 )
    {
        if ( threadIdx.x < offset )
        {
            my_rvec2[threadIdx.x][0] += my_rvec2[threadIdx.x + offset][0];
            my_rvec2[threadIdx.x][1] += my_rvec2[threadIdx.x + offset][1];
        }

        __syncthreads( );
    }

    if ( threadIdx.x == 0 )
    {
        results[blockIdx.x][0] = my_rvec2[0][0];
        results[blockIdx.x][1] = my_rvec2[0][1];
    }

#else
    HIP_DYNAMIC_SHARED( rvec2, svec2_data)
    unsigned int i;
    int offset;
    rvec2 x;

    i = blockIdx.x * blockDim.x + threadIdx.x;

    x[0] = 0.0;
    x[1] = 0.0;

    if ( i < n )
    {
        x[0] += input[i][0];
        x[1] += input[i][1];
    }

    svec2_data[threadIdx.x][0] = x[0];
    svec2_data[threadIdx.x][1] = x[1];
    __syncthreads( );

    for ( offset = blockDim.x / 2; offset > 0; offset >>= 1 )
    {
        if ( threadIdx.x < offset )
        {
            svec2_data[threadIdx.x][0] += svec2_data[threadIdx.x + offset][0];
            svec2_data[threadIdx.x][1] += svec2_data[threadIdx.x + offset][1];
        }

        __syncthreads( );
    }

    if ( threadIdx.x == 0 )
    {
        //rvec_Copy( results[blockIdx.x], svec_data[0] );
        results[blockIdx.x][0] += svec2_data[0][0];
        results[blockIdx.x][1] += svec2_data[0][1];
    }
#endif
}


CUDA_GLOBAL void k_dot( const real *a, const real *b, real *per_block_results,
        const size_t n )
{
#if defined(__SM_35__)
    HIP_DYNAMIC_SHARED( real, my_dot)
    real sdot;
    unsigned int i;
    int z, offset;

    i = blockIdx.x * blockDim.x + threadIdx.x;

    if ( i < n )
    {
        sdot = a[i] * b[i];
    }
    else
    {
        sdot = 0.0;
    }

    __syncthreads( );

    for ( z = 16; z >= 1; z /= 2 )
    {
        sdot += shfl( sdot, z );
    }

    if ( threadIdx.x % 32 == 0 )
    {
        my_dot[threadIdx.x >> 5] = sdot;
    }

    __syncthreads( );

    for ( offset = blockDim.x >> 6; offset > 0; offset >>= 1 )
    {
        if ( threadIdx.x < offset )
        {
            my_dot[threadIdx.x] += my_dot[threadIdx.x + offset];
        }

        __syncthreads( );
    }

    if ( threadIdx.x == 0 )
    {
        per_block_results[blockIdx.x] = my_dot[0];
    }

#else
    HIP_DYNAMIC_SHARED( real, sdot)
    unsigned int i;
    int offset;
    real x;

    i = blockIdx.x * blockDim.x + threadIdx.x;

    if ( i < n )
    {
        x = a[i] * b[i];
    }
    else
    {
        x = 0.0;
    }
    sdot[threadIdx.x] = x;
    __syncthreads( );

    for ( offset = blockDim.x / 2; offset > 0; offset >>= 1 )
    {
        if ( threadIdx.x < offset )
        {
            sdot[threadIdx.x] += sdot[threadIdx.x + offset];
        }

        __syncthreads( );
    }

    if ( threadIdx.x == 0 )
    {
        per_block_results[blockIdx.x] = sdot[0];
    }

#endif
}


CUDA_GLOBAL void k_norm( const real *input, real *per_block_results,
        const size_t n, int pass )
{
#if defined(__SM_35__)
    HIP_DYNAMIC_SHARED( real, my_norm)
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    real snorm = 0.0;

    if( i < n )
    {
        snorm = SQR( input[i] );
    }

    __syncthreads();

    for(int z = 16; z >=1; z/=2)
    {
        snorm += shfl ( snorm, z );
    }

    if (threadIdx.x % 32 == 0)
    {
        my_norm[threadIdx.x >> 5] = snorm;
    }

    __syncthreads( );

    for(int offset = blockDim.x >> 6; offset > 0; offset >>= 1)
    {
        if(threadIdx.x < offset)
        {
            my_norm[threadIdx.x] += my_norm[threadIdx.x + offset];
        }

        __syncthreads();
    }

    if(threadIdx.x == 0)
    {
        per_block_results[blockIdx.x] = my_norm[0];
    }

#else
    HIP_DYNAMIC_SHARED( real, snorm)
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    real x = 0;

    if(i < n)
    {
        x = SQR( input[i] );
    }

    snorm[threadIdx.x] = x;
    __syncthreads();

    for(int offset = blockDim.x / 2; offset > 0; offset >>= 1)
    {
        if(threadIdx.x < offset)
        {
            snorm[threadIdx.x] += snorm[threadIdx.x + offset];
        }

        __syncthreads();
    }

    if(threadIdx.x == 0)
    {
        per_block_results[blockIdx.x] = snorm[0];
    }
#endif
}


CUDA_GLOBAL void k_norm_rvec2( const rvec2 *input, rvec2 *per_block_results,
        const size_t n, int pass )
{
#if defined(__SM_35__)
    HIP_DYNAMIC_SHARED( rvec2, my_norm2)
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    rvec2 snorm2;
    snorm2[0] = snorm2[1] = 0;

    if(i < n)
    {
        if (pass == INITIAL)
        {
            snorm2[0] = SQR( input[i][0] );
            snorm2[1] = SQR( input[i][1] );
        }
        else
        {
            snorm2[0] = input[i][0];
            snorm2[1] = input[i][1];
        }
    }
    __syncthreads( );

    for(int z = 16; z >=1; z/=2)
    {
        snorm2[0] += shfl( snorm2[0], z );
        snorm2[1] += shfl( snorm2[1], z );
    }

    if (threadIdx.x % 32 == 0){
        my_norm2[threadIdx.x >> 5][0] = snorm2[0];
        my_norm2[threadIdx.x >> 5][1] = snorm2[1];
    }

    __syncthreads( );

    for(int offset = blockDim.x >> 6; offset > 0; offset >>= 1)
    {
        if(threadIdx.x < offset)
        {
            my_norm2[threadIdx.x][0] += my_norm2[threadIdx.x + offset][0];
            my_norm2[threadIdx.x][1] += my_norm2[threadIdx.x + offset][1];
        }

        __syncthreads( );
    }

    if(threadIdx.x == 0)
    {
        per_block_results[blockIdx.x][0] = my_norm2[0][0];
        per_block_results[blockIdx.x][1] = my_norm2[0][1];
    }

#else
    HIP_DYNAMIC_SHARED( rvec2, snorm2)
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    rvec2 x;
    x[0] = x[1] = 0;

    if( i < n )
    {
        if( pass == INITIAL )
        {
            x[0] = SQR( input[i][0] );
            x[1] = SQR( input[i][1] );
        }
        else
        {
            x[0] = input[i][0];
            x[1] = input[i][1];
        }
    }

    snorm2[threadIdx.x][0] = x[0];
    snorm2[threadIdx.x][1] = x[1];
    __syncthreads( );

    for(int offset = blockDim.x / 2; offset > 0; offset >>= 1)
    {
        if(threadIdx.x < offset)
        {
            snorm2[threadIdx.x][0] += snorm2[threadIdx.x + offset][0];
            snorm2[threadIdx.x][1] += snorm2[threadIdx.x + offset][1];
        }

        __syncthreads( );
    }

    if(threadIdx.x == 0)
    {
        per_block_results[blockIdx.x][0] = snorm2[0][0];
        per_block_results[blockIdx.x][1] = snorm2[0][1];
    }
#endif
}


CUDA_GLOBAL void k_dot_rvec2( const rvec2 *a, rvec2 *b, rvec2 *res,
        const size_t n )
{
#if defined(__SM_35__)
    HIP_DYNAMIC_SHARED( rvec2, my_dot2)
    rvec2 sdot2;

    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    sdot2[0] = sdot2[1] = 0;

    if(i < n) {
        sdot2[0] = a[i][0] * b[i][0];
        sdot2[1] = a[i][1] * b[i][1];
    }

    __syncthreads();

    for(int z = 16; z >=1; z/=2){
        sdot2[0] += shfl ( sdot2[0], z);
        sdot2[1] += shfl ( sdot2[1], z);
    }

    if (threadIdx.x % 32 == 0){
        my_dot2[threadIdx.x >> 5][0] = sdot2[0];
        my_dot2[threadIdx.x >> 5][1] = sdot2[1];
    }

    __syncthreads ();

    for(int offset = blockDim.x >> 6; offset > 0; offset >>= 1)
    {
        if(threadIdx.x < offset)
        {
            my_dot2[threadIdx.x][0] += my_dot2[threadIdx.x + offset][0];
            my_dot2[threadIdx.x][1] += my_dot2[threadIdx.x + offset][1];
        }

        __syncthreads();
    }

    if(threadIdx.x == 0) {
        res[blockIdx.x][0] = my_dot2[0][0];
        res[blockIdx.x][1] = my_dot2[0][1];
    }


#else
    HIP_DYNAMIC_SHARED( rvec2, sdot2)
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    rvec2 x;
    x[0] = x[1] = 0;

    if(i < n) {
        x[0] = a[i][0] * b[i][0];
        x[1] = a[i][1] * b[i][1];
    }

    sdot2[threadIdx.x][0] = x[0];
    sdot2[threadIdx.x][1] = x[1];
    __syncthreads();

    for(int offset = blockDim.x / 2; offset > 0; offset >>= 1)
    {
        if(threadIdx.x < offset)
        {
            sdot2[threadIdx.x][0] += sdot2[threadIdx.x + offset][0];
            sdot2[threadIdx.x][1] += sdot2[threadIdx.x + offset][1];
        }

        __syncthreads();
    }

    if(threadIdx.x == 0) {
        res[blockIdx.x][0] = sdot2[0][0];
        res[blockIdx.x][1] = sdot2[0][1];
    }
#endif
}


//////////////////////////////////////////////////
//vector functions
//////////////////////////////////////////////////
CUDA_GLOBAL void k_vector_sum( real* dest, real c, real* v, real d, real* y,
        int n )
{
    int i;

    i = blockIdx.x * blockDim.x + threadIdx.x;

    if ( i >= n )
    {
        return;
    }

    dest[i] = c * v[i] + d * y[i];

  //  printf("Index: %d, Dest %f, c %f, v %f, y %f\n",i, dest[i], c,v[i], y[i]);

}


CUDA_GLOBAL void k_vector_mul( real* dest, real* v, real* y, int n )
{
    int i;

    i = blockIdx.x * blockDim.x + threadIdx.x;

    if ( i >= n )
    {
        return;
    }

    dest[i] = v[i] * y[i];

}


CUDA_GLOBAL void k_rvec2_mul( rvec2* dest, rvec2* v, rvec2* y, int n )
{
    int i;

    i = blockIdx.x * blockDim.x + threadIdx.x;

    if ( i >= n )
    {
        return;
    }

    dest[i][0] = v[i][0] * y[i][0];
    dest[i][1] = v[i][1] * y[i][1];
}


CUDA_GLOBAL void k_rvec2_pbetad( rvec2 *dest, rvec2 *a, 
        real beta0, real beta1, rvec2 *b, int n )
{
    int i;

    i = blockIdx.x * blockDim.x + threadIdx.x;

    if ( i >= n )
    {
        return;
    }

    dest[i][0] = a[i][0] + beta0 * b[i][0];
    dest[i][1] = a[i][1] + beta1 * b[i][1];
}
