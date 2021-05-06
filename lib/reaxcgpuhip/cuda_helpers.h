#ifndef __CUDA_HELPERS__
#define __CUDA_HELPERS__

#include "reaxc_types.h"


CUDA_DEVICE static inline int cuda_strcmp( const char * a,
        const char * b, int len )
{
    int i;
    char *src, *dst;

    src = (char *) a;
    dst = (char *) b;

    for ( i = 0; i < len; i++ )
    {
        if ( *dst == '\0' )
        {
            return FALSE;
        }

        if ( *src != *dst )
        {
            return TRUE;
        }

        src++;
        dst++;
    }

    return FALSE;
}


CUDA_DEVICE static inline real myatomicAdd( real* address, real val )
{
    unsigned long long int* address_as_ull =
        (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do
    {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val + __longlong_as_double(assumed)));
    }
    while (assumed != old);

    return __longlong_as_double(old);
}


CUDA_DEVICE static inline void atomic_rvecAdd( rvec ret, rvec v )
{
    myatomicAdd( &ret[0], v[0] );
    myatomicAdd( &ret[1], v[1] );
    myatomicAdd( &ret[2], v[2] );
}


CUDA_DEVICE static inline void atomic_rvecScaledAdd( rvec ret, real c, rvec v )
{
    myatomicAdd( &ret[0], c * v[0] );
    myatomicAdd( &ret[1], c * v[1] );
    myatomicAdd( &ret[2], c * v[2] );
}

#endif
