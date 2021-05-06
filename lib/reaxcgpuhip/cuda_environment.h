
#ifndef __CUDA_ENVIRONMENT_H__
#define __CUDA_ENVIRONMENT_H__

#include "reaxc_types.h"


#ifdef __cplusplus
extern "C"  {
#endif

void Setup_Cuda_Environment( int, int, int );
void Cleanup_Cuda_Environment( );

#ifdef __cplusplus
}
#endif


#endif
