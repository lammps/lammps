
#ifndef __CUDA_LOOKUP_H__
#define __CUDA_LOOKUP_H__

#include "reaxc_types.h"


#ifdef __cplusplus
extern "C" {
#endif

void copy_LR_table_to_device( reax_system *, control_params *,
        storage *, int * );

#ifdef __cplusplus
}
#endif


#endif
