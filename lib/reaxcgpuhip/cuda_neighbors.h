
#ifndef __CUDA_NEIGHBORS_H__
#define __CUDA_NEIGHBORS_H__

#include "reaxc_types.h"


#ifdef __cplusplus
extern "C" {
#endif

int Cuda_Generate_Neighbor_Lists( reax_system *, simulation_data *, storage *, reax_list ** );

void Cuda_Estimate_Neighbors( reax_system * );

#ifdef __cplusplus
}
#endif


#endif
