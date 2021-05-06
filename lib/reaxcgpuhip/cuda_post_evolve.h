
#ifndef __CUDA_POST_EVOLVE_H__
#define __CUDA_POST_EVOLVE_H__

#include "reaxc_types.h"


#ifdef __cplusplus
extern "C" {
#endif

void post_evolve_velocities( reax_system *, simulation_data * );

#ifdef __cplusplus
}
#endif


#endif
