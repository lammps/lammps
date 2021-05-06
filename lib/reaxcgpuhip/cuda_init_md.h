
#ifndef __CUDA_INIT_MD_H__
#define __CUDA_INIT_MD_H__

#include "reaxc_types.h"


#ifdef __cplusplus
extern "C" {
#endif

void Cuda_Initialize( reax_system*, control_params*, simulation_data*,
        storage*, reax_list**,reax_list*, output_controls*, mpi_datatypes* );
void Cuda_Write_Reax_Lists(reax_system *system, reax_list**, reax_list*);


#ifdef __cplusplus
}
#endif


#endif
