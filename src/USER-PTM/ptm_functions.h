#ifndef PTM_FUNCTIONS_H
#define PTM_FUNCTIONS_H

#include <stdint.h>
#include <stdbool.h>
#include "ptm_initialize_data.h"
#include "ptm_constants.h"


//------------------------------------
//    function declarations
//------------------------------------
#ifdef __cplusplus
extern "C" {
#endif


int ptm_index(	ptm_local_handle_t local_handle, int32_t flags, int num_points, double (*atomic_positions)[3], int32_t* atomic_numbers, bool topological_ordering,	//inputs
		int32_t* p_type, int32_t* p_alloy_type, double* p_scale, double* p_rmsd, double* q, double* F, double* F_res, double* U, double* P, int8_t* mapping, double* p_interatomic_distance, double* p_lattice_constant);	//outputs


#ifdef __cplusplus
}
#endif

#endif

