
#include "cuda_post_evolve.h"

#include "cuda_utils.h"

#include "vector.h"


CUDA_GLOBAL void ker_post_evolve( reax_atom *my_atoms, 
        simulation_data *data, int n )
{
    rvec diff, cross;
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= n)
    {
        return;
    }

    //for( i = 0; i < system->n; i++ ) { 
    /* remove translational vel */
    rvec_ScaledAdd( my_atoms[i].v, -1., data->vcm );

    /* remove rotational */
    rvec_ScaledSum( diff, 1., my_atoms[i].x, -1., data->xcm );
    rvec_Cross( cross, data->avcm, diff );
    rvec_ScaledAdd( my_atoms[i].v, -1., cross );
    //}  
}


void post_evolve_velocities( reax_system *system, simulation_data *data )
{
    int blocks;

    blocks = system->n / DEF_BLOCK_SIZE + 
        ((system->n % DEF_BLOCK_SIZE) == 0 ? 0 : 1);
    hipLaunchKernelGGL(ker_post_evolve, dim3(blocks), dim3(DEF_BLOCK_SIZE ), 0, 0, system->d_my_atoms, (simulation_data *)data->d_simulation_data, system->n);
    hipDeviceSynchronize( );
    cudaCheckError( );
}
