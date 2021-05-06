
#include "cuda_lookup.h"

#include "cuda_utils.h"

#include "index_utils.h"



CUDA_GLOBAL void kvdW_coulomb_energy(
		LR_lookup_table *t_LR,  int N, int num_atom_types)
{

	int i = blockIdx.x * blockDim.x + threadIdx.x;

	if ( i >= N )
	{
		return;
	}

	if (i >= 20)
	{
		return;
	}


	LR_lookup_table *t;
	t = &t_LR[3];

	//printf("t %f \n", t->inv_dx);


}

void copy_LR_table_to_device( reax_system *system, control_params *control,
		storage *workspace, int *aggregated )
{
	int i, j;
	int num_atom_types;
	LR_data *d_y;
	cubic_spline_coef *temp;

	num_atom_types = system->reax_param.num_atom_types;

	//printf("atom types %d \n",num_atom_types);

	cuda_malloc( (void **) &workspace->d_LR,
			sizeof(LR_lookup_table) * ( num_atom_types * num_atom_types ),
			TRUE, "LR_lookup:table" );
	copy_host_device( workspace->LR, workspace->d_LR,
			sizeof(LR_lookup_table) * (num_atom_types * num_atom_types),
			hipMemcpyHostToDevice, "LR_lookup:table");

	/*copy_host_device( workspace->LR, workspace->d_LR,
			sizeof(LR_lookup_table) * (num_atom_types * num_atom_types),
			hipMemcpyDeviceToHost, "LR_lookup:table");



	for( i = 0; i < num_atom_types; ++i ) {
		if (aggregated[i]) {
			for( j = i; j < num_atom_types; ++j ) {
				if (aggregated[j]) {
					printf("%d,%d\n",i,j);
					printf("%f\n",workspace->LR[i][j].inv_dx);

				}
			}
		}
	}*/

	/*int blocks = ((system->N * VDW_KER_THREADS_PER_ATOM) / DEF_BLOCK_SIZE)
		        																+ (((system->N * VDW_KER_THREADS_PER_ATOM) % DEF_BLOCK_SIZE == 0) ? 0 : 1);

	hipLaunchKernelGGL(kvdW_coulomb_energy, dim3(blocks), dim3(DEF_BLOCK_SIZE ), 0, 0,
			workspace->d_LR, system->N,
			system->reax_param.num_atom_types);
	hipDeviceSynchronize();

	exit(0);*/




	for( i = 0; i < num_atom_types; ++i )
	{
		if ( aggregated[i] )
		{
			for( j = i; j < num_atom_types; ++j )
			{
				if ( aggregated[j] )
				{
					cuda_malloc((void **) &d_y,
							sizeof(LR_data) * (control->tabulate + 1) , FALSE, "LR_lookup:d_y");
					copy_host_device( workspace->LR[index_lr(i, j, num_atom_types)].y, d_y,
							sizeof(LR_data) * (control->tabulate + 1), hipMemcpyHostToDevice, "LR_lookup:y");
					copy_host_device ( &d_y, &workspace->d_LR[ index_lr(i, j, num_atom_types) ].y,
							sizeof(LR_data *), hipMemcpyHostToDevice, "LR_lookup:y" );




					cuda_malloc( (void **) &temp, sizeof(cubic_spline_coef) * (control->tabulate + 1), FALSE, "LR_lookup:h" );
					copy_host_device( workspace->LR[index_lr(i, j, num_atom_types)].H, temp,
							sizeof(cubic_spline_coef) * (control->tabulate + 1), hipMemcpyHostToDevice, "LR_lookup:h" );
					copy_host_device( &temp, &workspace->d_LR[ index_lr(i, j, num_atom_types) ].H,
							sizeof(cubic_spline_coef *), hipMemcpyHostToDevice, "LR_lookup:h" );



					cuda_malloc( (void **) &temp, sizeof(cubic_spline_coef) * (control->tabulate + 1), FALSE, "LR_lookup:vdW" );
					copy_host_device( workspace->LR[index_lr(i, j, num_atom_types)].vdW, temp,
							sizeof(cubic_spline_coef) * (control->tabulate + 1), hipMemcpyHostToDevice, "LR_lookup:vdW" );
					copy_host_device( &temp, &workspace->d_LR[ index_lr(i, j, num_atom_types) ].vdW,
							sizeof(cubic_spline_coef *), hipMemcpyHostToDevice, "LR_lookup:vdW" );


					cuda_malloc( (void **) &temp, sizeof(cubic_spline_coef) * (control->tabulate + 1), FALSE, "LR_lookup:CEvd" );
					copy_host_device( workspace->LR[index_lr(i, j, num_atom_types)].CEvd, temp,
							sizeof(cubic_spline_coef) * (control->tabulate + 1), hipMemcpyHostToDevice, "LR_lookup:CEvd" );
					copy_host_device( &temp, &workspace->d_LR[ index_lr(i, j, num_atom_types) ].CEvd,
							sizeof(cubic_spline_coef *), hipMemcpyHostToDevice, "LR_lookup:CDvd");



					cuda_malloc( (void **) &temp, sizeof(cubic_spline_coef) * (control->tabulate + 1), FALSE, "LR_lookup:ele" );
					copy_host_device( workspace->LR[index_lr(i, j, num_atom_types)].ele, temp,
							sizeof(cubic_spline_coef) * (control->tabulate + 1), hipMemcpyHostToDevice, "LR_lookup:ele" );
					copy_host_device( &temp, &workspace->d_LR[ index_lr(i, j, num_atom_types) ].ele,
							sizeof(cubic_spline_coef *), hipMemcpyHostToDevice, "LR_lookup:ele" );

					cuda_malloc( (void **) &temp, sizeof(cubic_spline_coef) * (control->tabulate + 1), FALSE, "LR_lookup:ceclmb" );
					copy_host_device( workspace->LR[index_lr(i, j, num_atom_types)].CEclmb, temp,
							sizeof(cubic_spline_coef) * (control->tabulate + 1), hipMemcpyHostToDevice, "LR_lookup:ceclmb" );
					copy_host_device( &temp, &workspace->d_LR[ index_lr(i, j, num_atom_types) ].CEclmb,
							sizeof(cubic_spline_coef *), hipMemcpyHostToDevice, "LR_lookup:ceclmb" );

				}
			}
		}
	}





	fprintf( stderr, "Copy of the LR Lookup Table to the device complete ... \n" );

	//exit(0);
}

