#include "cuda_allocate.h"
#include "cuda_forces.h"
#include "cuda_list.h"
#include "cuda_neighbors.h"
#include "cuda_utils.h"
#include "cuda_fix_qeq.h"
#include "cuda_reduction.h"
#include "hip/hip_runtime.h"
#include "index_utils.h"
#include "tool_box.h"
#include "vector.h"

extern "C"
{




CUDA_DEVICE real Init_Charge_Matrix_Entry(real *workspace_Tap,
		int i, int j, real r_ij, real gamma)
{
	real Tap,denom;

	Tap = workspace_Tap[7] * r_ij + workspace_Tap[6];
	Tap = Tap * r_ij + workspace_Tap[5];
	Tap = Tap * r_ij + workspace_Tap[4];
	Tap = Tap * r_ij + workspace_Tap[3];
	Tap = Tap * r_ij + workspace_Tap[2];
	Tap = Tap * r_ij + workspace_Tap[1];
	Tap = Tap * r_ij + workspace_Tap[0];

	denom = r_ij * r_ij * r_ij + gamma;
	denom = POW(denom, 1.0 / 3.0 );

	return Tap * EV_to_KCALpMOL / denom;
}


/* Compute the charge matrix entries and store the matrix in full format
 * using the far neighbors list (stored in full format) and according to
 * the full shell communication method */
CUDA_GLOBAL void k_init_cm_full_fs( reax_atom *my_atoms,
		reax_list far_nbrs_list, sparse_matrix H, float nonb_cut, int n, double *d_Tap, double *gamma, int small)
{



	int i, j, pj;
	int start_i, end_i;
	int type_i, type_j;
	int cm_top;
	int num_cm_entries;
	real r_ij;
	two_body_parameters *twbp;
	reax_atom *atom_i, *atom_j;
	far_neighbor_data *nbr_pj;
	double shld;
	double dx, dy, dz;

	int flag = 0;
	int flag3 = 0;

	i = blockIdx.x * blockDim.x + threadIdx.x;

	if ( i >= n)
	{
		return;
	}

	cm_top = H.start[i];

	atom_i = &my_atoms[i];
	type_i = atom_i->type;
	start_i = Cuda_Start_Index( i, &far_nbrs_list );
	end_i = Cuda_End_Index(i, &far_nbrs_list );



	for ( pj = start_i; pj < end_i; ++pj )
	{
		nbr_pj = &far_nbrs_list.select.far_nbr_list[pj];
		j = nbr_pj->nbr;
		atom_j = &my_atoms[j];
		type_j = atom_j->type;



		//if(i == 0)
			//	printf("%f\n", nbr_pj->d);


		flag = 0;


		if ( nbr_pj->d  <= nonb_cut)
		{
			//if ( i == 500)
				//printf("%d,%d,%f,%f\n", i, j, nbr_pj->d,nonb_cut);
			r_ij =  nbr_pj->d;
			H.entries[cm_top].j = j;
			shld = pow( gamma[type_i] * gamma[type_j], -1.5);
			H.entries[cm_top].val = Init_Charge_Matrix_Entry(d_Tap,
					i, H.entries[cm_top].j, r_ij, shld);
			//printf("%d,%.3f\n",cm_top,H.entries[cm_top].val);
			++cm_top;

		}

	}
	__syncthreads();


	H.end[i] = cm_top;
	num_cm_entries = cm_top - H.start[i];


}
/* Compute the distances and displacement vectors for entries
 * in the far neighbors list if it's a NOT re-neighboring step */
CUDA_GLOBAL void k_init_distance( reax_atom *my_atoms, reax_list far_nbrs_list, int N )
{
	int i, j, pj;
	int start_i, end_i;
	reax_atom *atom_i, *atom_j;
	far_neighbor_data *nbr_pj;

	i = blockIdx.x * blockDim.x + threadIdx.x;

	if ( i >= N )
	{
		return;
	}

	atom_i = &my_atoms[i];
	start_i = Cuda_Start_Index( i, &far_nbrs_list );
	end_i =   Cuda_End_Index( i, &far_nbrs_list );

	/* update distance and displacement vector between atoms i and j (i-j) */
	for ( pj = start_i; pj < end_i; ++pj )
	{
		nbr_pj = &far_nbrs_list.select.far_nbr_list[pj];
		j = nbr_pj->nbr;
		atom_j = &my_atoms[j];

		if ( i < j )
		{
			nbr_pj->dvec[0] = atom_j->x[0] - atom_i->x[0];
			nbr_pj->dvec[1] = atom_j->x[1] - atom_i->x[1];
			nbr_pj->dvec[2] = atom_j->x[2] - atom_i->x[2];
		}
		else
		{
			nbr_pj->dvec[0] = atom_i->x[0] - atom_j->x[0];
			nbr_pj->dvec[1] = atom_i->x[1] - atom_j->x[1];
			nbr_pj->dvec[2] = atom_i->x[2] - atom_j->x[2];
		}

		//if(i == 0)
			//printf("%d,%f,%f,%f\n",j,atom_j->x[0],atom_j->x[1],atom_j->x[2]);
		nbr_pj->d = rvec_Norm(nbr_pj->dvec);
	}
}

void Cuda_Allocate_Hist_ST(fix_qeq_gpu *qeq_gpu,int nmax)
{
	cuda_malloc( (void **) &qeq_gpu->s_hist, sizeof(rvec4) * nmax, TRUE, "b" );
	cuda_malloc( (void **) &qeq_gpu->t_hist, sizeof(rvec4) * nmax, TRUE, "x" );
}


void  CudaAllocateStorageForFixQeq(int nmax, int dual_enabled, fix_qeq_gpu *qeq_gpu)
{
	cuda_malloc( (void **) &qeq_gpu->s, sizeof(double) * nmax, TRUE,
			"Cuda_Allocate_Matrix::start" );
	cuda_malloc( (void **) &qeq_gpu->t, sizeof(double) * nmax, TRUE,
			"Cuda_Allocate_Matrix::start" );
	cuda_malloc( (void **) &qeq_gpu->Hdia_inv, sizeof(double) * nmax, TRUE,
			"Cuda_Allocate_Matrix::start" );
	cuda_malloc( (void **) &qeq_gpu->b_s, sizeof(double) * nmax, TRUE,
			"Cuda_Allocate_Matrix::start" );
	cuda_malloc( (void **) &qeq_gpu->b_t, sizeof(double) * nmax, TRUE,
			"Cuda_Allocate_Matrix::start" );
	cuda_malloc( (void **) &qeq_gpu->b_prc, sizeof(double) * nmax, TRUE,
			"Cuda_Allocate_Matrix::start" );
	cuda_malloc( (void **) &qeq_gpu->b_prm, sizeof(double) * nmax, TRUE,
			"Cuda_Allocate_Matrix::start" );

	int size = nmax;
	if (dual_enabled)
	{
		size*= 2;
	}

	cuda_malloc( (void **) &qeq_gpu->p, sizeof(double) * size, TRUE,
			"Cuda_Allocate_Matrix::start" );
	cuda_malloc( (void **) &qeq_gpu->q, sizeof(double) * size, TRUE,
			"Cuda_Allocate_Matrix::start" );
	cuda_malloc( (void **) &qeq_gpu->r, sizeof(double) * size, TRUE,
			"Cuda_Allocate_Matrix::start" );
	cuda_malloc( (void **) &qeq_gpu->d, sizeof(double) * size, TRUE,
			"Cuda_Allocate_Matrix::start" );
}


void  CudaInitStorageForFixQeq(fix_qeq_gpu *qeq_gpu, double *Hdia_inv, double *b_s,double *b_t,double *b_prc,double *b_prm,double *s,double *t, int NN)
{
	copy_host_device( Hdia_inv, qeq_gpu->Hdia_inv, sizeof(double) * NN,
			hipMemcpyHostToDevice, "Cuda_CG::q:get" );
	copy_host_device( b_s, qeq_gpu->b_s, sizeof(double) * NN,
			hipMemcpyHostToDevice, "Cuda_CG::q:get" );
	copy_host_device( b_t, qeq_gpu->b_t, sizeof(double) * NN,
			hipMemcpyHostToDevice, "Cuda_CG::q:get" );
	copy_host_device( b_prc, qeq_gpu->b_prc, sizeof(double) * NN,
			hipMemcpyHostToDevice, "Cuda_CG::q:get" );
	copy_host_device( b_prm, qeq_gpu->b_prm, sizeof(double) * NN,
			hipMemcpyHostToDevice, "Cuda_CG::q:get" );
	copy_host_device( s, qeq_gpu->s, sizeof(double) * NN,
			hipMemcpyHostToDevice, "Cuda_CG::q:get" );
	copy_host_device(t, qeq_gpu->t, sizeof(double) * NN,
			hipMemcpyHostToDevice, "Cuda_CG::q:get");


}

void  Cuda_Init_Fix_Atoms(reax_system *system,fix_qeq_gpu *qeq_gpu)

{
	cuda_malloc( (void **) &qeq_gpu->d_fix_my_atoms, sizeof(reax_atom) * system->N, TRUE,
			"Cuda_Allocate_Matrix::start" );
	copy_host_device(qeq_gpu->fix_my_atoms, qeq_gpu->d_fix_my_atoms, sizeof(reax_atom) * system->N,
			hipMemcpyHostToDevice, "Sync_Atoms::system->my_atoms");
}

void Cuda_Free_Memory(fix_qeq_gpu *qeq_gpu)
{
	hipFree(qeq_gpu->d_fix_my_atoms);
	hipFree(qeq_gpu->d_cm_entries);
	hipFree(qeq_gpu->d_max_cm_entries);
	hipFree(qeq_gpu->d_total_cm_entries);
	hipFree(qeq_gpu->H.start);
	hipFree(qeq_gpu->H.end);
	hipFree(qeq_gpu->H.entries);

	//printf("Freeing memory\n");
}







void  Cuda_Calculate_H_Matrix(reax_list **lists,  reax_system *system, fix_qeq_gpu *qeq_gpu,control_params *control, int n, int small)
{

	int blocks;
	blocks = (n) / DEF_BLOCK_SIZE +
			(((n % DEF_BLOCK_SIZE) == 0) ? 0 : 1);

	hipLaunchKernelGGL(k_init_distance, dim3(blocks), dim3(DEF_BLOCK_SIZE), 0, 0,  qeq_gpu->d_fix_my_atoms, *(lists[FAR_NBRS]), n );
	hipDeviceSynchronize();

	//printf("nonb %f, %d\n",control->nonb_cut,n );


	//printf("Blocks %d , blocks size %d\n", blocks, DEF_BLOCK_SIZE);
	//printf("N %d, h n %d \n",system->N, qeq_gpu->H.n);

	hipLaunchKernelGGL(k_init_cm_full_fs , dim3(blocks), dim3(DEF_BLOCK_SIZE), 0, 0,  qeq_gpu->d_fix_my_atoms,
			*(lists[FAR_NBRS]),qeq_gpu->H, control->nonb_cut, n, qeq_gpu->d_Tap,qeq_gpu->gamma,small);
	hipDeviceSynchronize();


}

void Cuda_Init_Taper(fix_qeq_gpu *qeq_gpu,double *Tap, int numTap)
{
	cuda_malloc( (void **) &qeq_gpu->d_Tap, sizeof(double)*numTap, TRUE,
			"Cuda_Allocate_Matrix::start");
	copy_host_device(Tap, qeq_gpu->d_Tap, sizeof(double) * numTap,
			hipMemcpyHostToDevice, "Cuda_CG::q:get");

}


void Cuda_Allocate_Matrix( sparse_matrix *H, int n, int m )
{
	H->m = m;
	H->n = n;

	cuda_malloc( (void **) &H->start, sizeof(int) * n, TRUE,
			"Cuda_Allocate_Matrix::start" );
	cuda_malloc( (void **) &H->end, sizeof(int) * n, TRUE,
			"Cuda_Allocate_Matrix::end" );
	cuda_malloc( (void **) &H->entries, sizeof(sparse_matrix_entry) * m, TRUE,
			"Cuda_Allocate_Matrix::entries" );

}


void Cuda_Deallocate_Matrix( sparse_matrix *H )
{
	cuda_free( H->start, "Cuda_Deallocate_Matrix::start" );
	cuda_free( H->end, "Cuda_Deallocate_Matrix::end" );
	cuda_free( H->entries, "Cuda_Deallocate_Matrix::entries" );
}




CUDA_GLOBAL void k_estimate_cm_entries_storage(reax_atom *my_atoms,
		control_params *control, reax_list far_nbrs,
		int n,
		int *cm_entries, int *max_cm_entries)
{
	int i, j, pj;
	int start_i, end_i;
	int type_i, type_j;
	int local;
	int  num_cm_entries;
	real cutoff;
	far_neighbor_data *nbr_pj;
	reax_atom *atom_i, *atom_j;

	i = blockIdx.x * blockDim.x + threadIdx.x;

	if ( i >= n )
	{
		return;
	}


	num_cm_entries = 0;


	atom_i = &my_atoms[i];
	type_i = atom_i->type;
	start_i = Cuda_Start_Index( i, &far_nbrs );
	end_i = Cuda_End_Index( i, &far_nbrs );


	cutoff = control->nonb_cut;

    ++num_cm_entries;


	for ( pj = start_i; pj < end_i; ++pj )
	{
		nbr_pj = &far_nbrs.select.far_nbr_list[pj];
		if (nbr_pj->d <= control->nonb_cut)
		{
			 ++num_cm_entries;

		}
	}

    __syncthreads( );


	cm_entries[i] = num_cm_entries;
	max_cm_entries[i] = MAX( (int)(num_cm_entries * SAFE_ZONE), MIN_CM_ENTRIES );

}

void Cuda_Estimate_CMEntries_Storages( reax_system *system, control_params *control, reax_list **lists, fix_qeq_gpu *qeq_gpu,int n)
{
	int blocks;

	cuda_malloc( (void **) &qeq_gpu->d_cm_entries,
			n * sizeof(int), TRUE, "system:d_cm_entries" );
	cudaCheckError();

	cuda_malloc( (void **) &qeq_gpu->d_max_cm_entries,
			n * sizeof(int), TRUE, "system:d_cm_entries" );
	cuda_malloc( (void **) &qeq_gpu->d_total_cm_entries,
			n * sizeof(int), TRUE, "system:d_cm_entries" );


	blocks = n / DEF_BLOCK_SIZE +
			(((n % DEF_BLOCK_SIZE == 0)) ? 0 : 1);

	//printf("nn %d, sys n %d \n",n,system->n);

	hipLaunchKernelGGL(k_estimate_cm_entries_storage, dim3(blocks), dim3(DEF_BLOCK_SIZE), 0, 0,  qeq_gpu->d_fix_my_atoms,
			(control_params *)control->d_control_params,
			*(lists[FAR_NBRS]), n,
			qeq_gpu->d_cm_entries,qeq_gpu->d_max_cm_entries);
	hipDeviceSynchronize();
	cudaCheckError();


	//TB:: Should max_cm or cm entries be used for calculating total_cm_entries
	Cuda_Reduction_Sum(qeq_gpu->d_max_cm_entries, qeq_gpu->d_total_cm_entries, n);
	copy_host_device( &system->total_cm_entries, qeq_gpu->d_total_cm_entries, sizeof(int),
			hipMemcpyDeviceToHost, "Cuda_Estimate_Storages::d_total_cm_entries" );

	//printf("Total cm entries %d \n", system->total_cm_entries);
}

void Cuda_Init_Sparse_Matrix_Indices( reax_system *system, fix_qeq_gpu *qeq_gpu, int n)
{
	int blocks;

	/* init indices */
	Cuda_Scan_Excl_Sum(qeq_gpu->d_max_cm_entries, qeq_gpu->H.start, n);

	/* init end_indices */
	blocks = n / DEF_BLOCK_SIZE
			+ ((n % DEF_BLOCK_SIZE == 0) ? 0 : 1);
	hipLaunchKernelGGL(k_init_end_index, dim3(blocks), dim3(DEF_BLOCK_SIZE ), 0, 0,  qeq_gpu->d_cm_entries, qeq_gpu->H.start, qeq_gpu->H.end, n);
	hipDeviceSynchronize( );
	cudaCheckError();
}


CUDA_GLOBAL void k_init_matvec_fix(fix_qeq_gpu d_qeq_gpu,int n, single_body_parameters
		*sbp,reax_atom *my_atoms)
{
	int i;
	int type_i;
	fix_qeq_gpu *qeq_gpu;
	qeq_gpu = &d_qeq_gpu;


	i = blockIdx.x * blockDim.x + threadIdx.x;
	if ( i >= n)
	{
		return;
	}
	reax_atom *atom;
	atom = &my_atoms[i];
	type_i = atom->type;


	qeq_gpu->Hdia_inv[i] = 1. / qeq_gpu->eta[type_i];
	qeq_gpu->b_s[i] = -qeq_gpu->chi[type_i];
	qeq_gpu->b_t[i] = -1.0;



	qeq_gpu->t[i] = qeq_gpu->t_hist[i][2] + 3 * ( qeq_gpu->t_hist[i][0] - qeq_gpu->t_hist[i][1]);
	/* cubic extrapolation for s & t from previous solutions */
	qeq_gpu->s[i] = 4*(qeq_gpu->s_hist[i][0]+qeq_gpu->s_hist[i][2])-(6*qeq_gpu->s_hist[i][1]+qeq_gpu->s_hist[i][3]);
}

void  Cuda_Init_Matvec_Fix(int n, fix_qeq_gpu *qeq_gpu, reax_system *system)
{
	int blocks;

	blocks = n / DEF_BLOCK_SIZE
			+ (( n % DEF_BLOCK_SIZE == 0 ) ? 0 : 1);

	hipLaunchKernelGGL(k_init_matvec_fix, dim3(blocks), dim3(DEF_BLOCK_SIZE ), 0, 0, *(qeq_gpu),n,system->reax_param.d_sbp,qeq_gpu->d_fix_my_atoms);
	hipDeviceSynchronize();
	cudaCheckError();
}

void  Cuda_Copy_Pertype_Parameters_To_Device(double *chi,double *eta,double *gamma,int ntypes,fix_qeq_gpu *qeq_gpu)
{
	cuda_malloc( (void **) &qeq_gpu->gamma, sizeof(double)*(ntypes+1), TRUE,
			"Cuda_Allocate_Matrix::start");
	copy_host_device(gamma, qeq_gpu->gamma, sizeof(double) * (ntypes+1),
			hipMemcpyHostToDevice, "Cuda_CG::q:get");
	cuda_malloc( (void **) &qeq_gpu->chi, sizeof(double)*(ntypes+1), TRUE,
			"Cuda_Allocate_Matrix::start");
	copy_host_device(chi, qeq_gpu->chi, sizeof(double) * (ntypes+1),
			hipMemcpyHostToDevice, "Cuda_CG::q:get");
	cuda_malloc( (void **) &qeq_gpu->eta, sizeof(double)*(ntypes+1), TRUE,
			"Cuda_Allocate_Matrix::start");
	copy_host_device(eta, qeq_gpu->eta, sizeof(double) * (ntypes+1),
			hipMemcpyHostToDevice, "Cuda_CG::q:get");

}

void  Cuda_Copy_From_Device_Comm_Fix(double *buf, double *x, int n, int offset)
{
	copy_host_device(buf, x+offset, sizeof(double) * n,
			hipMemcpyDeviceToHost, "Cuda_CG::x:get" );
	//printf("Copy from device  fix \n");
}


CUDA_GLOBAL void k_update_buf(double *dev_buf, double *x, int nn, int offset)
{
	int i, c, col;

	i = blockIdx.x * blockDim.x + threadIdx.x;

	if ( i >= nn )
	{
		return;
	}


	x[i+offset] = dev_buf[i];
}


void  Cuda_Copy_To_Device_Comm_Fix(double *buf,double *x,int nn,int offset)
{

	double *dev_buf;
	cuda_malloc( (void **) &dev_buf, sizeof(double)*nn, TRUE,
			"Cuda_Allocate_Matrix::start");
	copy_host_device(buf,dev_buf,sizeof(double)*nn,hipMemcpyHostToDevice, "Cuda_CG::x:get");


	int blocks;

	blocks = nn / DEF_BLOCK_SIZE
			+ (( nn % DEF_BLOCK_SIZE == 0 ) ? 0 : 1);
	//printf("Blocks %d \n",blocks);

	hipLaunchKernelGGL(k_update_buf, dim3(blocks), dim3(DEF_BLOCK_SIZE ), 0, 0, dev_buf,x,nn, offset);
	hipDeviceSynchronize();

	hipFree(dev_buf);

	/*copy_host_device(buf, x+offset, sizeof(double) * nn,
				hipMemcpyHostToDevice, "Cuda_CG::x:get" );*/

}


CUDA_GLOBAL void k_update_q(double *temp_buf, double *q, int nn)
{
	int i, c, col;

	i = blockIdx.x * blockDim.x + threadIdx.x;

	if ( i >= nn )
	{
		return;
	}


	q[i] = q[i] +  temp_buf[i];

	//printf("m: %d %f\n",i, q[i]);

}

void  Cuda_UpdateQ_And_Copy_To_Device_Comm_Fix(double *buf,fix_qeq_gpu *qeq_gpu,int nn)
{
	double *temp_buf;
	cuda_malloc( (void **) &temp_buf, sizeof(double)*nn, TRUE,
			"Cuda_Allocate_Matrix::start");
	copy_host_device(buf, temp_buf, sizeof(double) * nn,
			hipMemcpyHostToDevice, "Cuda_CG::q:get");

	int blocks;

	blocks = nn / DEF_BLOCK_SIZE
			+ (( nn % DEF_BLOCK_SIZE == 0 ) ? 0 : 1);
	//printf("Blocks %d \n",blocks);


	hipLaunchKernelGGL(k_update_q, dim3(blocks), dim3(DEF_BLOCK_SIZE ), 0, 0, temp_buf,qeq_gpu->q,nn);
	hipDeviceSynchronize();

}



CUDA_GLOBAL void k_matvec_csr_fix( sparse_matrix H, real *vec, real *results,
		int num_rows)
{

	int i, c, col;
	real results_row;
	real val;

	i = blockIdx.x * blockDim.x + threadIdx.x;

	if ( i >= num_rows )
	{
		return;
	}



	results_row = results[i];


	int iter = 0;
	for ( c = H.start[i]; c < H.end[i]; c++ )
	{
		col = H.entries [c].j;
		val = H.entries[c].val;
		results_row += val * vec[col];
		iter++;


	}

	__syncthreads();


	results[i] = results_row;
}

CUDA_GLOBAL void k_init_q(reax_atom *my_atoms, double *q, double *x,double *eta, int nn, int NN)
{

	int i;
	int type_i;
	reax_atom *atom;


	i = blockIdx.x * blockDim.x + threadIdx.x;

	if ( i >= NN)
	{
		return;
	}

	if (i < nn ) {
		atom = &my_atoms[i];
		type_i = atom->type;


		q[i] = eta[type_i] * x[i];

		//printf("i %d, eta %f, x %f, q%f\n ", i, eta[type_i],x[i],q[i]);
	}
	else
	{
		q[i] = 0.0;

	}



}


void Cuda_Sparse_Matvec_Compute(sparse_matrix *H,double *x, double *q, double *eta, reax_atom *d_fix_my_atoms, int nn, int NN)
{

	int blocks;

	blocks = NN / DEF_BLOCK_SIZE
			+ (( NN % DEF_BLOCK_SIZE == 0 ) ? 0 : 1);
	//printf("Blocks %d \n",blocks);


	hipLaunchKernelGGL(k_init_q, dim3(blocks), dim3(DEF_BLOCK_SIZE ), 0, 0, d_fix_my_atoms,q,x,eta,nn,NN);
	hipDeviceSynchronize();



	blocks = nn / DEF_BLOCK_SIZE
			+ (( nn % DEF_BLOCK_SIZE == 0 ) ? 0 : 1);
	//printf("Blocks %d \n",blocks);


	hipLaunchKernelGGL(k_matvec_csr_fix, dim3(blocks), dim3(DEF_BLOCK_SIZE), 0 , 0, *H, x, q, nn);
	hipDeviceSynchronize();
	cudaCheckError();

	//printf("\n\n");
}

void Cuda_Vector_Sum_Fix( real *res, real a, real *x, real b, real *y, int count )
{
	//res = ax + by
	//use the cublas here
	int blocks;

	blocks = (count / DEF_BLOCK_SIZE)+ ((count % DEF_BLOCK_SIZE == 0) ? 0 : 1);

	hipLaunchKernelGGL(k_vector_sum, dim3(blocks), dim3(DEF_BLOCK_SIZE ), 0, 0,  res, a, x, b, y, count );
	hipDeviceSynchronize();
	cudaCheckError();
}

void Cuda_CG_Preconditioner_Fix(real *res, real *a, real *b, int count)
{

	int blocks;

	blocks = (count / DEF_BLOCK_SIZE) + ((count % DEF_BLOCK_SIZE == 0) ? 0 : 1);


	hipLaunchKernelGGL(k_vector_mul, dim3(blocks), dim3(DEF_BLOCK_SIZE ), 0, 0,  res, a, b, count );
	hipDeviceSynchronize( );
	cudaCheckError( );

}

void  Cuda_Copy_Vector_To_Device(real *host_vector, real *device_vector, int nn)
{
	copy_host_device( host_vector, device_vector, sizeof(real) * nn,
			hipMemcpyHostToDevice, "Cuda_CG::b:get" );

}


void  Cuda_Copy_Vector_From_Device(real *host_vector, real *device_vector, int nn)
{
	copy_host_device( host_vector, device_vector, sizeof(real) * nn,
			hipMemcpyDeviceToHost, "Cuda_CG::b:get" );
}


int  compute_nearest_pow_2_fix( int blocks)
{

	int result = 0;
	result = (int) EXP2( CEIL( LOG2((double) blocks)));
	return result;
}





float  Cuda_Calculate_Local_S_Sum(int nn,fix_qeq_gpu *qeq_gpu)
{
	int blocks;
	real *output;
	//cuda malloc this
	cuda_malloc((void **) &output, sizeof(real), TRUE,
			"Cuda_Allocate_Matrix::start");
	double my_acc;


	blocks = nn / DEF_BLOCK_SIZE
			+ (( nn % DEF_BLOCK_SIZE == 0 ) ? 0 : 1);


	Cuda_Reduction_Sum(qeq_gpu->s, output, nn);

	my_acc = 0;

	copy_host_device( &my_acc, output,
			sizeof(real), hipMemcpyDeviceToHost, "charges:x" );

	return my_acc;
}

float  Cuda_Calculate_Local_T_Sum(int nn,fix_qeq_gpu *qeq_gpu)
{
	int blocks;
	real *output;
	//cuda malloc this
	cuda_malloc((void **) &output, sizeof(real), TRUE,
			"Cuda_Allocate_Matrix::start");
	double my_acc;


	blocks = nn / DEF_BLOCK_SIZE
			+ (( nn % DEF_BLOCK_SIZE == 0 ) ? 0 : 1);

	Cuda_Reduction_Sum(qeq_gpu->t, output, nn);



	my_acc = 0;

	copy_host_device( &my_acc, output,
			sizeof(real), hipMemcpyDeviceToHost, "charges:x" );



	return my_acc;

}

CUDA_GLOBAL void k_update_q_and_backup_st(double *q, double *s, double *t, reax_atom *my_atoms,rvec4 *s_hist, rvec4 *t_hist, double u, int nn,reax_atom *sys_my_atoms)
{

	int i;
	int type_i;


	i = blockIdx.x * blockDim.x + threadIdx.x;

	reax_atom *atom;
	atom = &my_atoms[i];

	reax_atom *atom2;
	atom2 = &sys_my_atoms[i];


	if ( i >= nn)
	{
		return;
	}

	q[i]  = atom->q  = sys_my_atoms[i].q = s[i] - u*t[i];

	//printf("S[%d] %f,T[%d] %f,Q[%d] %f \n",i, s[i],i,t[i],i,sys_my_atoms[i].q);


	s_hist[i][3] = s_hist[i][2];
	s_hist[i][2] = s_hist[i][1];
	s_hist[i][1] = s_hist[i][0];
	s_hist[i][0] = s[i];


	t_hist[i][3] = t_hist[i][2];
	t_hist[i][2] = t_hist[i][1];
	t_hist[i][1] = t_hist[i][0];
	t_hist[i][0] = t[i];
}

void  Cuda_Update_Q_And_Backup_ST(int nn, fix_qeq_gpu *qeq_gpu, double u,reax_system *system)
{
	int blocks;
	blocks = nn / DEF_BLOCK_SIZE
			+ (( nn % DEF_BLOCK_SIZE == 0 ) ? 0 : 1);

	hipLaunchKernelGGL(k_update_q_and_backup_st, dim3(blocks), dim3(DEF_BLOCK_SIZE), 0, 0,   qeq_gpu->q,qeq_gpu->s,qeq_gpu->t,qeq_gpu->d_fix_my_atoms,qeq_gpu->s_hist,qeq_gpu->t_hist,u, nn, system->d_my_atoms);
	hipDeviceSynchronize();
	cudaCheckError();

	copy_host_device(qeq_gpu->fix_my_atoms, qeq_gpu->d_fix_my_atoms, sizeof(reax_atom) * system->N,
			hipMemcpyDeviceToHost, "Sync_Atoms::system->my_atoms");
}

void  CudaFreeFixQeqParams(fix_qeq_gpu *qeq_gpu)
{
	   cuda_free(qeq_gpu->s, "S");
       cuda_free(qeq_gpu->t, "T");
       cuda_free(qeq_gpu->Hdia_inv, "Hdia");
       cuda_free(qeq_gpu->b_s, "B_S");
       cuda_free(qeq_gpu->b_t, "B_T");
       cuda_free(qeq_gpu->b_prc, "PRC");
       cuda_free(qeq_gpu->b_prm, "PRM");

       cuda_free(qeq_gpu->p, "P");
       cuda_free(qeq_gpu->q, "Q");
       cuda_free(qeq_gpu->r, "R");
       cuda_free(qeq_gpu->d, "D");
}

void CudaFreeHMatrix(fix_qeq_gpu *qeq_gpu)
{
	   cuda_free(&qeq_gpu->H.start, "start");
       cuda_free(&qeq_gpu->H.end, "end");
       cuda_free(&qeq_gpu->H.entries, "entries");
}


}
