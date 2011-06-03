/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator 

   Original Version:
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov 

   See the README file in the top-level LAMMPS directory. 

   ----------------------------------------------------------------------- 

   USER-CUDA Package and associated modifications:
   https://sourceforge.net/projects/lammpscuda/ 

   Christian Trott, christian.trott@tu-ilmenau.de
   Lars Winterfeld, lars.winterfeld@tu-ilmenau.de
   Theoretical Physics II, University of Technology Ilmenau, Germany 

   See the README file in the USER-CUDA directory. 

   This software is distributed under the GNU General Public License.
------------------------------------------------------------------------- */

#ifdef CUDA_USE_BINNING
#include <stdio.h>
#define MY_PREFIX binning
#include "cuda_shared.h"
#include "cuda_common.h"
#include "crm_cuda_utils.cu"
#include "binning_cu.h"
#include "binning_kernel.cu"

void Cuda_PreBinning(cuda_shared_data* sdata)
{
	// initialize only on first call
	short init = 0;
	if(! init)
	{
		init = 1;
		int cuda_dummy_type = sdata->atom.ntypes + 1;
		X_FLOAT outside[3] =
		{
			(sdata->domain.subhi[0] - sdata->domain.sublo[0])/1000.0,
			(sdata->domain.subhi[1] - sdata->domain.sublo[1])/1000.0,
			(sdata->domain.subhi[2] - sdata->domain.sublo[2])/1000.0
		};
		cudaMemcpyToSymbol("binned_size_all"    , & sdata->atom.binned_type.dim[0]  , sizeof(unsigned) );
		cudaMemcpyToSymbol("cuda_dummy_type"    , & cuda_dummy_type                 , sizeof(int)      );
		cudaMemcpyToSymbol("outside"            , & outside                         , sizeof(X_FLOAT)*3);
		cudaMemcpyToSymbol(MY_CONST(binned_type), & sdata->atom.binned_type.dev_data, sizeof(int*)     );
		cudaMemcpyToSymbol(MY_CONST(binned_x)   , & sdata->atom.binned_x   .dev_data, sizeof(X_FLOAT*) );
		cudaMemcpyToSymbol(MY_CONST(binned_tag) , & sdata->atom.binned_tag .dev_data, sizeof(int*)     );
		cudaMemcpyToSymbol(MY_CONST(subhi)      ,   sdata->domain.subhi             , sizeof(X_FLOAT)*3);
		// bin_nmax == blockDim.x
		
		// printf("# CUDA: MY_CONST(binned_type) = %s\n", MY_CONST(binned_type));
		// int* p = pre_binning_binned_type; // pre_binning_binned_type is defined here!!
	}
	
	dim3 grid(sdata->domain.bin_dim[0], sdata->domain.bin_dim[1] * sdata->domain.bin_dim[2], 1);
	dim3 threads(sdata->domain.bin_nmax, 1, 1);
	
	MYDBG(printf("# CUDA: Cuda_PreBinning: pre binning grid = (%u, %u, %u)\n", grid.x, grid.y, grid.z);)
	MYDBG(printf("# CUDA: Cuda_PreBinning: pre binning threads = (%u, %u, %u)\n", threads.x, threads.y, threads.z);	)
	PreBinning_Kernel<<<grid, threads>>> ();
	cudaThreadSynchronize();
    MYDBG(printf("ERROR-CUDA pre_binning: %s\n",cudaGetErrorString(cudaGetLastError())));
	CUT_CHECK_ERROR("Cuda_PreBinning: binning Kernel execution failed");
}

void Cuda_Binning(cuda_shared_data* sdata)
{
	MYDBG(	// check assumption in debug mode
		if(sdata->atom.x.dim[1] != 3)
		{
			printf("# CUDA: Cuda_Binning: binning error: atom array dimensions not Nx3\n");
			return;
		}
	)
	
	// initialize only on first call
	short init = 0;
	if(! init)
	{
		init = 0;
		X_FLOAT const_rez_bin_size[3] = 
		{
			(1.0 * sdata->domain.bin_dim[0]-4.0) / (sdata->domain.subhi[0] - sdata->domain.sublo[0]),
			(1.0 * sdata->domain.bin_dim[1]-4.0) / (sdata->domain.subhi[1] - sdata->domain.sublo[1]),
			(1.0 * sdata->domain.bin_dim[2]-4.0) / (sdata->domain.subhi[2] - sdata->domain.sublo[2])
		};
		cudaMemcpyToSymbol("bin_error_count"        , & sdata->atom.bin_error_count.dev_data, sizeof(X_FLOAT)*1);
		cudaMemcpyToSymbol("rez_bin_size"           , & const_rez_bin_size                  , sizeof(X_FLOAT)*3);
		cudaMemcpyToSymbol(MY_CONST(bin_count_all)  , & sdata->atom.bin_count_all  .dev_data, sizeof(unsigned*));
		cudaMemcpyToSymbol(MY_CONST(bin_count_local), & sdata->atom.bin_count_local.dev_data, sizeof(unsigned*));
		cudaMemcpyToSymbol(MY_CONST(bin_dim)        ,   sdata->domain.bin_dim               , sizeof(int3)     );
		cudaMemcpyToSymbol(MY_CONST(bin_nmax)       , & sdata->domain.bin_nmax              , sizeof(unsigned) );
		cudaMemcpyToSymbol(MY_CONST(binned_f)       , & sdata->atom.binned_f       .dev_data, sizeof(F_FLOAT*) );
		cudaMemcpyToSymbol(MY_CONST(binned_q)       , & sdata->atom.binned_q       .dev_data, sizeof(F_FLOAT*) );
		cudaMemcpyToSymbol(MY_CONST(binned_rmass)   , & sdata->atom.binned_rmass   .dev_data, sizeof(V_FLOAT*) );
		cudaMemcpyToSymbol(MY_CONST(binned_tag)     , & sdata->atom.binned_tag     .dev_data, sizeof(int*)     );
		cudaMemcpyToSymbol(MY_CONST(binned_type)    , & sdata->atom.binned_type    .dev_data, sizeof(int*)     );
		cudaMemcpyToSymbol(MY_CONST(binned_v)       , & sdata->atom.binned_v       .dev_data, sizeof(V_FLOAT*) );
		cudaMemcpyToSymbol(MY_CONST(binpos)         , & sdata->atom.binpos         .dev_data, sizeof(int*));
		cudaMemcpyToSymbol(MY_CONST(f)              , & sdata->atom.f              .dev_data, sizeof(F_FLOAT*) );
		cudaMemcpyToSymbol(MY_CONST(natoms)         , & sdata->atom.nall                    , sizeof(unsigned) );
		cudaMemcpyToSymbol(MY_CONST(nghost)         , & sdata->atom.nghost                  , sizeof(unsigned) );
		cudaMemcpyToSymbol(MY_CONST(nlocal)         , & sdata->atom.nlocal                  , sizeof(unsigned) );
		cudaMemcpyToSymbol(MY_CONST(nmax)           , & sdata->atom.nmax                    , sizeof(unsigned) );
		cudaMemcpyToSymbol(MY_CONST(q)              , & sdata->atom.q              .dev_data, sizeof(F_FLOAT*) );
		cudaMemcpyToSymbol(MY_CONST(rmass)          , & sdata->atom.rmass          .dev_data, sizeof(V_FLOAT*) );
		cudaMemcpyToSymbol(MY_CONST(sublo)          ,   sdata->domain.sublo                 , sizeof(X_FLOAT)*3);
		cudaMemcpyToSymbol(MY_CONST(tag)            , & sdata->atom.tag            .dev_data, sizeof(int*)     );
		cudaMemcpyToSymbol(MY_CONST(type)           , & sdata->atom.type           .dev_data, sizeof(int*)     );
		cudaMemcpyToSymbol(MY_CONST(v)              , & sdata->atom.v              .dev_data, sizeof(V_FLOAT*) );
	}
	
	dim3 grid((unsigned)(1 + sdata->atom.nlocal/64.0), 1, 1);
	MYDBG( printf("# CUDA: Cuda_Binning: grid dim.x = %u (nlocal: %i)\n", grid.x,sdata->atom.nlocal); )
	dim3 threads(64, 1, 1);
	
	cudaMemset((int*) (sdata->atom.bin_count_all.dev_data),0,sizeof(int)*(sdata->domain.bin_dim[0])*(sdata->domain.bin_dim[1])*(sdata->domain.bin_dim[2]));
	cudaMemset((int*) (sdata->atom.bin_count_local.dev_data),0,sizeof(int)*(sdata->domain.bin_dim[0])*(sdata->domain.bin_dim[1])*(sdata->domain.bin_dim[2]));
	cudaMemset(sdata->atom.bin_error_count.dev_data,0,sizeof(int)*1);
	int binning_error_l[1];
	
	
	Binning_Kernel<<<grid, threads>>> (
		(X_FLOAT*) (sdata->atom.       x.dev_data),
		(X_FLOAT*) (sdata->atom.binned_x.dev_data),
		sdata->atom.q_flag,
		0,
		sdata->atom.rmass_flag
	);
	cudaThreadSynchronize();
	cudaMemcpy((void*) binning_error_l,(void*) sdata->atom.bin_error_count.dev_data,1*sizeof(int),cudaMemcpyDeviceToHost);
	if(binning_error_l[0]!=0) 
	{
		printf("CUDA-ERROR: binning local: could not bin %i atoms\n",binning_error_l[0]);
	}
	CUT_CHECK_ERROR("Cuda_Binning: binning Kernel execution failed");
	
	grid.x=(unsigned)(1 + (sdata->atom.nall-sdata->atom.nlocal)/32.0);
	MYDBG( printf("# CUDA: Cuda_Binning Ghost: grid dim.x = %u\n", grid.x); )
	
	
	Binning_Kernel<<<grid, threads>>> (
		(X_FLOAT*) (sdata->atom.       x.dev_data),
		(X_FLOAT*) (sdata->atom.binned_x.dev_data),
		sdata->atom.q_flag,
		sdata->atom.nlocal,
		sdata->atom.rmass_flag
	);
	cudaThreadSynchronize();
	cudaMemcpy((void*) binning_error_l,(void*) sdata->atom.bin_error_count.dev_data,1*sizeof(int),cudaMemcpyDeviceToHost);
	if(binning_error_l[0]!=0) printf("CUDA-ERROR: binning ghost: could not bin %i atoms\n",binning_error_l[0]);
}

void Cuda_ReverseBinning(cuda_shared_data* sdata)
{
	// initialize only on first call
	short init = 0;
	if(! init)
	{
		init = 0;
		cudaMemcpyToSymbol(MY_CONST(bin_count_all)  , & sdata->atom.bin_count_all  .dev_data, sizeof(unsigned*));
		cudaMemcpyToSymbol(MY_CONST(bin_count_local), & sdata->atom.bin_count_local.dev_data, sizeof(unsigned*));
		cudaMemcpyToSymbol(MY_CONST(bin_dim)        ,   sdata->domain.bin_dim               , sizeof(int3)     );
		cudaMemcpyToSymbol(MY_CONST(binned_f)       , & sdata->atom.binned_f       .dev_data, sizeof(F_FLOAT*) );
		cudaMemcpyToSymbol(MY_CONST(binned_q)       , & sdata->atom.binned_q       .dev_data, sizeof(F_FLOAT*) );
		cudaMemcpyToSymbol(MY_CONST(binned_tag)     , & sdata->atom.binned_tag     .dev_data, sizeof(int*)     );
		cudaMemcpyToSymbol(MY_CONST(binned_type)    , & sdata->atom.binned_type    .dev_data, sizeof(int*)     );
		cudaMemcpyToSymbol(MY_CONST(binned_v)       , & sdata->atom.binned_v       .dev_data, sizeof(V_FLOAT*) );
		cudaMemcpyToSymbol(MY_CONST(f)              , & sdata->atom.f              .dev_data, sizeof(F_FLOAT*) );
		cudaMemcpyToSymbol(MY_CONST(natoms)         , & sdata->atom.nall                    , sizeof(unsigned) );
		cudaMemcpyToSymbol(MY_CONST(nmax)           , & sdata->atom.nmax                    , sizeof(unsigned) );
		cudaMemcpyToSymbol(MY_CONST(q)              , & sdata->atom.q              .dev_data, sizeof(F_FLOAT*) );
		cudaMemcpyToSymbol(MY_CONST(tag)            , & sdata->atom.tag            .dev_data, sizeof(int*)     );
		cudaMemcpyToSymbol(MY_CONST(type)           , & sdata->atom.type           .dev_data, sizeof(int*)     );
		cudaMemcpyToSymbol(MY_CONST(v)              , & sdata->atom.v              .dev_data, sizeof(V_FLOAT*) );
	}
	
	dim3 grid((unsigned)(1 + sdata->atom.nlocal/32.0), 1, 1);
	MYDBG( printf("# CUDA: Cuda_ReverseBinning: grid dim.x = %u (nlocal: %i)\n", grid.x,sdata->atom.nlocal); )
	dim3 threads(32, 1, 1);

	ReverseBinning_Kernel<<<grid, threads>>> (
		(X_FLOAT*) (sdata->atom.       x.dev_data),
		(X_FLOAT*) (sdata->atom.binned_x.dev_data),
		sdata->atom.q_flag
	);
	cudaThreadSynchronize();
	CUT_CHECK_ERROR("Cuda_Binning: reverse binning Kernel execution failed");
}

#endif
