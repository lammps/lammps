#include "hip/hip_runtime.h"
/*----------------------------------------------------------------------
  PuReMD - Purdue ReaxFF Molecular Dynamics Program

  Copyright (2010) Purdue University
  Hasan Metin Aktulga, haktulga@cs.purdue.edu
  Joseph Fogarty, jcfogart@mail.usf.edu
  Sagar Pandit, pandit@usf.edu
  Ananth Y Grama, ayg@cs.purdue.edu

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of 
  the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
  See the GNU General Public License for more details:
  <http://www.gnu.org/licenses/>.
  ----------------------------------------------------------------------*/

#include "cuda_nonbonded.h"

#include "cuda_list.h"
#include "cuda_utils.h"
#include "cuda_reduction.h"
#include "cuda_shuffle.h"

#include "index_utils.h"
#include "vector.h"


//CUDA_GLOBAL void __launch_bounds__ (960) k_vdW_coulomb_energy(    
CUDA_GLOBAL void k_vdW_coulomb_energy( reax_atom *my_atoms, 
		two_body_parameters *tbp, global_parameters gp, control_params *control,
		storage p_workspace, reax_list p_far_nbrs, int n, int N, int num_atom_types,
		real *data_e_vdW, real *data_e_ele, rvec *data_ext_press )
{
#if defined(__SM_35__)
	real sh_vdw;
	real sh_ele;
	rvec sh_force;
#else
	HIP_DYNAMIC_SHARED( real, _vdw)
	HIP_DYNAMIC_SHARED( real, _ele)
	HIP_DYNAMIC_SHARED( rvec, _force)
	real *sh_vdw;
	real *sh_ele;
	rvec *sh_force;
#endif
	int i, j, pj, natoms;
	int start_i, end_i, orig_i, orig_j;
	real p_vdW1, p_vdW1i;
	real powr_vdW1, powgi_vdW1;
	real tmp, r_ij, fn13, exp1, exp2;
	real Tap, dTap, dfn13, CEvd, CEclmb, de_core;
	real dr3gamij_1, dr3gamij_3;
	real e_ele, e_vdW, e_core;
	rvec temp, ext_press;
	two_body_parameters *twbp;
	far_neighbor_data *nbr_pj;
	reax_list *far_nbrs;
	storage *workspace = &p_workspace;
	int thread_id;
	int warpid;
	int laneid;

	thread_id = blockIdx.x * blockDim.x + threadIdx.x;
	warpid = thread_id / VDW_KER_THREADS_PER_ATOM;
	laneid = thread_id & (VDW_KER_THREADS_PER_ATOM -1);
#if defined(__SM_35__)
	sh_vdw = 0.0;
	sh_ele = 0.0;
	rvec_MakeZero ( sh_force );
#else
	sh_vdw = _vdw;
	sh_ele = _vdw + blockDim.x;
	sh_force = (rvec *)( _vdw + 2*blockDim.x);

	sh_vdw[threadIdx.x] = 0.0;
	sh_ele[threadIdx.x] = 0.0;
	rvec_MakeZero ( sh_force [threadIdx.x] );
#endif
	//i = blockIdx.x * blockDim.x + threadIdx.x;
	//if (i >= N) return;
	i = warpid;

	if ( i < N )
	{
		natoms = n;
		far_nbrs = &p_far_nbrs;
		p_vdW1 = gp.l[28];
		p_vdW1i = 1.0 / p_vdW1;
		e_core = 0;
		e_vdW = 0;

		data_e_vdW[i] = 0;
		data_e_ele[i] = 0;

		//for( i = 0; i < natoms; ++i ) {
		start_i = Cuda_Start_Index(i, far_nbrs);
		end_i = Cuda_End_Index(i, far_nbrs);
		orig_i = my_atoms[i].orig_id;
		//fprintf( stderr, "i:%d, start_i: %d, end_i: %d\n", i, start_i, end_i );

		//for( pj = start_i; pj < end_i; ++pj )
		pj = start_i + laneid;
		while ( pj < end_i )
		{

			nbr_pj = &far_nbrs->select.far_nbr_list[pj];
			j = nbr_pj->nbr;
			orig_j  = my_atoms[j].orig_id;

			if( nbr_pj->d <= control->nonb_cut &&
					(((i < j) && (i < natoms) && (j < natoms || orig_i < orig_j))
							|| ((i > j) && (i < natoms) && (j < natoms))
							|| (i > j && i >= natoms && j < natoms && orig_j < orig_i)))
			{ // ji with j >= n
				r_ij = nbr_pj->d;
				twbp = &tbp[ index_tbp(my_atoms[i].type, my_atoms[j].type, num_atom_types) ];

				/* Calculate Taper and its derivative */
				// Tap = nbr_pj->Tap;   -- precomputed during compte_H
				Tap = workspace->Tap[7] * r_ij + workspace->Tap[6];
				Tap = Tap * r_ij + workspace->Tap[5];
				Tap = Tap * r_ij + workspace->Tap[4];
				Tap = Tap * r_ij + workspace->Tap[3];
				Tap = Tap * r_ij + workspace->Tap[2];
				Tap = Tap * r_ij + workspace->Tap[1];
				Tap = Tap * r_ij + workspace->Tap[0];

				dTap = 7 * workspace->Tap[7] * r_ij + 6 * workspace->Tap[6];
				dTap = dTap * r_ij + 5*workspace->Tap[5];
				dTap = dTap * r_ij + 4*workspace->Tap[4];
				dTap = dTap * r_ij + 3*workspace->Tap[3];
				dTap = dTap * r_ij + 2*workspace->Tap[2];
				dTap += workspace->Tap[1] / r_ij;

				/* shielding vdWaals Calculations */
				if ( gp.vdw_type == 1 || gp.vdw_type == 3 )
				{
					powr_vdW1 = POW( r_ij, p_vdW1 );
					powgi_vdW1 = POW( 1.0 / twbp->gamma_w, p_vdW1 );

					fn13 = POW( powr_vdW1 + powgi_vdW1, p_vdW1i );
					exp1 = EXP( twbp->alpha * (1.0 - fn13 / twbp->r_vdW) );
					exp2 = EXP( 0.5 * twbp->alpha * (1.0 - fn13 / twbp->r_vdW) );

					e_vdW = twbp->D * (exp1 - 2.0 * exp2);

					//data_e_vdW[i] += Tap * e_vdW;
					//data_e_vdW[i] += Tap * e_vdW / 2.0;
#if defined(__SM_35__)
					sh_vdw += Tap * e_vdW / 2.0;
#else
					sh_vdw[threadIdx.x] += Tap * e_vdW / 2.0;
#endif

					dfn13 = POW( powr_vdW1 + powgi_vdW1, p_vdW1i - 1.0) *
							POW( r_ij, p_vdW1 - 2.0 );

					CEvd = dTap * e_vdW -
							Tap * twbp->D * (twbp->alpha / twbp->r_vdW) * (exp1 - exp2) * dfn13;
				}
				/* no shielding */
				else
				{
					exp1 = EXP( twbp->alpha * (1.0 - r_ij / twbp->r_vdW) );
					exp2 = EXP( 0.5 * twbp->alpha * (1.0 - r_ij / twbp->r_vdW) );

					e_vdW = twbp->D * (exp1 - 2.0 * exp2);

					//data_e_vdW[i] += Tap * e_vdW;
					//data_e_vdW[i] += Tap * e_vdW / 2.0;
#if defined(__SM_35__)
					sh_vdw += Tap * e_vdW / 2.0;
#else
					sh_vdw[threadIdx.x] += Tap * e_vdW / 2.0;
#endif

					CEvd = dTap * e_vdW -
							Tap * twbp->D * (twbp->alpha / twbp->r_vdW) * (exp1 - exp2);
				}

				/* inner wall */
				if ( gp.vdw_type == 2 || gp.vdw_type == 3 )
				{
					e_core = twbp->ecore * EXP(twbp->acore * (1.0-(r_ij/twbp->rcore)));

					//data_e_vdW[i] += Tap * e_core;
					//data_e_vdW[i] += Tap * e_core / 2.0;
#if defined(__SM_35__)
					sh_vdw += Tap * e_core / 2.0;
#else
					sh_vdw[ threadIdx.x ] += Tap * e_core / 2.0;
#endif

					de_core = -(twbp->acore/twbp->rcore) * e_core;
					CEvd += dTap * e_core + Tap * de_core;
				}

				/*Coulomb Calculations*/
				dr3gamij_1 = r_ij * r_ij * r_ij + twbp->gamma;
				dr3gamij_3 = POW( dr3gamij_1, 1.0 / 3.0 );

				tmp = Tap / dr3gamij_3;
				//data_e_ele[i] += e_ele = C_ele * my_atoms[i].q * my_atoms[j].q * tmp;
				e_ele = C_ele * my_atoms[i].q * my_atoms[j].q * tmp;


				//if(my_atoms[thread_id].orig_id == 8)
					//printf("%f,%f,%f,%f,%d,%d\n", tmp,C_ele,my_atoms[i].q, my_atoms[j].q,i,j);
				   //   printf("%d,%d,%f,%f,%f\n",thread_id , my_atoms[thread_id].orig_id,workspace->f[thread_id][0], workspace->f[thread_id][1], workspace->f[thread_id][2]);




				//data_e_ele[i] += e_ele;
				//data_e_ele[i] += e_ele  / 2.0;
#if defined(__SM_35__)
				sh_ele += e_ele  / 2.0;
#else
				sh_ele[ threadIdx.x ] += e_ele  / 2.0;
#endif

				CEclmb = C_ele * my_atoms[i].q * my_atoms[j].q *
						( dTap -  Tap * r_ij / dr3gamij_1 ) / dr3gamij_3;

				// fprintf( fout, "%5d %5d %10.6f %10.6f\n",
				//   MIN( system->my_atoms[i].orig_id, system->my_atoms[j].orig_id ),
				//   MAX( system->my_atoms[i].orig_id, system->my_atoms[j].orig_id ),
				//   CEvd, CEclmb );

				if ( control->virial == 0 )
				{
					if ( i < j )
					{
						//rvec_ScaledAdd( workspace->f[i], -(CEvd + CEclmb), nbr_pj->dvec );
#if defined (__SM_35__)
						rvec_ScaledAdd( sh_force, -(CEvd + CEclmb), nbr_pj->dvec );
#else
						rvec_ScaledAdd( sh_force[ threadIdx.x ], -(CEvd + CEclmb), nbr_pj->dvec );
#endif
					}
					else
					{
						//rvec_ScaledAdd( workspace->f[i], +(CEvd + CEclmb), nbr_pj->dvec );
#if defined (__SM_35__)
						rvec_ScaledAdd( sh_force , +(CEvd + CEclmb), nbr_pj->dvec );
#else
						rvec_ScaledAdd( sh_force[ threadIdx.x ], +(CEvd + CEclmb), nbr_pj->dvec );
#endif
						//rvec_ScaledAdd( workspace->f[j], +(CEvd + CEclmb), nbr_pj->dvec );
					}
				}
				/* NPT, iNPT or sNPT */
				else
				{
					/* for pressure coupling, terms not related to bond order
                       derivatives are added directly into pressure vector/tensor */
					rvec_Scale( temp, CEvd + CEclmb, nbr_pj->dvec );

					rvec_ScaledAdd( workspace->f[i], -1., temp );
					rvec_Add( workspace->f[j], temp );

					rvec_iMultiply( ext_press, nbr_pj->rel_box, temp );
					rvec_Add( data_ext_press [i], ext_press );

					// fprintf( stderr, "nonbonded(%d,%d): rel_box (%f %f %f)
					//   force(%f %f %f) ext_press (%12.6f %12.6f %12.6f)\n",
					//   i, j, nbr_pj->rel_box[0], nbr_pj->rel_box[1], nbr_pj->rel_box[2],
					//   temp[0], temp[1], temp[2],
					//   data->ext_press[0], data->ext_press[1], data->ext_press[2] );
				}

#ifdef TEST_ENERGY
				// fprintf( out_control->evdw,
				// "%12.9f%12.9f%12.9f%12.9f%12.9f%12.9f%12.9f%12.9f\n",
				// workspace->Tap[7],workspace->Tap[6],workspace->Tap[5],
				// workspace->Tap[4],workspace->Tap[3],workspace->Tap[2],
				// workspace->Tap[1], Tap );
				//fprintf( out_control->evdw, "%6d%6d%24.15e%24.15e%24.15e\n",
				fprintf( out_control->evdw, "%6d%6d%12.4f%12.4f%12.4f\n",
						system->my_atoms[i].orig_id, system->my_atoms[j].orig_id,
						r_ij, e_vdW, data->my_en.e_vdW );
				//fprintf(out_control->ecou,"%6d%6d%24.15e%24.15e%24.15e%24.15e%24.15e\n",
				fprintf( out_control->ecou, "%6d%6d%12.4f%12.4f%12.4f%12.4f%12.4f\n",
						system->my_atoms[i].orig_id, system->my_atoms[j].orig_id,
						r_ij, system->my_atoms[i].q, system->my_atoms[j].q,
						e_ele, data->my_en.e_ele );
#endif

#ifdef TEST_FORCES
				rvec_ScaledAdd( workspace->f_vdw[i], -CEvd, nbr_pj->dvec );
				rvec_ScaledAdd( workspace->f_vdw[j], +CEvd, nbr_pj->dvec );
				rvec_ScaledAdd( workspace->f_ele[i], -CEclmb, nbr_pj->dvec );
				rvec_ScaledAdd( workspace->f_ele[j], +CEclmb, nbr_pj->dvec );
#endif
			}

			pj += VDW_KER_THREADS_PER_ATOM;

		}
		//  }
	} // if i < N

#if defined( __SM_35__)
	for ( int x = VDW_KER_THREADS_PER_ATOM >> 1; x >= 1; x/=2 )
	{
		sh_vdw += shfl( sh_vdw, x);
		sh_ele += shfl( sh_ele, x );
		sh_force[0] += shfl( sh_force[0], x );
		sh_force[1] += shfl( sh_force[1], x );
		sh_force[2] += shfl( sh_force[2], x );
	}

	if ( laneid == 0 )
	{
		data_e_vdW[i] += sh_vdw;
		data_e_ele[i] += sh_ele;
		rvec_Add( workspace->f[i], sh_force );
	}

#else

	__syncthreads( );

	if (laneid < 16)
	{
		sh_vdw[threadIdx.x] += sh_vdw[threadIdx.x + 16];
		sh_ele[threadIdx.x] += sh_ele[threadIdx.x + 16];
		rvec_Add( sh_force[threadIdx.x], sh_force[threadIdx.x + 16] );
	}
	__syncthreads( );
	if (laneid < 8)
	{
		sh_vdw[threadIdx.x] += sh_vdw[threadIdx.x + 8];
		sh_ele[threadIdx.x] += sh_ele[threadIdx.x + 8];
		rvec_Add( sh_force[threadIdx.x], sh_force[threadIdx.x + 8] );
	}
	__syncthreads( );
	if (laneid < 4)
	{
		sh_vdw[threadIdx.x] += sh_vdw[threadIdx.x + 4];
		sh_ele[threadIdx.x] += sh_ele[threadIdx.x + 4];
		rvec_Add( sh_force[threadIdx.x], sh_force[threadIdx.x + 4] );
	}
	__syncthreads( );
	if (laneid < 2)
	{
		sh_vdw[threadIdx.x] += sh_vdw[threadIdx.x + 2];
		sh_ele[threadIdx.x] += sh_ele[threadIdx.x + 2];
		rvec_Add( sh_force[threadIdx.x], sh_force[threadIdx.x + 2] );
	}
	__syncthreads( );
	if (laneid < 1)
	{
		sh_vdw[threadIdx.x] += sh_vdw[threadIdx.x + 1];
		sh_ele[threadIdx.x] += sh_ele[threadIdx.x + 1];
		rvec_Add( sh_force[threadIdx.x], sh_force[threadIdx.x + 1] );
	}
	__syncthreads( );
	if (laneid == 0)
	{
		data_e_vdW[i] += sh_vdw[threadIdx.x];
		data_e_ele[i] += sh_ele[threadIdx.x];
		rvec_Add( workspace->f[i], sh_force[ threadIdx.x ] );
	}
#endif

	// if(thread_id < 20)
	  //    printf("%d,%d,%f,%f,%f\n",thread_id , my_atoms[thread_id].orig_id,workspace->f[thread_id][0], workspace->f[thread_id][1], workspace->f[thread_id][2]);


}


CUDA_GLOBAL void k_tabulated_vdW_coulomb_energy( reax_atom *my_atoms, 
		global_parameters gp, control_params *control,
		storage p_workspace, reax_list p_far_nbrs,
		LR_lookup_table *t_LR, int n, int N, int num_atom_types,
		int step, int prev_steps, int energy_update_freq,
		real *data_e_vdW, real *data_e_ele, rvec *data_ext_press )
{
	int i, j, pj, r, natoms, steps, update_freq, update_energies;
	int type_i, type_j, tmin, tmax;
	int start_i, end_i, orig_i, orig_j;
	real r_ij, base, dif;
	real e_vdW, e_ele;
	real CEvd, CEclmb;
	rvec temp, ext_press;
	far_neighbor_data *nbr_pj;
	reax_list *far_nbrs;
	LR_lookup_table *t;
	storage *workspace;

	i = blockIdx.x * blockDim.x + threadIdx.x;

	if ( i >= N )
	{
		return;
	}

	workspace = &p_workspace;
	natoms = n;
	far_nbrs = &p_far_nbrs;
	steps = step - prev_steps;
	update_freq = energy_update_freq;
	update_energies = update_freq > 0 && steps % update_freq == 0;
	e_ele = e_vdW = 0;
	data_e_vdW[i] = 0;
	data_e_ele[i] = 0;

	//for( i = 0; i < natoms; ++i ) {
	type_i = my_atoms[i].type;
	start_i = Cuda_Start_Index(i,far_nbrs);
	end_i = Cuda_End_Index(i,far_nbrs);
	orig_i = my_atoms[i].orig_id;

	for ( pj = start_i; pj < end_i; ++pj )
	{
		nbr_pj = &far_nbrs->select.far_nbr_list[pj];
		j = nbr_pj->nbr;
		orig_j  = my_atoms[j].orig_id;

		//if ( nbr_pj->d <= control->nonb_cut && (j < natoms || orig_i < orig_j) ) {
		if ( nbr_pj->d <= control->nonb_cut &&
				(((i < j) && (i < natoms) && (j < natoms || orig_i < orig_j))
						|| ((i > j) && (i < natoms) && (j < natoms))
						|| (i > j && i >= natoms && j < natoms && orig_j < orig_i)))
		{ // ji with j >= n
			j = nbr_pj->nbr;
			type_j = my_atoms[j].type;
			r_ij   = nbr_pj->d;
			tmin  = MIN( type_i, type_j );
			tmax  = MAX( type_i, type_j );

			t = &t_LR[ index_lr(tmin, tmax, num_atom_types) ];

			//table = &LR[type_i][type_j];

			/* Cubic Spline Interpolation */
			r = (int)(r_ij * t->inv_dx);
			if( r == 0 )
			{
				++r;
			}
			base = (real)(r+1) * t->dx;
			dif = r_ij - base;
			//fprintf(stderr, "r: %f, i: %d, base: %f, dif: %f\n", r, i, base, dif);

			if ( update_energies )
			{
				e_vdW = ((t->vdW[r].d*dif + t->vdW[r].c)*dif + t->vdW[r].b)*dif +
						t->vdW[r].a;

				e_ele = ((t->ele[r].d*dif + t->ele[r].c)*dif + t->ele[r].b)*dif +
						t->ele[r].a;
				e_ele *= my_atoms[i].q * my_atoms[j].q;

				//data_e_vdW[i] += e_vdW;
				data_e_vdW[i] += e_vdW / 2.0;
				//data_e_ele[i] += e_ele;
				data_e_ele[i] += e_ele / 2.0;
			}

			CEvd = ((t->CEvd[r].d * dif + t->CEvd[r].c) * dif + t->CEvd[r].b) * dif +
					t->CEvd[r].a;

			CEclmb = ((t->CEclmb[r].d * dif + t->CEclmb[r].c) * dif + t->CEclmb[r].b) * dif +
					t->CEclmb[r].a;
			CEclmb *= my_atoms[i].q * my_atoms[j].q;

			if( control->virial == 0 )
			{
				if ( i < j )
				{
					rvec_ScaledAdd( workspace->f[i], -(CEvd + CEclmb), nbr_pj->dvec );
				}
				else
				{
					rvec_ScaledAdd( workspace->f[i], +(CEvd + CEclmb), nbr_pj->dvec );

				}
				//rvec_ScaledAdd( workspace->f[i], -(CEvd + CEclmb), nbr_pj->dvec );
				//rvec_ScaledAdd( workspace->f[j], +(CEvd + CEclmb), nbr_pj->dvec );
			}
			/* NPT, iNPT or sNPT */
			else
			{
				/* for pressure coupling, terms not related to bond order derivatives
		                   are added directly into pressure vector/tensor */
				rvec_Scale( temp, CEvd + CEclmb, nbr_pj->dvec );

				rvec_ScaledAdd( workspace->f[i], -1., temp );
				rvec_Add( workspace->f[j], temp );

				rvec_iMultiply( ext_press, nbr_pj->rel_box, temp );
				rvec_Add( data_ext_press [i], ext_press );
			}


		}
	}

	/*if(i < 5)
	{
		printf("%f,%f,%f\n", workspace->f[i][0],workspace->f[i][1],workspace->f[i][2]);
	}*/
	//  }

}


CUDA_GLOBAL void k_pol_energy( reax_atom *my_atoms, 
		single_body_parameters *sbp, int n, real *data_e_pol )
{
	int i, type_i;
	real q;

	i = blockIdx.x * blockDim.x + threadIdx.x;

	if ( i >= n )
	{
		return;
	}

	q = my_atoms[i].q;
	type_i = my_atoms[i].type;

	data_e_pol[i] =
			KCALpMOL_to_EV * (sbp[type_i].chi * q +
					(sbp[type_i].eta / 2.) * SQR(q));
}


static void Cuda_Compute_Polarization_Energy( reax_system *system, storage *workspace,
		simulation_data *data )
{
	int blocks;
	real *spad = (real *) workspace->scratch;

	cuda_memset( spad, 0, sizeof(real) * 2 * system->n, "pol_energy" );

	blocks = system->n / DEF_BLOCK_SIZE +
			((system->n % DEF_BLOCK_SIZE == 0) ? 0 : 1);

	hipLaunchKernelGGL(k_pol_energy, dim3(blocks), dim3(DEF_BLOCK_SIZE ), 0, 0,  system->d_my_atoms, system->reax_param.d_sbp,
			system->n, spad );
	cudaCheckError( );

	Cuda_Reduction_Sum( spad,
			&((simulation_data *)data->d_simulation_data)->my_en.e_pol,
			system->n );
}


void Cuda_NonBonded_Energy( reax_system *system, control_params *control, 
		storage *workspace, simulation_data *data, reax_list **lists,
		output_controls *out_control, bool isTabulated )
{
	int blocks, rblocks, update_energy;
	int size = (2 * system->N + 2 * system->N ) * sizeof(real) +
			2 * system->N * sizeof(rvec);
	//printf("Size %d \n", size);
	rvec *spad_rvec;
	real *spad = (real *) workspace->scratch;

	update_energy = (out_control->energy_update_freq > 0
			&& data->step % out_control->energy_update_freq == 0) ? TRUE : FALSE;
	rblocks = system->N / DEF_BLOCK_SIZE + ((system->N % DEF_BLOCK_SIZE == 0) ? 0 : 1);
	blocks = ((system->N * VDW_KER_THREADS_PER_ATOM) / DEF_BLOCK_SIZE)
        												+ (((system->N * VDW_KER_THREADS_PER_ATOM) % DEF_BLOCK_SIZE == 0) ? 0 : 1);

	cuda_memset( spad, 0, size, "pol_energy" );

	if ( !isTabulated )
	{
		hipLaunchKernelGGL(k_vdW_coulomb_energy, dim3(blocks), dim3(DEF_BLOCK_SIZE), DEF_BLOCK_SIZE * (2 * sizeof(real) + sizeof(rvec)) , 0,  system->d_my_atoms, system->reax_param.d_tbp,
				system->reax_param.d_gp, (control_params *)control->d_control_params,
				*(workspace->d_workspace), *(lists[FAR_NBRS]),
				system->n, system->N, system->reax_param.num_atom_types,
				spad, spad + 2 * system->N, (rvec *)(spad + 4 * system->N));
		hipDeviceSynchronize( );
		cudaCheckError( );
	}
	else
	{
		hipLaunchKernelGGL(k_tabulated_vdW_coulomb_energy, dim3(blocks), dim3(DEF_BLOCK_SIZE ), 0, 0,  system->d_my_atoms, system->reax_param.d_gp,
				(control_params *)control->d_control_params,
				*(workspace->d_workspace), *(lists[FAR_NBRS]),
				workspace->d_LR, system->n, system->N,
				system->reax_param.num_atom_types,
				data->step, data->prev_steps,
				out_control->energy_update_freq,
				spad, spad + 2 * system->N,
				(rvec *)(spad + 4 * system->N));
		hipDeviceSynchronize( );
		cudaCheckError();
	}

	/* reduction for vdw */
	if ( update_energy == TRUE )
	{
		Cuda_Reduction_Sum( spad,
				&((simulation_data *)data->d_simulation_data)->my_en.e_vdW,
				system->N );
	}

	/* reduction for ele */
	if ( update_energy == TRUE )
	{
		Cuda_Reduction_Sum( spad + 2 * system->N,
				&((simulation_data *)data->d_simulation_data)->my_en.e_ele,
				system->N );
	}

	/* reduction for ext_press */
	spad_rvec = (rvec *) (spad + 4 * system->N);
	hipLaunchKernelGGL(k_reduction_rvec, dim3(rblocks), dim3(DEF_BLOCK_SIZE), sizeof(rvec) * DEF_BLOCK_SIZE , 0,  spad_rvec, spad_rvec + system->N, system->N);
	hipDeviceSynchronize( );
	cudaCheckError( );

	hipLaunchKernelGGL(k_reduction_rvec, dim3(1), dim3(control->blocks_pow_2_n), sizeof(rvec) * control->blocks_pow_2_n, 0,  spad_rvec + system->N, &((simulation_data *)data->d_simulation_data)->my_ext_press, rblocks);
	hipDeviceSynchronize( );
	cudaCheckError( );

	if ( update_energy == TRUE )
	{
		Cuda_Compute_Polarization_Energy( system, workspace, data );
	}
}
