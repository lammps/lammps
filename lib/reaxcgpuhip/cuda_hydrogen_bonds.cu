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

#include "cuda_hydrogen_bonds.h"

#include "cuda_valence_angles.h"
#include "cuda_helpers.h"
#include "cuda_list.h"
#include "cuda_shuffle.h"

#include "index_utils.h"
#include "vector.h"


CUDA_GLOBAL void Cuda_Hydrogen_Bonds( reax_atom *my_atoms, single_body_parameters *sbp, 
		hbond_parameters *d_hbp, global_parameters gp, control_params *control,
		storage p_workspace, reax_list p_bonds, reax_list p_hbonds, int n,
		int num_atom_types, real *data_e_hb, rvec *data_ext_press, int rank, int step )
{
	int i, j, k, pi, pk;
	int type_i, type_j, type_k;
	int start_j, end_j, hb_start_j, hb_end_j;
	int hblist[MAX_BONDS];
	int itr, top;
	ivec rel_jk;
	real r_jk, theta, cos_theta, sin_xhz4, cos_xhz1, sin_theta2;
	real e_hb, exp_hb2, exp_hb3, CEhb1, CEhb2, CEhb3;
	rvec dcos_theta_di, dcos_theta_dj, dcos_theta_dk;
	rvec dvec_jk, force, ext_press;
	hbond_parameters *hbp;
	bond_order_data *bo_ij;
	bond_data *pbond_ij;
	far_neighbor_data *nbr_jk;
	reax_list *bonds, *hbonds;
	bond_data *bond_list;
	hbond_data *hbond_list, *hbond_jk;
	storage *workspace;

	j = blockIdx.x * blockDim.x + threadIdx.x;

	if ( j >= n )
	{
		return;
	}




	bonds = &p_bonds;
	bond_list = bonds->select.bond_list;
	hbonds = &p_hbonds;
	hbond_list = hbonds->select.hbond_list;
	workspace = &p_workspace;

	/* loops below discover the Hydrogen bonds between i-j-k triplets.
	 * here j is H atom and there has to be some bond between i and j.
	 * Hydrogen bond is between j and k.
	 * so in this function i->X, j->H, k->Z when we map
	 * variables onto the ones in the handout. */
	//for( j = 0; j < system->n; ++j )
	if ( sbp[ my_atoms[j].type ].p_hbond == H_ATOM )
	{
		type_j = my_atoms[j].type;
		start_j = Cuda_Start_Index( j, bonds );
		end_j = Cuda_End_Index( j, bonds );
		hb_start_j = Cuda_Start_Index( my_atoms[j].Hindex, hbonds );
		hb_end_j = Cuda_End_Index( my_atoms[j].Hindex, hbonds );

		top = 0;
		for ( pi = start_j; pi < end_j; ++pi )
		{
			pbond_ij = &bond_list[pi];
			i = pbond_ij->nbr;
			bo_ij = &pbond_ij->bo_data;
			type_i = my_atoms[i].type;

			if ( sbp[type_i].p_hbond == H_BONDING_ATOM
					&& bo_ij->BO >= HB_THRESHOLD )
			{
				hblist[top++] = pi;
			}
		}

		for ( pk = hb_start_j; pk < hb_end_j; ++pk )
		{
			/* set k's varibles */
			k = hbond_list[pk].nbr;
			type_k = my_atoms[k].type;

			nbr_jk = hbond_list[pk].ptr;
			r_jk = nbr_jk->d;
			rvec_Scale( dvec_jk, hbond_list[pk].scl, nbr_jk->dvec );

			hbond_jk = &hbond_list[pk];
			rvec_MakeZero( hbond_jk->hb_f);

			/* find matching hbond to atom k */
			for ( itr = 0; itr < top; ++itr )
			{
				pi = hblist[itr];
				pbond_ij = &bonds->select.bond_list[pi];
				i = pbond_ij->nbr;

				if ( my_atoms[i].orig_id != my_atoms[k].orig_id )
				{
					bo_ij = &pbond_ij->bo_data;
					type_i = my_atoms[i].type;
					hbp = &d_hbp[ index_hbp(type_i, type_j, type_k, num_atom_types) ];

				}
			}

		}



	}
}


//CUDA_GLOBAL void __launch_bounds__ (256, 4) Cuda_Hydrogen_Bonds_MT ( reax_atom *my_atoms,
CUDA_GLOBAL void Cuda_Hydrogen_Bonds_MT( reax_atom *my_atoms, single_body_parameters *sbp,
		hbond_parameters *d_hbp, global_parameters gp, control_params *control,
		storage p_workspace, reax_list p_bonds, reax_list p_hbonds, int n,
		int num_atom_types, real *data_e_hb, rvec *data_ext_press )
{
#if defined( __SM_35__)
	real sh_hb;
	real sh_cdbo;
	rvec sh_atomf;
	rvec sh_hf;
#else
	HIP_DYNAMIC_SHARED( real, t_hb)
	HIP_DYNAMIC_SHARED( rvec, t__f)
	HIP_DYNAMIC_SHARED( rvec, t_cdbo)
	HIP_DYNAMIC_SHARED( rvec, t_hf)
	real *sh_hb = t_hb;
	real *sh_cdbo = t_hb + blockDim.x;
	rvec *sh_atomf = (rvec *)(sh_cdbo + blockDim.x);
	rvec *sh_hf = (rvec *)(sh_atomf + blockDim.x);
#endif
	int __THREADS_PER_ATOM__, thread_id, group_id, lane_id;
	int i, j, k, pi, pk;
	int type_i, type_j, type_k;
	int start_j, end_j, hb_start_j, hb_end_j;
	//TODO: re-write and remove
	int hblist[MAX_BONDS];
	int itr, top;
	int loopcount, count;
	ivec rel_jk;
	real r_jk, theta, cos_theta, sin_xhz4, cos_xhz1, sin_theta2;
	real e_hb, exp_hb2, exp_hb3, CEhb1, CEhb2, CEhb3;
	rvec dcos_theta_di, dcos_theta_dj, dcos_theta_dk;
	rvec dvec_jk, force, ext_press;
	hbond_parameters *hbp;
	bond_order_data *bo_ij;
	bond_data *pbond_ij;
	far_neighbor_data *nbr_jk;
	reax_list *bonds, *hbonds;
	bond_data *bond_list;
	hbond_data *hbond_list, *hbond_jk;
	storage *workspace;

	__THREADS_PER_ATOM__ = HB_KER_THREADS_PER_ATOM;
	thread_id = blockIdx.x * blockDim.x + threadIdx.x;
	group_id = thread_id / __THREADS_PER_ATOM__;
	lane_id = thread_id & (__THREADS_PER_ATOM__ - 1);

	if ( group_id >= n )
	{
		return;
	}

	workspace = &p_workspace;
	bonds = &p_bonds;
	bond_list = bonds->select.bond_list;
	hbonds = &p_hbonds;
	hbond_list = hbonds->select.hbond_list;
	j = group_id;

	/* loops below discover the Hydrogen bonds between i-j-k triplets.
       here j is H atom and there has to be some bond between i and j.
       Hydrogen bond is between j and k.
       so in this function i->X, j->H, k->Z when we map 
       variables onto the ones in the handout.*/
	//for( j = 0; j < system->n; ++j )

#if defined( __SM_35__)
	sh_hb = 0;
	rvec_MakeZero( sh_atomf );
#else
	sh_hb[threadIdx.x] = 0;
	rvec_MakeZero( sh_atomf[threadIdx.x] );
#endif

	/* j has to be of type H */
	if ( sbp[ my_atoms[j].type ].p_hbond == H_ATOM )
	{
		/* set j's variables */
		type_j = my_atoms[j].type;
		start_j = Cuda_Start_Index(j, bonds);
		end_j = Cuda_End_Index(j, bonds);
		hb_start_j = Cuda_Start_Index( my_atoms[j].Hindex, hbonds );
		hb_end_j = Cuda_End_Index( my_atoms[j].Hindex, hbonds );

		top = 0;
		for ( pi = start_j; pi < end_j; ++pi )
		{
			pbond_ij = &bond_list[pi];
			i = pbond_ij->nbr;
			bo_ij = &pbond_ij->bo_data;
			type_i = my_atoms[i].type;

			if ( sbp[type_i].p_hbond == H_BONDING_ATOM
					&& bo_ij->BO >= HB_THRESHOLD )
			{
				hblist[top++] = pi;
			}
		}

		//        fprintf( stderr, "j: %d, top: %d, hb_start_j: %d, hb_end_j:%d\n",
		//                j, top, hb_start_j, hb_end_j );

		for ( itr = 0; itr < top; ++itr )
		{
			pi = hblist[itr];
			pbond_ij = &bonds->select.bond_list[pi];
			i = pbond_ij->nbr;

#if defined( __SM_35__)
			rvec_MakeZero( sh_hf );
			sh_cdbo = 0;
#else
			rvec_MakeZero( sh_hf[threadIdx.x] );
			sh_cdbo[threadIdx.x] = 0;
#endif

			//for( pk = hb_start_j; pk < hb_end_j; ++pk ) {
			loopcount = (hb_end_j - hb_start_j) / HB_KER_THREADS_PER_ATOM +
					(((hb_end_j - hb_start_j) % HB_KER_THREADS_PER_ATOM == 0) ? 0 : 1);

			count = 0;
			pk = hb_start_j + lane_id;
			while ( count < loopcount )
			{
				/* only allow threads with an actual hbond */
				if ( pk < hb_end_j )
				{
					hbond_jk = &hbond_list[pk];

					/* set k's varibles */
					k = hbond_list[pk].nbr;
					type_k = my_atoms[k].type;
					nbr_jk = hbond_list[pk].ptr;
					r_jk = nbr_jk->d;
					rvec_Scale( dvec_jk, hbond_list[pk].scl, nbr_jk->dvec );
				}
				else
				{
					k = -1;
				}

				if ( my_atoms[i].orig_id != my_atoms[k].orig_id && k != -1 )
				{
					bo_ij = &pbond_ij->bo_data;
					type_i = my_atoms[i].type;
					hbp = &d_hbp[ index_hbp(type_i,type_j,type_k,num_atom_types) ];

					Calculate_Theta( pbond_ij->dvec, pbond_ij->d, dvec_jk, r_jk,
							&theta, &cos_theta );
					/* the derivative of cos(theta) */
					Calculate_dCos_Theta( pbond_ij->dvec, pbond_ij->d, dvec_jk, r_jk,
							&dcos_theta_di, &dcos_theta_dj, &dcos_theta_dk );

					/* hydrogen bond energy */
					sin_theta2 = SIN( theta / 2.0 );
					sin_xhz4 = SQR(sin_theta2);
					sin_xhz4 *= sin_xhz4;
					cos_xhz1 = ( 1.0 - cos_theta );
					exp_hb2 = EXP( -hbp->p_hb2 * bo_ij->BO );
					exp_hb3 = EXP( -hbp->p_hb3 * ( hbp->r0_hb / r_jk +
							r_jk / hbp->r0_hb - 2.0 ) );

					e_hb = hbp->p_hb1 * (1.0 - exp_hb2) * exp_hb3 * sin_xhz4;
#if defined( __SM_35__)
					sh_hb += e_hb;
#else
					sh_hb[threadIdx.x] += e_hb;
#endif

					CEhb1 = hbp->p_hb1 * hbp->p_hb2 * exp_hb2 * exp_hb3 * sin_xhz4;
					CEhb2 = -hbp->p_hb1/2.0 * (1.0 - exp_hb2) * exp_hb3 * cos_xhz1;
					CEhb3 = -hbp->p_hb3 *
							(-hbp->r0_hb / SQR(r_jk) + 1.0 / hbp->r0_hb) * e_hb;

					/*fprintf( stdout,
                      "%6d%6d%6d%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f\n",
                      system->my_atoms[i].orig_id, system->my_atoms[j].orig_id, 
                      system->my_atoms[k].orig_id, 
                      r_jk, theta, hbp->p_hb1, exp_hb2, hbp->p_hb3, hbp->r0_hb, 
                      exp_hb3, sin_xhz4, e_hb ); */

					/* hydrogen bond forces */
#if defined( __SM_35__)
					sh_cdbo += CEhb1; // dbo term
#else
					sh_cdbo[threadIdx.x] += CEhb1; // dbo term
#endif

					if ( control->virial == 0 )
					{
						// dcos terms
#if defined( __SM_35__)
						rvec_ScaledAdd( sh_hf, +CEhb2, dcos_theta_di );
#else
						rvec_ScaledAdd( sh_hf[threadIdx.x], +CEhb2, dcos_theta_di );
#endif

#if defined( __SM_35__)
						rvec_ScaledAdd( sh_atomf, +CEhb2, dcos_theta_dj );
#else
						rvec_ScaledAdd( sh_atomf[threadIdx.x], +CEhb2, dcos_theta_dj );
#endif

						rvec_ScaledAdd( hbond_jk->hb_f, +CEhb2, dcos_theta_dk );

						// dr terms
#if defined( __SM_35__)
						rvec_ScaledAdd( sh_atomf, -CEhb3/r_jk, dvec_jk );
#else
						rvec_ScaledAdd( sh_atomf[threadIdx.x], -CEhb3/r_jk, dvec_jk );
#endif

						rvec_ScaledAdd( hbond_jk->hb_f, +CEhb3/r_jk, dvec_jk );
					}
					else
					{
						/* for pressure coupling, terms that are not related to bond order
                           derivatives are added directly into pressure vector/tensor */
						rvec_Scale( force, +CEhb2, dcos_theta_di ); // dcos terms
						//rvec_Add( workspace->f[i], force );
						rvec_Add( pbond_ij->hb_f, force );
						rvec_iMultiply( ext_press, pbond_ij->rel_box, force );
						rvec_ScaledAdd( data_ext_press [j], 1.0, ext_press );

						rvec_ScaledAdd( workspace->f[j], +CEhb2, dcos_theta_dj );

						ivec_Scale( rel_jk, hbond_list[pk].scl, nbr_jk->rel_box );
						rvec_Scale( force, +CEhb2, dcos_theta_dk );
						//rvec_Add( workspace->f[k], force );
						rvec_Add( hbond_jk->hb_f, force );
						rvec_iMultiply( ext_press, rel_jk, force );
						rvec_ScaledAdd( data_ext_press[j], 1.0, ext_press );
						// dr terms
						rvec_ScaledAdd( workspace->f[j], -CEhb3/r_jk, dvec_jk );

						rvec_Scale( force, CEhb3/r_jk, dvec_jk );
						//rvec_Add( workspace->f[k], force );
						rvec_Add( hbond_jk->hb_f, force );
						rvec_iMultiply( ext_press, rel_jk, force );
						rvec_ScaledAdd( data_ext_press[j], 1.0, ext_press );
					}

				} //orid id end

				pk += __THREADS_PER_ATOM__;
				count++;

			} //for itr loop end

			//Reduction here
#if defined( __SM_35__)
			for ( int s = __THREADS_PER_ATOM__ >> 1; s >= 1; s/=2 )
			{
				sh_cdbo += shfl( sh_cdbo, s);
				sh_hf[0] += shfl( sh_hf[0], s);
				sh_hf[1] += shfl( sh_hf[1], s);
				sh_hf[2] += shfl( sh_hf[2], s);
			}
			//end of the shuffle
			if ( lane_id == 0 )
			{
				bo_ij->Cdbo += sh_cdbo ;
				rvec_Add( pbond_ij->hb_f, sh_hf );
			}
#else
			if ( lane_id < 16 )
			{
				sh_cdbo[threadIdx.x] += sh_cdbo[threadIdx.x + 16];
				rvec_Add( sh_hf [threadIdx.x], sh_hf[threadIdx.x + 16] );
			}
			if ( lane_id < 8 )
			{
				sh_cdbo[threadIdx.x] += sh_cdbo[threadIdx.x + 8];
				rvec_Add( sh_hf [threadIdx.x], sh_hf[threadIdx.x + 8] );
			}
			if ( lane_id < 4 )
			{
				sh_cdbo[threadIdx.x] += sh_cdbo[threadIdx.x + 4];
				rvec_Add( sh_hf [threadIdx.x], sh_hf[threadIdx.x + 4] );
			}
			if ( lane_id < 2 )
			{
				sh_cdbo[threadIdx.x] += sh_cdbo[threadIdx.x + 2];
				rvec_Add( sh_hf [threadIdx.x], sh_hf[threadIdx.x + 2] );
			}
			if ( lane_id < 1 )
			{
				sh_cdbo[threadIdx.x] += sh_cdbo[threadIdx.x + 1];
				rvec_Add( sh_hf [threadIdx.x], sh_hf[threadIdx.x + 1] );

				bo_ij->Cdbo += sh_cdbo[threadIdx.x];
				rvec_Add( pbond_ij->hb_f, sh_hf[threadIdx.x] );
			}
#endif
		} // for loop hbonds end
	} //if Hbond check end

#if defined( __SM_35__)
	for ( int s = __THREADS_PER_ATOM__ >> 1; s >= 1; s/=2 )
	{
		sh_hb += shfl( sh_hb, s);
		sh_atomf[0] += shfl( sh_atomf[0], s);
		sh_atomf[1] += shfl( sh_atomf[1], s);
		sh_atomf[2] += shfl( sh_atomf[2], s);
	}
	if ( lane_id == 0 )
	{
		data_e_hb[j] += sh_hb;
		rvec_Add( workspace->f[j], sh_atomf );
	}
#else
	if ( lane_id < 16 )
	{
		sh_hb[threadIdx.x] += sh_hb[threadIdx.x + 16];
		rvec_Add ( sh_atomf [threadIdx.x], sh_atomf[threadIdx.x + 16] );
	}
	if ( lane_id < 8 )
	{
		sh_hb[threadIdx.x] += sh_hb[threadIdx.x + 8];
		rvec_Add ( sh_atomf [threadIdx.x], sh_atomf[threadIdx.x + 8] );
	}
	if ( lane_id < 4 )
	{
		sh_hb[threadIdx.x] += sh_hb[threadIdx.x + 4];
		rvec_Add ( sh_atomf [threadIdx.x], sh_atomf[threadIdx.x + 4] );
	}
	if ( lane_id < 2 )
	{
		sh_hb[threadIdx.x] += sh_hb[threadIdx.x + 2];
		rvec_Add ( sh_atomf [threadIdx.x], sh_atomf[threadIdx.x + 2] );
	}
	if ( lane_id < 1 )
	{
		sh_hb[threadIdx.x] += sh_hb[threadIdx.x + 1];
		rvec_Add ( sh_atomf [threadIdx.x], sh_atomf[threadIdx.x + 1] );

		data_e_hb[j] += sh_hb[threadIdx.x];
		rvec_Add( workspace->f[j], sh_atomf[threadIdx.x] );
	}
#endif
}


CUDA_GLOBAL void Cuda_Hydrogen_Bonds_PostProcess( reax_atom *atoms,
		storage p_workspace, reax_list p_bonds, int N )
{
	int i, pj;
	storage *workspace;
	bond_data *pbond;
	bond_data *sym_index_bond;
	reax_list *bonds;

	i = blockIdx.x * blockDim.x + threadIdx.x;

	if ( i >= N )
	{
		return;
	}

	workspace = &p_workspace;
	bonds = &p_bonds;

	for ( pj = Cuda_Start_Index(i, bonds); pj < Cuda_End_Index(i, bonds); ++pj )
	{
		pbond = &bonds->select.bond_list[pj];
		sym_index_bond = &bonds->select.bond_list[pbond->sym_index];

		//rvec_Add( atoms[i].f, sym_index_bond->hb_f );
		rvec_Add( workspace->f[i], sym_index_bond->hb_f );
	}


}


CUDA_GLOBAL void Cuda_Hydrogen_Bonds_HNbrs( reax_atom *atoms,
		storage p_workspace, reax_list p_hbonds )
{
#if defined(__SM_35__)
	rvec __f;
#else
	HIP_DYNAMIC_SHARED( rvec, __f)
#endif
	int i, pj;
	int start, end;
	storage *workspace;
	hbond_data *nbr_pj, *sym_index_nbr;
	reax_list *hbonds;

	i = blockIdx.x;
	workspace = &p_workspace;
	hbonds = &p_hbonds;

	start = Cuda_Start_Index( atoms[i].Hindex, hbonds );
	end = Cuda_End_Index( atoms[i].Hindex, hbonds );
	pj = start + threadIdx.x;
#if defined(__SM_35__)
	rvec_MakeZero( __f );
#else
	rvec_MakeZero( __f[threadIdx.x] );
#endif

	while ( pj < end )
	{
		nbr_pj = &hbonds->select.hbond_list[pj];

		if (nbr_pj->sym_index == -1)
		{
			rvec hbf;
			rvec_MakeZero(hbf);

#if defined(__SM_35__) 
			rvec_Add( __f, hbf );
#else   
			rvec_Add( __f[threadIdx.x], hbf );
#endif

		}

		else
		{
			sym_index_nbr = &hbonds->select.hbond_list[ nbr_pj->sym_index ];

#if defined(__SM_35__)
			rvec_Add( __f, sym_index_nbr->hb_f );
#else
			rvec_Add( __f[threadIdx.x], sym_index_nbr->hb_f );
#endif
		}
		pj += blockDim.x;
	}

	__syncthreads( );

#if defined(__SM_35__)
	for ( int s = 16; s >= 1; s /= 2 )
	{
		__f[0] += shfl( __f[0], s );
		__f[1] += shfl( __f[1], s );
		__f[2] += shfl( __f[2], s );
	}

	if ( threadIdx.x == 0 )
	{
		rvec_Add( workspace->f[i], __f );
	}
#else
	if ( threadIdx.x < 16 )
	{
		rvec_Add( __f[threadIdx.x], __f[threadIdx.x + 16] );
	}
	__syncthreads( );

	if ( threadIdx.x < 8 )
	{
		rvec_Add( __f[threadIdx.x], __f[threadIdx.x + 8] );
	}
	__syncthreads( );

	if ( threadIdx.x < 4 )
	{
		rvec_Add( __f[threadIdx.x], __f[threadIdx.x + 4] );
	}
	__syncthreads( );

	if ( threadIdx.x < 2 )
	{
		rvec_Add( __f[threadIdx.x], __f[threadIdx.x + 2] );
	}
	__syncthreads( );

	if ( threadIdx.x < 1 )
	{
		rvec_Add( __f[threadIdx.x], __f[threadIdx.x + 1] );
	}
	__syncthreads( );

	if ( threadIdx.x == 0 )
	{
		//rvec_Add( atoms[i].f, __f[0] );

		rvec_Add( workspace->f[i], __f[0] );
		/*if(i < 5)
		{
			printf("%f,%f,%f\n", workspace->f[i][0],workspace->f[i][1],workspace->f[i][2]);
		}*/
	}
#endif

}


CUDA_GLOBAL void Cuda_Hydrogen_Bonds_HNbrs_BL( reax_atom *atoms,
		storage p_workspace, reax_list p_hbonds, int N )
{
#if defined(__SM_35__)
	rvec __f;
	int s;
#else
	HIP_DYNAMIC_SHARED( rvec, __f)
#endif
	int i, pj;
	int start, end;
	storage *workspace;
	hbond_data *nbr_pj, *sym_index_nbr;
	reax_list *hbonds;
	int __THREADS_PER_ATOM__;
	int thread_id;
	int group_id;
	int lane_id;

	__THREADS_PER_ATOM__ = HB_POST_PROC_KER_THREADS_PER_ATOM;
	thread_id = blockIdx.x * blockDim.x + threadIdx.x;
	group_id = thread_id / __THREADS_PER_ATOM__;
	lane_id = thread_id & (__THREADS_PER_ATOM__ - 1);

	if ( group_id >= N )
	{
		return;
	}

	workspace = &p_workspace;
	hbonds = &p_hbonds;
	i = group_id;
	start = Cuda_Start_Index( atoms[i].Hindex, hbonds );
	end = Cuda_End_Index( atoms[i].Hindex, hbonds );
	pj = start + lane_id;
#if defined(__SM_35__)
	rvec_MakeZero( __f );
#else
	rvec_MakeZero( __f[threadIdx.x] );
#endif

	while ( pj < end )
	{
		nbr_pj = &hbonds->select.hbond_list[pj];

		sym_index_nbr = &hbonds->select.hbond_list[ nbr_pj->sym_index ];
#if defined(__SM_35__)
		rvec_Add( __f, sym_index_nbr->hb_f );
#else
		rvec_Add( __f[threadIdx.x], sym_index_nbr->hb_f );
#endif

		pj += __THREADS_PER_ATOM__;
	}

	__syncthreads( );

#if defined(__SM_35__)
	for ( s = __THREADS_PER_ATOM__ >> 1; s >= 1; s /= 2 )
	{
		__f[0] += shfl( __f[0], s );
		__f[1] += shfl( __f[1], s );
		__f[2] += shfl( __f[2], s );
	}

	if ( lane_id == 0 )
	{
		rvec_Add( workspace->f[i], __f );
	}
#else
	if ( lane_id < 16 )
	{
		rvec_Add( __f[threadIdx.x], __f[threadIdx.x + 16] );
	}
	__syncthreads( );

	if ( lane_id < 8 )
	{
		rvec_Add( __f[threadIdx.x], __f[threadIdx.x + 8] );
	}
	__syncthreads( );

	if ( lane_id < 4 )
	{
		rvec_Add( __f[threadIdx.x], __f[threadIdx.x + 4] );
	}
	__syncthreads( );

	if ( lane_id < 2 )
	{
		rvec_Add( __f[threadIdx.x], __f[threadIdx.x + 2] );
	}
	__syncthreads( );

	if ( lane_id < 1 )
	{
		rvec_Add( __f[threadIdx.x], __f[threadIdx.x + 1] );
	}
	__syncthreads( );

	if ( lane_id == 0 )
	{
		rvec_Add( workspace->f[i], __f[threadIdx.x] );
	}
#endif
}
