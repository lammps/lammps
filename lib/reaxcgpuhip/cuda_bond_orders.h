
#ifndef __CUDA_BOND_ORDERS_H__
#define __CUDA_BOND_ORDERS_H__

#include "reaxc_types.h"

#include "vector.h"
#include "cuda_list.h"

extern "C" {

void Cuda_Total_Forces( reax_system *, control_params *,
		simulation_data *, storage *, reax_list ** );
void Cuda_Total_Forces_PURE( reax_system *, storage * );

}

CUDA_GLOBAL void Cuda_Calculate_BO_init( reax_atom *,
		single_body_parameters *, storage , int );

CUDA_GLOBAL void Cuda_Calculate_BO( reax_atom *, global_parameters,
		single_body_parameters *, two_body_parameters *,
		storage , reax_list , int , int );

CUDA_GLOBAL void Cuda_Update_Uncorrected_BO( storage , reax_list , int );

CUDA_GLOBAL void Cuda_Update_Workspace_After_BO( reax_atom *, global_parameters ,
		single_body_parameters *, storage , int );

CUDA_DEVICE static inline int Cuda_BOp( reax_list bonds, real bo_cut,
		int i, int btop_i, far_neighbor_data *nbr_pj,
		single_body_parameters *sbp_i, single_body_parameters *sbp_j,
		two_body_parameters *twbp, rvec *dDeltap_self, real *total_bond_order )
{

	int j, btop_j;
	real r2, C12, C34, C56;
	real Cln_BOp_s, Cln_BOp_pi, Cln_BOp_pi2;
	real BO, BO_s, BO_pi, BO_pi2;
	bond_data *ibond, *jbond;
	bond_order_data *bo_ij, *bo_ji;
	rvec bo_ij_dln_BOp_s;
	rvec bo_ij_dln_BOp_pi;
	rvec bo_ij_dln_BOp_pi2;
	rvec bo_ij_dBOp;

	j = nbr_pj->nbr;
	r2 = SQR( nbr_pj->d );


	if ( sbp_i->r_s > 0.0 && sbp_j->r_s > 0.0 )
	{
		C12 = twbp->p_bo1 * POW( nbr_pj->d / twbp->r_s, twbp->p_bo2 );
		BO_s = (1.0 + bo_cut) * EXP( C12 );
	}
	else
	{
		C12 = 0.0;
		BO_s = 0.0;
	}

	if ( sbp_i->r_pi > 0.0 && sbp_j->r_pi > 0.0 )
	{
		C34 = twbp->p_bo3 * POW( nbr_pj->d / twbp->r_p, twbp->p_bo4 );
		BO_pi = EXP( C34 );
	}
	else
	{
		C34 = 0.0;
		BO_pi = 0.0;
	}

	if ( sbp_i->r_pi_pi > 0.0 && sbp_j->r_pi_pi > 0.0 )
	{
		C56 = twbp->p_bo5 * POW( nbr_pj->d / twbp->r_pp, twbp->p_bo6 );
		BO_pi2 = EXP( C56 );
	}
	else
	{
		C56 = 0.0;
		BO_pi2 = 0.0;
	}

	/* Initially BO values are the uncorrected ones, page 1 */
	BO = BO_s + BO_pi + BO_pi2;



	if ( BO >= bo_cut )
	{
		/****** bonds i-j and j-i ******/

		/* Bond Order page2-3, derivative of total bond order prime */
		Cln_BOp_s = twbp->p_bo2 * C12 / r2;
		Cln_BOp_pi = twbp->p_bo4 * C34 / r2;
		Cln_BOp_pi2 = twbp->p_bo6 * C56 / r2;




		if ( i < j )
		{

			ibond = &( bonds.select.bond_list[btop_i] );
			ibond->nbr = j;
			ibond->d = nbr_pj->d;
			rvec_Copy( ibond->dvec, nbr_pj->dvec );
			btop_j = Cuda_End_Index( j, &bonds );

			/*if(i == 0 || i == 14)
			{
				printf("when i %d, J %d, btop %d,%d\n", i, j, btop_i,btop_j);
			}*/

			/*if(i == 14 )
			{
				bond_data *temp_bond;
				temp_bond = &( bonds.select.bond_list[350]);

				printf("%f,%f,%f\n", temp_bond->dvec[0], temp_bond->dvec[1], temp_bond->dvec[2]);
				//printf("Less %d,%d,%d,%f,%f,%f\n", btop_i,i,j,nbr_pj->dvec[0],nbr_pj->dvec[1],nbr_pj->dvec[2]);
			}*/
			ivec_Copy( ibond->rel_box, nbr_pj->rel_box );

			//ibond->dbond_index = btop_i;
			//ibond->sym_index = btop_j;

			bo_ij = &( ibond->bo_data );
			bo_ij->BO = BO;
			bo_ij->BO_s = BO_s;
			bo_ij->BO_pi = BO_pi;
			bo_ij->BO_pi2 = BO_pi2;

			/* Only dln_BOp_xx wrt. dr_i is stored here, note that
			 * dln_BOp_xx/dr_i = -dln_BOp_xx/dr_j and all others are 0 */
			rvec_Scale( bo_ij->dln_BOp_s,
					-bo_ij->BO_s * Cln_BOp_s, ibond->dvec );
			rvec_Scale( bo_ij->dln_BOp_pi,
					-bo_ij->BO_pi * Cln_BOp_pi, ibond->dvec );
			rvec_Scale( bo_ij->dln_BOp_pi2,
					-bo_ij->BO_pi2 * Cln_BOp_pi2, ibond->dvec );

			/* Only dBOp wrt. dr_i is stored here, note that
			 * dBOp/dr_i = -dBOp/dr_j and all others are 0 */
			rvec_Scale( bo_ij->dBOp, -(bo_ij->BO_s * Cln_BOp_s +
					bo_ij->BO_pi * Cln_BOp_pi + bo_ij->BO_pi2 *
					Cln_BOp_pi2), ibond->dvec );

			rvec_Add( dDeltap_self[i], bo_ij->dBOp );

			bo_ij->BO_s -= bo_cut;
			bo_ij->BO -= bo_cut;



			//currently total_BOp
			total_bond_order[i] += bo_ij->BO;



			bo_ij->Cdbo = bo_ij->Cdbopi = bo_ij->Cdbopi2 = 0.0;

			//CUDA Specific
			ibond->ae_CdDelta = 0;
			ibond->va_CdDelta = 0;
			rvec_MakeZero( ibond->va_f );
			ibond->ta_CdDelta = 0;
			ibond->ta_Cdbo = 0;
			rvec_MakeZero( ibond->ta_f );
			rvec_MakeZero( ibond->hb_f );
			rvec_MakeZero( ibond->tf_f );
		}
		else
		{

			//btop_j = End_Index( j, bonds );
			btop_j = btop_i;

			jbond = &(bonds.select.bond_list[btop_j]);
			//jbond->nbr = i;
			jbond->nbr = j;
			jbond->d = nbr_pj->d;
			rvec_Scale( jbond->dvec, -1, nbr_pj->dvec );
			ivec_Scale( jbond->rel_box, -1, nbr_pj->rel_box );

			/*if(i == 14)
			{
				bond_data *temp_bond;
				temp_bond = &( bonds.select.bond_list[350]);

				printf("%f,%f,%f\n", temp_bond->dvec[0], temp_bond->dvec[1], temp_bond->dvec[2]);
				//printf(" Greater %d,%d,%d,%f,%f,%f\n",btop_j, i,j,jbond->dvec[0],jbond->dvec[1],jbond->dvec[2]);
			}*/


			//jbond->dbond_index = btop_i;
			//jbond->sym_index = btop_i;

			//Set_End_Index( j, btop_j + 1, bonds );

			bo_ji = &( jbond->bo_data );
			bo_ji->BO = BO;
			bo_ji->BO_s = BO_s;
			bo_ji->BO_pi = BO_pi;
			bo_ji->BO_pi2 = BO_pi2;

			/* Only dln_BOp_xx wrt. dr_i is stored here, note that
            dln_BOp_xx/dr_i = -dln_BOp_xx/dr_j and all others are 0 */

			rvec_Scale( bo_ij_dln_BOp_s, -BO_s * Cln_BOp_s, nbr_pj->dvec );
			rvec_Scale( bo_ij_dln_BOp_pi, -BO_pi * Cln_BOp_pi, nbr_pj->dvec );
			rvec_Scale( bo_ij_dln_BOp_pi2,
					-BO_pi2 * Cln_BOp_pi2, nbr_pj->dvec );
			rvec_Scale( bo_ji->dln_BOp_s, -1., bo_ij_dln_BOp_s );
			rvec_Scale( bo_ji->dln_BOp_pi, -1., bo_ij_dln_BOp_pi );
			rvec_Scale( bo_ji->dln_BOp_pi2, -1., bo_ij_dln_BOp_pi2 );

			/* Only dBOp wrt. dr_i is stored here, note that
            dBOp/dr_i = -dBOp/dr_j and all others are 0 */
			rvec_Scale( bo_ij_dBOp, -(BO_s * Cln_BOp_s +
					BO_pi * Cln_BOp_pi +
					BO_pi2 * Cln_BOp_pi2), nbr_pj->dvec );
			rvec_Scale( bo_ji->dBOp, -1.0, bo_ij_dBOp );

			rvec_Add( dDeltap_self[i], bo_ji->dBOp );

			bo_ji->BO_s -= bo_cut;
			bo_ji->BO -= bo_cut;

			total_bond_order[i] += bo_ji->BO; //currently total_BOp





			bo_ji->Cdbo = bo_ji->Cdbopi = bo_ji->Cdbopi2 = 0.0;

			//CUDA Specific
			jbond->ae_CdDelta = 0;
			jbond->va_CdDelta = 0;
			rvec_MakeZero( jbond->va_f );
			jbond->ta_CdDelta = 0;
			jbond->ta_Cdbo = 0;
			rvec_MakeZero( jbond->ta_f );
			rvec_MakeZero( jbond->hb_f );
			rvec_MakeZero( jbond->tf_f );
		}




		return TRUE;
	}

	return FALSE;
}

#endif
