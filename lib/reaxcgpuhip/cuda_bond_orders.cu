
#include "cuda_bond_orders.h"

#include "cuda_list.h"
#include "cuda_utils.h"
#include "cuda_reduction.h"

#include "index_utils.h"
#include "bond_orders.h"


CUDA_GLOBAL void Cuda_Calculate_BO_init( reax_atom *my_atoms, 
		single_body_parameters *sbp, storage p_workspace, int N )
{
	int i, type_i;
	single_body_parameters *sbp_i;
	storage *workspace;

	i = blockIdx.x * blockDim.x + threadIdx.x;

	if ( i >= N )
	{
		return;
	}

	workspace = &p_workspace;

	/* Calculate Deltaprime, Deltaprime_boc values */
	type_i = my_atoms[i].type;
	sbp_i = &sbp[type_i];
	workspace->Deltap[i] = workspace->total_bond_order[i] - sbp_i->valency;
	workspace->Deltap_boc[i] = workspace->total_bond_order[i]
														   - sbp_i->valency_val;
	workspace->total_bond_order[i] = 0;
}


CUDA_GLOBAL void Cuda_Calculate_BO( reax_atom *my_atoms, global_parameters gp, 
		single_body_parameters *sbp, two_body_parameters *tbp,
		storage p_workspace, reax_list p_bonds,
		int num_atom_types, int N )
{
	int i, j, pj, type_i, type_j;
	int start_i, end_i;
	//    int sym_index;
	real val_i, Deltap_i, Deltap_boc_i;
	real val_j, Deltap_j, Deltap_boc_j;
	real f1, f2, f3, f4, f5, f4f5, exp_f4, exp_f5;
	real exp_p1i, exp_p2i, exp_p1j, exp_p2j;
	real temp, u1_ij, u1_ji, Cf1A_ij, Cf1B_ij, Cf1_ij, Cf1_ji;
	real Cf45_ij, Cf45_ji;
	real A0_ij, A1_ij, A2_ij, A2_ji, A3_ij, A3_ji;
	real p_boc1, p_boc2;
	single_body_parameters *sbp_i;
	two_body_parameters *twbp;
	bond_order_data *bo_ij;
	//    bond_order_data *bo_ji;
	storage *workspace;
	reax_list *bonds;

	i = blockIdx.x * blockDim.x + threadIdx.x;

	if ( i >= N )
	{
		return;
	}

	workspace = &p_workspace;
	bonds = &p_bonds;
	p_boc1 = gp.l[0];
	p_boc2 = gp.l[1];

	/* Corrected Bond Order calculations */
	//for( i = 0; i < system->N; ++i ) {
	type_i = my_atoms[i].type;
	sbp_i = &sbp[type_i];
	val_i = sbp_i->valency;
	Deltap_i = workspace->Deltap[i];
	Deltap_boc_i = workspace->Deltap_boc[i];
	start_i = Cuda_Start_Index( i, bonds );
	end_i = Cuda_End_Index( i, bonds );

	for( pj = start_i; pj < end_i; ++pj )
	{
		j = bonds->select.bond_list[pj].nbr;
		type_j = my_atoms[j].type;
		bo_ij = &bonds->select.bond_list[pj].bo_data;
		//TODO
		//if( i < j || workspace->bond_mark[j] > 3 ) {
		if( i < j )
		{
			twbp = &tbp[ index_tbp(type_i, type_j, num_atom_types)];

#ifdef TEST_FORCES
			Set_Start_Index( pj, top_dbo, dBOs );

			/* fprintf( stderr, "%6d%6d%12.6f%12.6f%12.6f\n",
               workspace->reverse_map[i], workspace->reverse_map[j],
               twbp->ovc, twbp->v13cor, bo_ij->BO ); */
#endif

			if ( twbp->ovc < 0.001 && twbp->v13cor < 0.001 )
			{
				/* There is no correction to bond orders nor to derivatives
				 * of bond order prime! So we leave bond orders unchanged and
				 * set derivative of bond order coefficients such that
				 * dBO = dBOp & dBOxx = dBOxxp in Add_dBO_to_Forces */
				bo_ij->C1dbo = 1.000000;
				bo_ij->C2dbo = 0.000000;
				bo_ij->C3dbo = 0.000000;

				bo_ij->C1dbopi = 1.0;
				bo_ij->C2dbopi = 0.0;
				bo_ij->C3dbopi = 0.0;
				bo_ij->C4dbopi = 0.0;

				bo_ij->C1dbopi2 = 1.0;
				bo_ij->C2dbopi2 = 0.0;
				bo_ij->C3dbopi2 = 0.0;
				bo_ij->C4dbopi2 = 0.0;

#ifdef TEST_FORCES
				pdbo = &dBOs->dbo_list[ top_dbo ];

				// compute dBO_ij/dr_i
				pdbo->wrt = i;
				rvec_Copy( pdbo->dBO, bo_ij->dBOp );
				rvec_Scale( pdbo->dBOpi, bo_ij->BO_pi, bo_ij->dln_BOp_pi );
				rvec_Scale( pdbo->dBOpi2, bo_ij->BO_pi2, bo_ij->dln_BOp_pi2);

				// compute dBO_ij/dr_j
				pdbo++;
				pdbo->wrt = j;
				rvec_Scale( pdbo->dBO, -1.0, bo_ij->dBOp );
				rvec_Scale( pdbo->dBOpi, -bo_ij->BO_pi, bo_ij->dln_BOp_pi );
				rvec_Scale(pdbo->dBOpi2, -bo_ij->BO_pi2, bo_ij->dln_BOp_pi2);

				top_dbo += 2;
#endif
			}
			else
			{
				val_j = sbp[type_j].valency;
				Deltap_j = workspace->Deltap[j];
				Deltap_boc_j = workspace->Deltap_boc[j];

				/* on page 1 */
				if ( twbp->ovc >= 0.001 )
				{
					/* Correction for overcoordination */
					exp_p1i = EXP( -p_boc1 * Deltap_i );
					exp_p2i = EXP( -p_boc2 * Deltap_i );
					exp_p1j = EXP( -p_boc1 * Deltap_j );
					exp_p2j = EXP( -p_boc2 * Deltap_j );

					f2 = exp_p1i + exp_p1j;
					f3 = -1.0 / p_boc2 * LOG( 0.5 * ( exp_p2i  + exp_p2j ) );
					f1 = 0.5 * ( ( val_i + f2 )/( val_i + f2 + f3 ) +
							( val_j + f2 )/( val_j + f2 + f3 ) );

					/* Now come the derivates */
					/* Bond Order pages 5-7, derivative of f1 */
					temp = f2 + f3;
					u1_ij = val_i + temp;
					u1_ji = val_j + temp;
					Cf1A_ij = 0.5 * f3 * (1.0 / SQR( u1_ij ) +
							1.0 / SQR( u1_ji ));
					Cf1B_ij = -0.5 * (( u1_ij - f3 ) / SQR( u1_ij ) +
							( u1_ji - f3 ) / SQR( u1_ji ));

					//Cf1_ij = -Cf1A_ij * p_boc1 * exp_p1i +
					//          Cf1B_ij * exp_p2i / ( exp_p2i + exp_p2j );
					Cf1_ij = 0.50 * ( -p_boc1 * exp_p1i / u1_ij -
							((val_i+f2) / SQR(u1_ij)) *
							( -p_boc1 * exp_p1i +
									exp_p2i / ( exp_p2i + exp_p2j ) ) +
									-p_boc1 * exp_p1i / u1_ji -
									((val_j+f2) / SQR(u1_ji)) *
									( -p_boc1 * exp_p1i +
											exp_p2i / ( exp_p2i + exp_p2j ) ));

					Cf1_ji = -Cf1A_ij * p_boc1 * exp_p1j +
							Cf1B_ij * exp_p2j / ( exp_p2i + exp_p2j );

					//fprintf( stderr, "\tCf1:%g  %g\n", Cf1_ij, Cf1_ji );
				}
				else
				{
					/* No overcoordination correction! */
					f1 = 1.0;
					Cf1_ij = Cf1_ji = 0.0;
				}

				if ( twbp->v13cor >= 0.001 )
				{
					/* Correction for 1-3 bond orders */
					exp_f4 = EXP(-(twbp->p_boc4 * SQR( bo_ij->BO ) -
							Deltap_boc_i) * twbp->p_boc3 + twbp->p_boc5);
					exp_f5 = EXP(-(twbp->p_boc4 * SQR( bo_ij->BO ) -
							Deltap_boc_j) * twbp->p_boc3 + twbp->p_boc5);

					f4 = 1.0 / (1.0 + exp_f4);
					f5 = 1.0 / (1.0 + exp_f5);
					f4f5 = f4 * f5;

					/* Bond Order pages 8-9, derivative of f4 and f5 */
					Cf45_ij = -f4 * exp_f4;
					Cf45_ji = -f5 * exp_f5;
				}
				else
				{
					f4 = f5 = f4f5 = 1.0;
					Cf45_ij = Cf45_ji = 0.0;
				}

				/* Bond Order page 10, derivative of total bond order */
				A0_ij = f1 * f4f5;
				A1_ij = -2 * twbp->p_boc3 * twbp->p_boc4 * bo_ij->BO *
						(Cf45_ij + Cf45_ji);
				A2_ij = Cf1_ij / f1 + twbp->p_boc3 * Cf45_ij;
				A2_ji = Cf1_ji / f1 + twbp->p_boc3 * Cf45_ji;
				A3_ij = A2_ij + Cf1_ij / f1;
				A3_ji = A2_ji + Cf1_ji / f1;

				/*fprintf( stderr, "\tBO: %f, A0: %f, A1: %f"
                  "A2_ij: %f A2_ji: %f, A3_ij: %f, A3_ji: %f\n",
                  bo_ij->BO, 
                  A0_ij, A1_ij, A2_ij, A2_ji, A3_ij, A3_ji );*/

				/* find corrected bond orders and their derivative coef */

				/*if(i == 18 || i == 2 || i == 15 || i == 4 || i== 14)
				{
					printf("i:%d,j:%d,boij:%f,aoij:%f\n", i,j,bo_ij->BO,A0_ij );
				}*/
				bo_ij->BO = bo_ij->BO * A0_ij;
				bo_ij->BO_pi = bo_ij->BO_pi * A0_ij *f1;
				bo_ij->BO_pi2 = bo_ij->BO_pi2* A0_ij *f1;
				bo_ij->BO_s = bo_ij->BO - ( bo_ij->BO_pi + bo_ij->BO_pi2 );

				bo_ij->C1dbo = A0_ij + bo_ij->BO * A1_ij;



				bo_ij->C2dbo = bo_ij->BO * A2_ij;
				bo_ij->C3dbo = bo_ij->BO * A2_ji;

				bo_ij->C1dbopi = f1 * f1 * f4 * f5;
				bo_ij->C2dbopi = bo_ij->BO_pi * A1_ij;
				bo_ij->C3dbopi = bo_ij->BO_pi * A3_ij;
				bo_ij->C4dbopi = bo_ij->BO_pi * A3_ji;

				bo_ij->C1dbopi2 = f1 * f1 * f4 * f5;
				bo_ij->C2dbopi2 = bo_ij->BO_pi2 * A1_ij;
				bo_ij->C3dbopi2 = bo_ij->BO_pi2 * A3_ij;
				bo_ij->C4dbopi2 = bo_ij->BO_pi2 * A3_ji;

				//CHANGE ORIGINAL
			}
			//CHANGE ORIGINAL

			/* neglect bonds that are < 1e-10 */
			if ( bo_ij->BO < 1e-10 )
			{
				bo_ij->BO = 0.0;
			}
			if ( bo_ij->BO_s < 1e-10 )
			{
				bo_ij->BO_s = 0.0;
			}
			if ( bo_ij->BO_pi < 1e-10 )
			{
				bo_ij->BO_pi = 0.0;
			}
			if ( bo_ij->BO_pi2 < 1e-10 )
			{
				bo_ij->BO_pi2 = 0.0;
			}



			workspace->total_bond_order[i] += bo_ij->BO; //now keeps total_BO



			// printf("%d,%f\n",i, workspace->total_bond_order[i]);

			/* fprintf( stderr, "%d %d\t%g %g %g %g\n"
               "Cdbo:\t%g %g %g\n"
               "Cdbopi:\t%g %g %g %g\n"
               "Cdbopi2:%g %g %g %g\n\n", 
               i+1, j+1, 
               bonds->bond_list[ pj ].d, 
               bo_ij->BO,bo_ij->BO_pi, bo_ij->BO_pi2, 
               bo_ij->C1dbo, bo_ij->C2dbo, bo_ij->C3dbo,
               bo_ij->C1dbopi, bo_ij->C2dbopi, 
               bo_ij->C3dbopi, bo_ij->C4dbopi,
               bo_ij->C1dbopi2,bo_ij->C2dbopi2, 
               bo_ij->C3dbopi2, bo_ij->C4dbopi2 ); */

			/* fprintf( stderr, "%d %d  BO:%f BO_s:%f BO_pi:%f BO_pi2:%f\n",
               i+1,j+1,bo_ij->BO,bo_ij->BO_s,bo_ij->BO_pi,bo_ij->BO_pi2 );*/

#ifdef TEST_FORCES
			Set_End_Index( pj, top_dbo, dBOs );
			Add_dBO( system, lists, i, pj, 1.0, workspace->dDelta );
#endif

			//CHANGE ORIGINAL
			//}
			//CHANGE ORIGINAL
			/*
               else {
            // We only need to update bond orders from bo_ji
            //   everything else is set in uncorrected_bo calculations
            sym_index = bonds->bond_list[pj].sym_index;
            bo_ji = &bonds->bond_list[ sym_index ].bo_data;
            bo_ij->BO = bo_ji->BO;
            bo_ij->BO_s = bo_ji->BO_s;
            bo_ij->BO_pi = bo_ji->BO_pi;
            bo_ij->BO_pi2 = bo_ji->BO_pi2;

            workspace->total_bond_order[i] += bo_ij->BO;// now keeps total_BO
#ifdef TEST_FORCES
            Add_dBO( system, lists, j, sym_index, 1.0, workspace->dDelta );
#endif
}
			 */
		}
	}
	//if(i == 6)
	//{


}


CUDA_GLOBAL void Cuda_Update_Uncorrected_BO( storage p_workspace,
		reax_list p_bonds, int N )
{
	int i, j, pj;
	int start_i, end_i;
	int sym_index;
	storage *workspace;
	reax_list *bonds;
	bond_order_data *bo_ij, *bo_ji;

	i = blockIdx.x * blockDim.x + threadIdx.x;

	if ( i >= N )
	{
		return;
	}

	workspace = &p_workspace;
	bonds = &p_bonds;
	start_i = Cuda_Start_Index( i, bonds );
	end_i = Cuda_End_Index( i, bonds );

	for( pj = start_i; pj < end_i; ++pj )
	{

		j = bonds->select.bond_list[pj].nbr;
		bo_ij = &bonds->select.bond_list[pj].bo_data;

		//if( (i >= j)  || (workspace->bond_mark [i] <= 3)) {
		if ( i >= j )
		{
			/* We only need to update bond orders from bo_ji
               everything else is set in uncorrected_bo calculations */
			sym_index = bonds->select.bond_list[pj].sym_index;
			bo_ji = &bonds->select.bond_list[ sym_index ].bo_data;
			bo_ij->BO = bo_ji->BO;
			bo_ij->BO_s = bo_ji->BO_s;
			bo_ij->BO_pi = bo_ji->BO_pi;
			bo_ij->BO_pi2 = bo_ji->BO_pi2;
			// now keeps total_BO
			workspace->total_bond_order[i] += bo_ij->BO;
		}
	}

}


CUDA_GLOBAL void Cuda_Update_Workspace_After_BO( reax_atom *my_atoms,
		global_parameters gp, single_body_parameters *sbp,
		storage p_workspace, int N )
{
	int j, type_j;
	real explp1, p_lp1;
	single_body_parameters *sbp_j;
	storage *workspace;

	j = blockIdx.x * blockDim.x + threadIdx.x;

	if ( j >= N )
	{
		return;
	}

	workspace = &p_workspace;
	p_lp1 = gp.l[15];

	/* Calculate some helper variables that are  used at many places
       throughout force calculations */
	//for( j = 0; j < system->N; ++j ){
	type_j = my_atoms[j].type;
	sbp_j = &sbp[ type_j ];

	workspace->Delta[j] = workspace->total_bond_order[j] - sbp_j->valency;
	workspace->Delta_e[j] = workspace->total_bond_order[j] - sbp_j->valency_e;
	workspace->Delta_boc[j] = workspace->total_bond_order[j]
														  - sbp_j->valency_boc;

	workspace->vlpex[j] = workspace->Delta_e[j] -
			2.0 * (int)(workspace->Delta_e[j]/2.0);
	explp1 = EXP(-p_lp1 * SQR(2.0 + workspace->vlpex[j]));
	workspace->nlp[j] = explp1 - (int)(workspace->Delta_e[j] / 2.0);
	workspace->Delta_lp[j] = sbp_j->nlp_opt - workspace->nlp[j];
	workspace->Clp[j] = 2.0 * p_lp1 * explp1 * (2.0 + workspace->vlpex[j]);
	/* Adri uses different dDelta_lp values than the ones in notes... */
	workspace->dDelta_lp[j] = workspace->Clp[j];
	//workspace->dDelta_lp[j] = workspace->Clp[j] + (0.5-workspace->Clp[j]) *
	//((FABS(workspace->Delta_e[j]/2.0 -
	//       (int)(workspace->Delta_e[j]/2.0)) < 0.1) ? 1 : 0 );

	if( sbp_j->mass > 21.0 )
	{
		workspace->nlp_temp[j] = 0.5 * (sbp_j->valency_e - sbp_j->valency);
		workspace->Delta_lp_temp[j] = sbp_j->nlp_opt - workspace->nlp_temp[j];
		workspace->dDelta_lp_temp[j] = 0.;
	}
	else
	{
		workspace->nlp_temp[j] = workspace->nlp[j];
		workspace->Delta_lp_temp[j] = sbp_j->nlp_opt - workspace->nlp_temp[j];
		workspace->dDelta_lp_temp[j] = workspace->Clp[j];
	}

}


CUDA_DEVICE void Cuda_Add_dBond_to_Forces_NPT( int i, int pj,
		simulation_data *data, storage *workspace, reax_list *bonds,
		rvec data_ext_press )
{
	bond_data *nbr_j, *nbr_k;
	bond_order_data *bo_ij, *bo_ji;
	dbond_coefficients coef;
	rvec temp, ext_press;
	ivec rel_box;
	int pk, k, j;

	/* Initializations */
	nbr_j = &bonds->select.bond_list[pj];
	j = nbr_j->nbr;

	//bo_ij = &nbr_j->bo_data;
	//bo_ji = &bonds->bond_list[ nbr_j->sym_index ].bo_data;



	if (i < j)
	{
		bo_ij = &nbr_j->bo_data;
		bo_ji = &bonds->select.bond_list[ nbr_j->sym_index ].bo_data;
	}
	else
	{
		bo_ji = &nbr_j->bo_data;
		bo_ij = &bonds->select.bond_list[ nbr_j->sym_index ].bo_data;
	}

	coef.C1dbo = bo_ij->C1dbo * (bo_ij->Cdbo + bo_ji->Cdbo);
	coef.C2dbo = bo_ij->C2dbo * (bo_ij->Cdbo + bo_ji->Cdbo);
	coef.C3dbo = bo_ij->C3dbo * (bo_ij->Cdbo + bo_ji->Cdbo);

	coef.C1dbopi = bo_ij->C1dbopi * (bo_ij->Cdbopi + bo_ji->Cdbopi);
	coef.C2dbopi = bo_ij->C2dbopi * (bo_ij->Cdbopi + bo_ji->Cdbopi);
	coef.C3dbopi = bo_ij->C3dbopi * (bo_ij->Cdbopi + bo_ji->Cdbopi);
	coef.C4dbopi = bo_ij->C4dbopi * (bo_ij->Cdbopi + bo_ji->Cdbopi);

	coef.C1dbopi2 = bo_ij->C1dbopi2 * (bo_ij->Cdbopi2 + bo_ji->Cdbopi2);
	coef.C2dbopi2 = bo_ij->C2dbopi2 * (bo_ij->Cdbopi2 + bo_ji->Cdbopi2);
	coef.C3dbopi2 = bo_ij->C3dbopi2 * (bo_ij->Cdbopi2 + bo_ji->Cdbopi2);
	coef.C4dbopi2 = bo_ij->C4dbopi2 * (bo_ij->Cdbopi2 + bo_ji->Cdbopi2);

	coef.C1dDelta = bo_ij->C1dbo * (workspace->CdDelta[i]+workspace->CdDelta[j]);
	coef.C2dDelta = bo_ij->C2dbo * (workspace->CdDelta[i]+workspace->CdDelta[j]);
	coef.C3dDelta = bo_ij->C3dbo * (workspace->CdDelta[i]+workspace->CdDelta[j]);

	/************************************
	 * forces related to atom i          *
	 * first neighbors of atom i         *
	 ************************************/
	if ( i < j )
	{
		for ( pk = Cuda_Start_Index(i, bonds); pk < Cuda_End_Index(i, bonds); ++pk )
		{
			nbr_k = &bonds->select.bond_list[pk];
			k = nbr_k->nbr;

			rvec_MakeZero( nbr_k->tf_f );

			rvec_Scale( temp, -coef.C2dbo, nbr_k->bo_data.dBOp );       /*2nd, dBO*/
			rvec_ScaledAdd( temp, -coef.C2dDelta, nbr_k->bo_data.dBOp );/*dDelta*/
			rvec_ScaledAdd( temp, -coef.C3dbopi, nbr_k->bo_data.dBOp ); /*3rd, dBOpi*/
			rvec_ScaledAdd( temp, -coef.C3dbopi2, nbr_k->bo_data.dBOp );/*3rd, dBOpi2*/

			/* force */
			rvec_Add( nbr_k->tf_f, temp );
			/* pressure */
			rvec_iMultiply( ext_press, nbr_k->rel_box, temp );
			rvec_Add( data_ext_press, ext_press );

			/* if( !ivec_isZero( nbr_k->rel_box ) )
               fprintf( stderr, "%3d %3d %3d: dvec[%10.6f %10.6f %10.6f]"
               "ext[%3d %3d %3d] f[%10.6f %10.6f %10.6f]\n",
               i+1, system->my_atoms[i].x[0], 
               system->my_atoms[i].x[1], system->my_atoms[i].x[2], 
               j+1, k+1, system->my_atoms[k].x[0], 
               system->my_atoms[k].x[1], system->my_atoms[k].x[2],
               nbr_k->dvec[0], nbr_k->dvec[1], nbr_k->dvec[2],
               nbr_k->rel_box[0], nbr_k->rel_box[1], nbr_k->rel_box[2],
               temp[0], temp[1], temp[2] ); */
		}

		/* then atom i itself  */
		rvec_Scale( temp, coef.C1dbo, bo_ij->dBOp );                      /*1st,dBO*/
		rvec_ScaledAdd( temp, coef.C2dbo, workspace->dDeltap_self[i] );   /*2nd,dBO*/
		rvec_ScaledAdd( temp, coef.C1dDelta, bo_ij->dBOp );               /*1st,dBO*/
		rvec_ScaledAdd( temp, coef.C2dDelta, workspace->dDeltap_self[i] );/*2nd,dBO*/
		rvec_ScaledAdd( temp, coef.C1dbopi, bo_ij->dln_BOp_pi );        /*1st,dBOpi*/
		rvec_ScaledAdd( temp, coef.C2dbopi, bo_ij->dBOp );              /*2nd,dBOpi*/
		rvec_ScaledAdd( temp, coef.C3dbopi, workspace->dDeltap_self[i]);/*3rd,dBOpi*/

		rvec_ScaledAdd( temp, coef.C1dbopi2, bo_ij->dln_BOp_pi2 );  /*1st,dBO_pi2*/
		rvec_ScaledAdd( temp, coef.C2dbopi2, bo_ij->dBOp );         /*2nd,dBO_pi2*/
		rvec_ScaledAdd( temp, coef.C3dbopi2, workspace->dDeltap_self[i] );/*3rd*/

		/* force */
		rvec_Add( workspace->f[i], temp );
		/* ext pressure due to i is dropped, counting force on j will be enough */
	}
	else
	{
		/******************************************************
		 * forces and pressure related to atom j               *
		 * first neighbors of atom j                           *
		 ******************************************************/
		for ( pk = Cuda_Start_Index(j, bonds); pk < Cuda_End_Index(j, bonds); ++pk )
		{
			nbr_k = &bonds->select.bond_list[pk];
			k = nbr_k->nbr;

			rvec_MakeZero (nbr_k->tf_f);

			rvec_Scale( temp, -coef.C3dbo, nbr_k->bo_data.dBOp );      /*3rd,dBO*/
			rvec_ScaledAdd( temp, -coef.C3dDelta, nbr_k->bo_data.dBOp);/*dDelta*/
			rvec_ScaledAdd( temp, -coef.C4dbopi, nbr_k->bo_data.dBOp); /*4th,dBOpi*/
			rvec_ScaledAdd( temp, -coef.C4dbopi2, nbr_k->bo_data.dBOp);/*4th,dBOpi2*/

			/* force */
			rvec_Add( nbr_k->tf_f, temp );
			/* pressure */
			if ( k != i )
			{
				ivec_Sum( rel_box, nbr_k->rel_box, nbr_j->rel_box ); //rel_box(k, i)
				rvec_iMultiply( ext_press, rel_box, temp );
				rvec_Add( data_ext_press, ext_press );

				/* if( !ivec_isZero( rel_box ) )
                   fprintf( stderr, "%3d %3d %3d: dvec[%10.6f %10.6f %10.6f]"
                   "ext[%3d %3d %3d] f[%10.6f %10.6f %10.6f]\n",
                   i+1, j+1, system->my_atoms[j].x[0], 
                   system->my_atoms[j].x[1], system->my_atoms[j].x[2], 
                   k+1, system->my_atoms[k].x[0], 
                   system->my_atoms[k].x[1], system->my_atoms[k].x[2],
                   nbr_k->dvec[0], nbr_k->dvec[1], nbr_k->dvec[2],
                   rel_box[0], rel_box[1], rel_box[2],
                   temp[0], temp[1], temp[2] ); */
			}
		}

		/* then atom j itself */
		rvec_Scale( temp, -coef.C1dbo, bo_ij->dBOp );                    /*1st, dBO*/
		rvec_ScaledAdd( temp, coef.C3dbo, workspace->dDeltap_self[j] );  /*2nd, dBO*/
		rvec_ScaledAdd( temp, -coef.C1dDelta, bo_ij->dBOp );             /*1st, dBO*/
		rvec_ScaledAdd( temp, coef.C3dDelta, workspace->dDeltap_self[j]);/*2nd, dBO*/

		rvec_ScaledAdd( temp, -coef.C1dbopi, bo_ij->dln_BOp_pi );       /*1st,dBOpi*/
		rvec_ScaledAdd( temp, -coef.C2dbopi, bo_ij->dBOp );             /*2nd,dBOpi*/
		rvec_ScaledAdd( temp, coef.C4dbopi, workspace->dDeltap_self[j]);/*3rd,dBOpi*/

		rvec_ScaledAdd( temp, -coef.C1dbopi2, bo_ij->dln_BOp_pi2 );    /*1st,dBOpi2*/
		rvec_ScaledAdd( temp, -coef.C2dbopi2, bo_ij->dBOp );           /*2nd,dBOpi2*/
		rvec_ScaledAdd( temp,coef.C4dbopi2,workspace->dDeltap_self[j]);/*3rd,dBOpi2*/

		/* force */
		rvec_Add( workspace->f[j], temp );
		/* pressure */
		rvec_iMultiply( ext_press, nbr_j->rel_box, temp );
		rvec_Add( data->my_ext_press, ext_press );

		/* if( !ivec_isZero( nbr_j->rel_box ) )
           fprintf( stderr, "%3d %3d %3d: dvec[%10.6f %10.6f %10.6f]" 
           "ext[%3d %3d %3d] f[%10.6f %10.6f %10.6f]\n",
           i+1, system->my_atoms[i].x[0], system->my_atoms[i].x[1], 
           system->my_atoms[i].x[2], 
           j+1,system->my_atoms[j].x[0], system->my_atoms[j].x[1], 
           system->my_atoms[j].x[2],
           j+1, nbr_j->dvec[0], nbr_j->dvec[1], nbr_j->dvec[2],
           nbr_j->rel_box[0], nbr_j->rel_box[1], nbr_j->rel_box[2],
           temp[0], temp[1], temp[2] ); */
	}


}


CUDA_DEVICE void Cuda_Add_dBond_to_Forces( int i, int pj,
		storage *workspace, reax_list *bonds, reax_atom* my_atoms)
{
	bond_data *nbr_j, *nbr_k;
	bond_order_data *bo_ij, *bo_ji;
	dbond_coefficients coef;
	int pk, j;
	rvec tf_f;

	rvec_MakeZero( tf_f );

	/* Initializations */
	nbr_j = &bonds->select.bond_list[pj];
	j = nbr_j->nbr;

	//bo_ij = &nbr_j->bo_data;
	//bo_ji = &bonds->bond_list[ nbr_j->sym_index ].bo_data;

	if ( i < j )
	{
		bo_ij = &nbr_j->bo_data;
		bo_ji = &bonds->select.bond_list[ nbr_j->sym_index ].bo_data;
	}
	else
	{
		bo_ji = &nbr_j->bo_data;
		bo_ij = &bonds->select.bond_list[ nbr_j->sym_index ].bo_data;
	}


	coef.C1dbo = bo_ij->C1dbo * (bo_ij->Cdbo + bo_ji->Cdbo);
	coef.C2dbo = bo_ij->C2dbo * (bo_ij->Cdbo + bo_ji->Cdbo);
	coef.C3dbo = bo_ij->C3dbo * (bo_ij->Cdbo + bo_ji->Cdbo);

	coef.C1dbopi = bo_ij->C1dbopi * (bo_ij->Cdbopi + bo_ji->Cdbopi);
	coef.C2dbopi = bo_ij->C2dbopi * (bo_ij->Cdbopi + bo_ji->Cdbopi);
	coef.C3dbopi = bo_ij->C3dbopi * (bo_ij->Cdbopi + bo_ji->Cdbopi);
	coef.C4dbopi = bo_ij->C4dbopi * (bo_ij->Cdbopi + bo_ji->Cdbopi);

	coef.C1dbopi2 = bo_ij->C1dbopi2 * (bo_ij->Cdbopi2 + bo_ji->Cdbopi2);
	coef.C2dbopi2 = bo_ij->C2dbopi2 * (bo_ij->Cdbopi2 + bo_ji->Cdbopi2);
	coef.C3dbopi2 = bo_ij->C3dbopi2 * (bo_ij->Cdbopi2 + bo_ji->Cdbopi2);
	coef.C4dbopi2 = bo_ij->C4dbopi2 * (bo_ij->Cdbopi2 + bo_ji->Cdbopi2);

	coef.C1dDelta = bo_ij->C1dbo * (workspace->CdDelta[i]+workspace->CdDelta[j]);
	coef.C2dDelta = bo_ij->C2dbo * (workspace->CdDelta[i]+workspace->CdDelta[j]);
	coef.C3dDelta = bo_ij->C3dbo * (workspace->CdDelta[i]+workspace->CdDelta[j]);

	if ( i < j )
	{
		for ( pk = Cuda_Start_Index(i, bonds); pk < Cuda_End_Index(i, bonds); ++pk )
		{
			nbr_k = &bonds->select.bond_list[pk];

			//if(pk == 323 || pk == 134 || pk == 153 || pk == 174)
			//printf("Updating pk %d,i %d, j %d, atom i %d, atom j %d \n",pk, i,j,my_atoms[i].orig_id,my_atoms[j].orig_id);


			rvec_MakeZero( tf_f );

			/*2nd,dBO*/
			rvec_ScaledAdd( tf_f, -coef.C2dbo, nbr_k->bo_data.dBOp );
			/*dDelta*/
			rvec_ScaledAdd( tf_f, -coef.C2dDelta, nbr_k->bo_data.dBOp );
			/*3rd, dBOpi*/
			rvec_ScaledAdd( tf_f, -coef.C3dbopi, nbr_k->bo_data.dBOp );
			/*3rd, dBOpi2*/
			rvec_ScaledAdd( tf_f, -coef.C3dbopi2, nbr_k->bo_data.dBOp );

			//Temp storage
			rvec_Add( nbr_k->tf_f, tf_f );
		}



		/*1st, dBO*/
		rvec_ScaledAdd( workspace->f[i], coef.C1dbo, bo_ij->dBOp );
		/*2nd, dBO*/
		rvec_ScaledAdd( workspace->f[i], coef.C2dbo, workspace->dDeltap_self[i] );




		/*1st, dBO*/
		rvec_ScaledAdd( workspace->f[i], coef.C1dDelta, bo_ij->dBOp );
		/*2nd, dBO*/
		rvec_ScaledAdd( workspace->f[i], coef.C2dDelta, workspace->dDeltap_self[i] );

		/*1st, dBOpi*/
		rvec_ScaledAdd( workspace->f[i], coef.C1dbopi, bo_ij->dln_BOp_pi );
		/*2nd, dBOpi*/
		rvec_ScaledAdd( workspace->f[i], coef.C2dbopi, bo_ij->dBOp );
		/*3rd, dBOpi*/
		rvec_ScaledAdd( workspace->f[i], coef.C3dbopi, workspace->dDeltap_self[i] );

		/*1st, dBO_pi2*/
		rvec_ScaledAdd( workspace->f[i], coef.C1dbopi2, bo_ij->dln_BOp_pi2 );
		/*2nd, dBO_pi2*/
		rvec_ScaledAdd( workspace->f[i], coef.C2dbopi2, bo_ij->dBOp );
		/*3rd, dBO_pi2*/
		rvec_ScaledAdd( workspace->f[i], coef.C3dbopi2, workspace->dDeltap_self[i] );
	}
	else
	{
		for ( pk = Cuda_Start_Index(i, bonds); pk < Cuda_End_Index(i, bonds); ++pk )
		{
			nbr_k = &bonds->select.bond_list[pk];
			rvec_MakeZero( tf_f );

			//if(pk == 323 || pk == 134 || pk == 153 || pk == 174)
			//printf("Updating pk %d,i %d, j %d, atom i %d, atom j %d \n",pk, i,j,my_atoms[i].orig_id,my_atoms[j].orig_id);


			/*3rd, dBO*/
			rvec_ScaledAdd( tf_f, -coef.C3dbo, nbr_k->bo_data.dBOp );
			/*dDelta*/
			rvec_ScaledAdd( tf_f, -coef.C3dDelta, nbr_k->bo_data.dBOp );
			/*4th, dBOpi*/
			rvec_ScaledAdd( tf_f, -coef.C4dbopi, nbr_k->bo_data.dBOp );
			/*4th, dBOpi2*/
			rvec_ScaledAdd( tf_f, -coef.C4dbopi2, nbr_k->bo_data.dBOp );

			//Temp Storage
			rvec_Add( nbr_k->tf_f, tf_f );

		}

		/*1st,dBO*/
		rvec_ScaledAdd( workspace->f[i], -coef.C1dbo, bo_ij->dBOp );
		/*2nd,dBO*/
		rvec_ScaledAdd( workspace->f[i], coef.C3dbo, workspace->dDeltap_self[i] );

		/*1st, dBO*/
		rvec_ScaledAdd( workspace->f[i], -coef.C1dDelta, bo_ij->dBOp );
		/*2nd, dBO*/
		rvec_ScaledAdd( workspace->f[i], coef.C3dDelta, workspace->dDeltap_self[i] );



		/*1st, dBOpi*/
		rvec_ScaledAdd( workspace->f[i], -coef.C1dbopi, bo_ij->dln_BOp_pi );
		/*2nd, dBOpi*/
		rvec_ScaledAdd( workspace->f[i], -coef.C2dbopi, bo_ij->dBOp );
		/*3rd, dBOpi*/
		rvec_ScaledAdd( workspace->f[i], coef.C4dbopi, workspace->dDeltap_self[i] );

		/*1st, dBOpi2*/
		rvec_ScaledAdd( workspace->f[i], -coef.C1dbopi2, bo_ij->dln_BOp_pi2 );
		/*2nd, dBOpi2*/
		rvec_ScaledAdd( workspace->f[i], -coef.C2dbopi2, bo_ij->dBOp );
		/*3rd, dBOpi2*/
		rvec_ScaledAdd( workspace->f[i], coef.C4dbopi2, workspace->dDeltap_self[i] );

	}


}


CUDA_DEVICE void Cuda_dbond_to_Forces_postprocess( int i, reax_atom *atoms,
		reax_list *bonds, storage *workspace )
{
	int pk;
	bond_data *nbr_k, *nbr_k_sym;


	for( pk = Cuda_Start_Index(i, bonds); pk < Cuda_End_Index(i, bonds); ++pk )
	{

		nbr_k = &bonds->select.bond_list[pk];
		nbr_k_sym = &bonds->select.bond_list [nbr_k->sym_index];

		workspace->f[i][0] += nbr_k_sym->tf_f[0];
		workspace->f[i][1] += nbr_k_sym->tf_f[1];
		workspace->f[i][2] += nbr_k_sym->tf_f[2];



		//ret[1] += v[1];
		//ret[2] += v[2];


		//nbr_k_sym->tf_f[0] = nbr_k_sym->tf_f[0] + 1;
		//nbr_k_sym->tf_f[1] = nbr_k_sym->tf_f[1] + 1;
		//nbr_k_sym->tf_f[2] = nbr_k_sym->tf_f[2] + 1;

		//printf("%d,%d,%f,%f,%f\n",pk,nbr_k->sym_index,nbr_k_sym->tf_f[0],nbr_k_sym->tf_f[1],nbr_k_sym->tf_f[2]);

		//rvec_Add( atoms[i].f, nbr_k_sym->tf_f );
		//rvec_Add( workspace->f[i], nbr_k_sym->tf_f );
	}


}


CUDA_GLOBAL void k_total_forces_postprocess( reax_atom *my_atoms,
		reax_list p_bonds, storage p_workspace, int N )
{
	int i;
	reax_list *bonds;
	storage *workspace;

	i = blockIdx.x * blockDim.x + threadIdx.x;

	if ( i >= N )
	{
		return;
	}

	bonds = &p_bonds;
	workspace = &p_workspace;

	Cuda_dbond_to_Forces_postprocess( i, my_atoms, bonds, workspace );

	//if( i < 20)
		//printf("Forces %d,%f,%f,%f\n",my_atoms[i].orig_id, workspace->f[i][0],workspace->f[i][1],workspace->f[i][2]);

}


CUDA_GLOBAL void k_total_forces( storage p_workspace, reax_list p_bonds, 
		control_params *control, simulation_data *data, rvec *data_ext_press,
		int N, reax_atom *my_atoms )
{
	int i, pj;
	reax_list *bonds;
	storage *workspace;

	i = blockIdx.x * blockDim.x + threadIdx.x;

	if ( i >= N )
	{
		return;
	}


	bonds = &p_bonds;
	workspace = &p_workspace;

	for ( pj = Cuda_Start_Index(i, bonds); pj < Cuda_End_Index(i, bonds); ++pj )
	{
		//if ( i < bonds->bond_list[pj].nbr ) {
		if ( control->virial == 0 )
		{
			//if (i < 20 && my_atoms[i].orig_id == 13)
			//printf("%d,%d\n",pj);
			Cuda_Add_dBond_to_Forces( i, pj, workspace, bonds,my_atoms );

		}
		else
		{
			Cuda_Add_dBond_to_Forces_NPT( i, pj, data, workspace, bonds,
					data_ext_press[i] );
		}
	}

	//if( i < 20)
	//printf("Forces %d,%f,%f,%f\n",my_atoms[i].orig_id, workspace->f[i][0],workspace->f[i][1],workspace->f[i][2]);
}


void Cuda_Total_Forces( reax_system *system, control_params *control, 
		simulation_data *data, storage *workspace, reax_list **lists )
{
	int blocks;
	rvec *spad_rvec = (rvec *) workspace->scratch;

	cuda_memset( spad_rvec, 0, system->N * 2 * sizeof(rvec),
			"total_forces:ext_press" );

	blocks = system->N / DEF_BLOCK_SIZE
			+ ((system->N % DEF_BLOCK_SIZE == 0) ? 0 : 1);

	//printf("virital %d \n", control->virial);

	hipLaunchKernelGGL(k_total_forces, dim3(blocks), dim3(DEF_BLOCK_SIZE ), 0, 0,  *(workspace->d_workspace), *(lists[BONDS]),
			(control_params *) control->d_control_params,
			(simulation_data *)data->d_simulation_data,
			spad_rvec, system->N, system->d_my_atoms );
	hipDeviceSynchronize( );
	cudaCheckError( );



	if ( control->virial != 0 )
	{
		//do the reduction here for ext press
		hipLaunchKernelGGL(k_reduction_rvec, dim3(blocks), dim3(DEF_BLOCK_SIZE), sizeof(rvec) * DEF_BLOCK_SIZE , 0,  spad_rvec, spad_rvec + system->N, system->N );
		hipDeviceSynchronize( );
		cudaCheckError( );

		hipLaunchKernelGGL(k_reduction_rvec, dim3(1), dim3(control->blocks_pow_2_n), sizeof(rvec) * control->blocks_pow_2_n, 0,  spad_rvec + system->N, &((simulation_data *)data->d_simulation_data)->my_ext_press, blocks );
		hipDeviceSynchronize( );
		cudaCheckError( );
	}


	//do the post processing for the atomic forces here
	hipLaunchKernelGGL(k_total_forces_postprocess, dim3(blocks), dim3(DEF_BLOCK_SIZE ), 0, 0,  system->d_my_atoms, *(lists[BONDS]), *(workspace->d_workspace), system->N );
	hipDeviceSynchronize( );
	cudaCheckError( );

}


CUDA_GLOBAL void k_total_forces_pure( reax_atom *my_atoms, int n, 
		storage p_workspace )
{
	int i;
	storage *workspace;

	i = blockIdx.x * blockDim.x + threadIdx.x;

	if ( i >= n )
	{
		return;
	}

	workspace = &p_workspace;

	rvec_Copy( my_atoms[i].f, workspace->f[i] );
}


void Cuda_Total_Forces_PURE( reax_system *system, storage *workspace )
{
	int blocks;

	blocks = system->n / DEF_BLOCK_SIZE
			+ ((system->n % DEF_BLOCK_SIZE == 0) ? 0 : 1);

	hipLaunchKernelGGL(k_total_forces_pure, dim3(blocks), dim3(DEF_BLOCK_SIZE ), 0, 0,  system->d_my_atoms, system->n, *(workspace->d_workspace) );
	hipDeviceSynchronize( );
	cudaCheckError( );
}
