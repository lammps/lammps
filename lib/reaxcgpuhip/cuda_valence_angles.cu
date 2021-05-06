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

#include "cuda_valence_angles.h"

#include "cuda_list.h"

#include "index_utils.h"
#include "vector.h"


/* Compute 3-body interactions, in which the main role is played by
   atom j, which sits in the middle of the other two atoms i and k. */
CUDA_GLOBAL void Cuda_Valence_Angles( reax_atom *my_atoms,
		global_parameters gp, single_body_parameters *sbp, three_body_header *d_thbh,
		control_params *control, storage p_workspace, reax_list p_bonds,
		reax_list p_thb_intrs, int n, int N, int num_atom_types,
		real *data_e_ang, real *data_e_pen, real *data_e_coa, rvec *my_ext_press )
{
	int i, j, pi, k, pk, t;
	int type_i, type_j, type_k;
	int start_j, end_j;
	//    int start_pk, end_pk;
	int cnt, num_thb_intrs;
	real temp, temp_bo_jt, pBOjt7;
	real p_val1, p_val2, p_val3, p_val4, p_val5;
	real p_val6, p_val7, p_val8, p_val9, p_val10;
	real p_pen1, p_pen2, p_pen3, p_pen4;
	real p_coa1, p_coa2, p_coa3, p_coa4;
	real trm8, expval6, expval7, expval2theta, expval12theta, exp3ij, exp3jk;
	real exp_pen2ij, exp_pen2jk, exp_pen3, exp_pen4, trm_pen34, exp_coa2;
	real dSBO1, dSBO2, SBO, SBO2, CSBO2, SBOp, prod_SBO, vlpadj;
	real CEval1, CEval2, CEval3, CEval4, CEval5, CEval6, CEval7, CEval8;
	real CEpen1, CEpen2, CEpen3;
	real e_ang, e_coa, e_pen;
	real CEcoa1, CEcoa2, CEcoa3, CEcoa4, CEcoa5;
	real Cf7ij, Cf7jk, Cf8j, Cf9j;
	real f7_ij, f7_jk, f8_Dj, f9_Dj;
	real Ctheta_0, theta_0, theta_00, theta, cos_theta, sin_theta;
	real BOA_ij, BOA_jk;
	rvec force, ext_press;
	three_body_header *thbh;
	three_body_parameters *thbp;
	three_body_interaction_data *p_ijk;
	//    three_body_interaction_data *p_kji;
	bond_data *pbond_ij, *pbond_jk, *pbond_jt;
	bond_order_data *bo_ij, *bo_jk, *bo_jt;
	reax_list *bonds;
	reax_list *thb_intrs;
	storage *workspace;

	j = blockIdx.x * blockDim.x + threadIdx.x;

	if ( j >= N )
	{
		return;
	}



	bonds = &p_bonds;
	thb_intrs =  &p_thb_intrs;
	workspace = &p_workspace;
	/* global parameters used in these calculations */
	p_val6 = gp.l[14];
	p_val8 = gp.l[33];
	p_val9 = gp.l[16];
	p_val10 = gp.l[17];
	//num_thb_intrs = j * THREE_BODY_OFFSET;
	type_j = my_atoms[j].type;
	start_j = Cuda_Start_Index( j, bonds );
	end_j = Cuda_End_Index( j, bonds );
	p_val3 = sbp[ type_j ].p_val3;
	p_val5 = sbp[ type_j ].p_val5;
	SBOp = 0.0;
	prod_SBO = 1.0;




	for( t = start_j; t < end_j; ++t )
	{
		bo_jt = &bonds->select.bond_list[t].bo_data;
		SBOp += (bo_jt->BO_pi + bo_jt->BO_pi2);
		temp = SQR( bo_jt->BO );
		temp *= temp;
		temp *= temp;
		prod_SBO *= EXP( -temp );
	}


	/* modifications to match Adri's code - 09/01/09 */
	if( workspace->vlpex[j] >= 0 )
	{

		vlpadj = 0;
		dSBO2 = prod_SBO - 1;
	}
	else
	{
		vlpadj = workspace->nlp[j];
		dSBO2 = (prod_SBO - 1.0) * (1.0 - p_val8 * workspace->dDelta_lp[j]);
	}



	SBO = SBOp + (1 - prod_SBO) * (-workspace->Delta_boc[j] - p_val8 * vlpadj);
	dSBO1 = -8 * prod_SBO * ( workspace->Delta_boc[j] + p_val8 * vlpadj );

	if( SBO <= 0 )
	{
		SBO2 = 0;
		CSBO2 = 0;
	}
	else if( SBO > 0 && SBO <= 1 )
	{
		SBO2 = POW( SBO, p_val9 );
		CSBO2 = p_val9 * POW( SBO, p_val9 - 1 );
	}
	else if( SBO > 1 && SBO < 2 )
	{
		SBO2 = 2 - POW( 2-SBO, p_val9 );
		CSBO2 = p_val9 * POW( 2 - SBO, p_val9 - 1 );
	}
	else
	{
		SBO2 = 2;
		CSBO2 = 0;
	}

	expval6 = EXP( p_val6 * workspace->Delta_boc[j] );

	for( pi = start_j; pi < end_j; ++pi )
	{
		num_thb_intrs = Cuda_Start_Index( pi, thb_intrs );

		pbond_ij = &bonds->select.bond_list[pi];
		bo_ij = &pbond_ij->bo_data;
		BOA_ij = bo_ij->BO - control->thb_cut;

		if ( BOA_ij > 0.0
				&& ( j < n || pbond_ij->nbr < n ) )
		{
			i = pbond_ij->nbr;


			type_i = my_atoms[i].type;

			/* first copy 3-body intrs from previously computed ones where i > k;
               in the second for-loop below, compute only new 3-body intrs where i < k */

			// The copy loop commented out because strange asynchronous issues started to surface
			// Each kernel now manually generates everything
			//            for( pk = start_j; pk < pi; ++pk )
			//            {
			//                start_pk = Cuda_Start_Index( pk, thb_intrs );
			//                end_pk = Cuda_End_Index( pk, thb_intrs );
			//
			//                for( t = start_pk; t < end_pk; ++t )
			//                {
			//                    if( thb_intrs->three_body_list[t].thb == i )
			//                    {
			//                        p_ijk = &thb_intrs->three_body_list[num_thb_intrs];
			//                        p_kji = &thb_intrs->three_body_list[t];
			//
			//                        p_ijk->thb = bonds->bond_list[pk].nbr;
			//                        p_ijk->pthb  = pk;
			//                        p_ijk->theta = p_kji->theta;
			//                        rvec_Copy( p_ijk->dcos_di, p_kji->dcos_dk );
			//                        rvec_Copy( p_ijk->dcos_dj, p_kji->dcos_dj );
			//                        rvec_Copy( p_ijk->dcos_dk, p_kji->dcos_di );
			//
			//                        ++num_thb_intrs;
			//                        break;
			//                    }
			//                }
			//            }

			/* and this is the second for loop mentioned above */
			//for( pk = pi+1; pk < end_j; ++pk ) {
			// Except that now the loop goes all the way from start_j to end_j
			for( pk = start_j; pk < end_j; ++pk )
			{
				if ( pk == pi )
				{
					continue;
				}

				pbond_jk = &bonds->select.bond_list[pk];
				bo_jk = &pbond_jk->bo_data;
				BOA_jk = bo_jk->BO - control->thb_cut;
				k = pbond_jk->nbr;
				type_k = my_atoms[k].type;
				p_ijk = &thb_intrs->select.three_body_list[num_thb_intrs];




				//CHANGE ORIGINAL
				//if ((BOA_jk <= 0) || ((j >= n) && (k >= n))) continue;
				if ( BOA_jk <= 0.0 )
				{

					continue;
				}
				//CHANGE ORIGINAL
				if(j == 0) {
					//printf("%d,%d,%d\n",i,j,k);

				}



				Calculate_Theta( pbond_ij->dvec, pbond_ij->d,
						pbond_jk->dvec, pbond_jk->d,
						&theta, &cos_theta );


				//printf("%d,%d,%f\n",i,j,theta);

				//printf("%f,%f,%f,%f,%f \n",pbond_ij->dvec[0],pbond_ij->dvec[1],pbond_ij->dvec[2],pbond_ij->d,theta);


				Calculate_dCos_Theta( pbond_ij->dvec, pbond_ij->d,
						pbond_jk->dvec, pbond_jk->d,
						&p_ijk->dcos_di, &p_ijk->dcos_dj,
						&p_ijk->dcos_dk );


				if(j == 0)
				{
					//printf("calc dos di %d,%d,%f,%f,%f,%f,%f,%f,%f,%f\n",pi,pk,pbond_ij->dvec[0],pbond_ij->dvec[1],pbond_ij->dvec[2],pbond_ij->d,pbond_jk->dvec[0],pbond_jk->dvec[1],pbond_jk->dvec[2],pbond_jk->d);
				}

				p_ijk->thb = k;
				p_ijk->pthb = pk;
				p_ijk->theta = theta;
				sin_theta = SIN( theta );

				if ( sin_theta < 1.0e-5 )
				{
					sin_theta = 1.0e-5;
				}

				++num_thb_intrs;

				if ( j < n && BOA_jk > 0.0
						&& bo_ij->BO * bo_jk->BO > SQR(control->thb_cut) )
				{
					thbh = &d_thbh[ index_thbp(type_i, type_j, type_k, num_atom_types) ];

					if(j == 0) {
						//printf("%d,%d,%d,%d\n",i,j,k,thbh->cnt);

					}

					for ( cnt = 0; cnt < thbh->cnt; ++cnt )
					{
						if ( FABS(thbh->prm[cnt].p_val1) > 0.001 )
						{
							thbp = &thbh->prm[cnt];

							/* ANGLE ENERGY */
							p_val1 = thbp->p_val1;
							p_val2 = thbp->p_val2;
							p_val4 = thbp->p_val4;
							p_val7 = thbp->p_val7;
							theta_00 = thbp->theta_00;


							//printf("Valence %f,%f,%f,%f,%f \n",p_val1,p_val2,p_val4,p_val7,theta_00);


							exp3ij = EXP( -p_val3 * POW( BOA_ij, p_val4 ) );
							f7_ij = 1.0 - exp3ij;
							Cf7ij = p_val3 * p_val4 * POW( BOA_ij, p_val4 - 1.0 ) * exp3ij;

							exp3jk = EXP( -p_val3 * POW( BOA_jk, p_val4 ) );
							f7_jk = 1.0 - exp3jk;
							Cf7jk = p_val3 * p_val4 * POW( BOA_jk, p_val4 - 1.0 ) * exp3jk;

							expval7 = EXP( -p_val7 * workspace->Delta_boc[j] );
							trm8 = 1.0 + expval6 + expval7;
							f8_Dj = p_val5 - ( (p_val5 - 1.0) * (2.0 + expval6) / trm8 );
							Cf8j = ( (1.0 - p_val5) / SQR(trm8) ) *
									( p_val6 * expval6 * trm8 -
											(2.0 + expval6) * ( p_val6*expval6 - p_val7*expval7 ) );

							theta_0 = 180.0 - theta_00 * (1.0 -
									EXP(-p_val10 * (2.0 - SBO2)));
							theta_0 = DEG2RAD( theta_0 );

							expval2theta  = EXP( -p_val2 * SQR(theta_0 - theta) );

							if ( p_val1 >= 0 )
							{
								expval12theta = p_val1 * (1.0 - expval2theta);
							}
							/* to avoid linear Me-H-Me angles (6/6/06) */
							else
							{
								expval12theta = p_val1 * -expval2theta;
							}

							CEval1 = Cf7ij * f7_jk * f8_Dj * expval12theta;


							CEval2 = Cf7jk * f7_ij * f8_Dj * expval12theta;
							CEval3 = Cf8j  * f7_ij * f7_jk * expval12theta;





							CEval4 = -2.0 * p_val1 * p_val2 * f7_ij * f7_jk * f8_Dj *
									expval2theta * (theta_0 - theta);

							Ctheta_0 = p_val10 * DEG2RAD(theta_00) *
									exp( -p_val10 * (2.0 - SBO2) );

							CEval5 = -CEval4 * Ctheta_0 * CSBO2;
							CEval6 = CEval5 * dSBO1;
							CEval7 = CEval5 * dSBO2;
							CEval8 = -CEval4 / sin_theta;



							if ( pk < pi )
							{
								e_ang = f7_ij * f7_jk * f8_Dj * expval12theta;
								data_e_ang[j] += e_ang;

							}
							/* END ANGLE ENERGY*/

							/* PENALTY ENERGY */
							p_pen1 = thbp->p_pen1;
							p_pen2 = gp.l[19];
							p_pen3 = gp.l[20];
							p_pen4 = gp.l[21];

							exp_pen2ij = EXP( -p_pen2 * SQR( BOA_ij - 2.0 ) );
							exp_pen2jk = EXP( -p_pen2 * SQR( BOA_jk - 2.0 ) );
							exp_pen3 = EXP( -p_pen3 * workspace->Delta[j] );
							exp_pen4 = EXP(  p_pen4 * workspace->Delta[j] );
							trm_pen34 = 1.0 + exp_pen3 + exp_pen4;
							f9_Dj = ( 2.0 + exp_pen3 ) / trm_pen34;
							Cf9j = ( -p_pen3 * exp_pen3 * trm_pen34 - (2.0 + exp_pen3)
									* ( -p_pen3 * exp_pen3 + p_pen4 * exp_pen4 ) )
                                																																												/ SQR( trm_pen34 );

							/* very important: since each kernel generates all interactions,
                               need to prevent all energies becoming duplicates */
							if ( pk < pi )
							{
								e_pen = p_pen1 * f9_Dj * exp_pen2ij * exp_pen2jk;
								data_e_pen[j] += e_pen;
							}

							CEpen1 = e_pen * Cf9j / f9_Dj;
							temp = -2.0 * p_pen2 * e_pen;
							CEpen2 = temp * (BOA_ij - 2.0);
							CEpen3 = temp * (BOA_jk - 2.0);
							/* END PENALTY ENERGY */

							/* COALITION ENERGY */
							p_coa1 = thbp->p_coa1;
							p_coa2 = gp.l[2];
							p_coa3 = gp.l[38];
							p_coa4 = gp.l[30];

							exp_coa2 = EXP( p_coa2 * workspace->Delta_boc[j] );

							/* similar to above comment regarding if statement */
							if ( pk < pi )
							{
								e_coa =
										p_coa1 / (1. + exp_coa2) *
										EXP( -p_coa3 * SQR(workspace->total_bond_order[i] - BOA_ij) ) *
										EXP( -p_coa3 * SQR(workspace->total_bond_order[k] - BOA_jk) ) *
										EXP( -p_coa4 * SQR(BOA_ij - 1.5) ) *
										EXP( -p_coa4 * SQR(BOA_jk - 1.5) );
								data_e_coa[j] += e_coa;
							}

							CEcoa1 = -2 * p_coa4 * (BOA_ij - 1.5) * e_coa;
							CEcoa2 = -2 * p_coa4 * (BOA_jk - 1.5) * e_coa;
							CEcoa3 = -p_coa2 * exp_coa2 * e_coa / (1 + exp_coa2);
							CEcoa4 = -2 * p_coa3 *
									(workspace->total_bond_order[i] - BOA_ij) * e_coa;
							CEcoa5 = -2 * p_coa3 *
									(workspace->total_bond_order[k] - BOA_jk) * e_coa;
							/* END COALITION ENERGY */

							/* FORCES */
							// we must again check for pk<pi for entire forces part
							if ( pk < pi )
							{
								//printf("Valence %f,%f,%f,%f \n",CEval1,CEpen2,CEcoa1,CEcoa4);

								bo_ij->Cdbo += (CEval1 + CEpen2 + (CEcoa1 - CEcoa4));



								bo_jk->Cdbo += (CEval2 + CEpen3 + (CEcoa2 - CEcoa5));

								workspace->CdDelta[j] += ((CEval3 + CEval7) + CEpen1 + CEcoa3);


								//                                workspace->CdDelta[i] += CEcoa4;
								//                                workspace->CdDelta[k] += CEcoa5;
								pbond_ij->va_CdDelta += CEcoa4;
								pbond_jk->va_CdDelta += CEcoa5;

								for ( t = start_j; t < end_j; ++t )
								{
									pbond_jt = &bonds->select.bond_list[t];
									bo_jt = &pbond_jt->bo_data;
									temp_bo_jt = bo_jt->BO;
									temp = CUBE( temp_bo_jt );
									pBOjt7 = temp * temp * temp_bo_jt;

									bo_jt->Cdbo += (CEval6 * pBOjt7);
									bo_jt->Cdbopi += CEval5;
									bo_jt->Cdbopi2 += CEval5;
								}

								if ( control->virial == 0 )
								{
									rvec_ScaledAdd( pbond_ij->va_f, CEval8, p_ijk->dcos_di);
									rvec_ScaledAdd( workspace->f[j], CEval8, p_ijk->dcos_dj );
									rvec_ScaledAdd( pbond_jk->va_f, CEval8, p_ijk->dcos_dk );



								}
								else
								{
									/* terms not related to bond order derivatives are
                                       added directly into forces and pressure vector/tensor */
									rvec_Scale( force, CEval8, p_ijk->dcos_di );
									//                                    rvec_Add( workspace->f[i], force );
									rvec_Add( pbond_ij->va_f, force );
									rvec_iMultiply( ext_press, pbond_ij->rel_box, force );
									//                                    rvec_Add( data->my_ext_press, ext_press );
									rvec_Add( my_ext_press[j], ext_press );

									rvec_ScaledAdd( workspace->f[j], CEval8, p_ijk->dcos_dj );

									rvec_Scale( force, CEval8, p_ijk->dcos_dk );
									//                                    rvec_Add( workspace->f[k], force );
									rvec_Add( pbond_jk->va_f, force );
									rvec_iMultiply( ext_press, pbond_jk->rel_box, force );
									rvec_Add( my_ext_press[j], ext_press );
								}
							}

#ifdef TEST_ENERGY
							/*fprintf( out_control->eval, "%12.8f%12.8f%12.8f%12.8f\n",
                              p_val3, p_val4, BOA_ij, BOA_jk );
                              fprintf(out_control->eval, "%13.8f%13.8f%13.8f%13.8f%13.8f\n",
                              workspace->Delta_e[j], workspace->vlpex[j],
                              dSBO1, dSBO2, vlpadj );
                              fprintf( out_control->eval, "%12.8f%12.8f%12.8f%12.8f\n",
                              f7_ij, f7_jk, f8_Dj, expval12theta );
                              fprintf( out_control->eval,
                              "%12.8f%12.8f%12.8f%12.8f%12.8f%12.8f%12.8f%12.8f\n",
                              CEval1, CEval2, CEval3, CEval4,
                              CEval5, CEval6, CEval7, CEval8 );

                              fprintf( out_control->eval,
                              "%12.8f%12.8f%12.8f\n%12.8f%12.8f%12.8f\n%12.8f%12.8f%12.8f\n",
                              p_ijk->dcos_di[0]/sin_theta, p_ijk->dcos_di[1]/sin_theta,
                              p_ijk->dcos_di[2]/sin_theta,
                              p_ijk->dcos_dj[0]/sin_theta, p_ijk->dcos_dj[1]/sin_theta,
                              p_ijk->dcos_dj[2]/sin_theta,
                              p_ijk->dcos_dk[0]/sin_theta, p_ijk->dcos_dk[1]/sin_theta,
                              p_ijk->dcos_dk[2]/sin_theta);

                              fprintf( out_control->eval,
                              "%6d%6d%6d%15.8f%15.8f\n",
                              system->my_atoms[i].orig_id,
                              system->my_atoms[j].orig_id,
                              system->my_atoms[k].orig_id,
                              RAD2DEG(theta), e_ang );*/

							fprintf( out_control->eval,
									//"%6d%6d%6d%24.15e%24.15e%24.15e%24.15e%24.15e%24.15e\n",
									"%6d%6d%6d%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f\n",
									system->my_atoms[i].orig_id,
									system->my_atoms[j].orig_id,
									system->my_atoms[k].orig_id,
									RAD2DEG(theta), theta_0, BOA_ij, BOA_jk,
									e_ang, data->my_en.e_ang );

							fprintf( out_control->epen,
									//"%6d%6d%6d%24.15e%24.15e%24.15e%24.15e%24.15e\n",
									"%6d%6d%6d%12.4f%12.4f%12.4f%12.4f%12.4f\n",
									system->my_atoms[i].orig_id,
									system->my_atoms[j].orig_id,
									system->my_atoms[k].orig_id,
									RAD2DEG(theta), BOA_ij, BOA_jk, e_pen,
									data->my_en.e_pen );

							fprintf( out_control->ecoa,
									//"%6d%6d%6d%24.15e%24.15e%24.15e%24.15e%24.15e\n",
									"%6d%6d%6d%12.4f%12.4f%12.4f%12.4f%12.4f\n",
									system->my_atoms[i].orig_id,
									system->my_atoms[j].orig_id,
									system->my_atoms[k].orig_id,
									RAD2DEG(theta), BOA_ij, BOA_jk,
									e_coa, data->my_en.e_coa );
#endif

#ifdef TEST_FORCES
							/* angle forces */
							Add_dBO( system, lists, j, pi, CEval1, workspace->f_ang );
							Add_dBO( system, lists, j, pk, CEval2, workspace->f_ang );
							Add_dDelta( system, lists, j,
									CEval3 + CEval7, workspace->f_ang );

							for( t = start_j; t < end_j; ++t )
							{
								pbond_jt = &bonds->bond_list[t];
								bo_jt = &pbond_jt->bo_data;
								temp_bo_jt = bo_jt->BO;
								temp = CUBE( temp_bo_jt );
								pBOjt7 = temp * temp * temp_bo_jt;

								Add_dBO( system, lists, j, t, pBOjt7 * CEval6,
										workspace->f_ang );
								Add_dBOpinpi2( system, lists, j, t, CEval5, CEval5,
										workspace->f_ang, workspace->f_ang );
							}

							rvec_ScaledAdd( workspace->f_ang[i], CEval8, p_ijk->dcos_di );
							rvec_ScaledAdd( workspace->f_ang[j], CEval8, p_ijk->dcos_dj );
							rvec_ScaledAdd( workspace->f_ang[k], CEval8, p_ijk->dcos_dk );
							/* end angle forces */

							/* penalty forces */
							Add_dDelta( system, lists, j, CEpen1, workspace->f_pen );
							Add_dBO( system, lists, j, pi, CEpen2, workspace->f_pen );
							Add_dBO( system, lists, j, pk, CEpen3, workspace->f_pen );
							/* end penalty forces */

							/* coalition forces */
							Add_dBO( system, lists, j, pi, CEcoa1 - CEcoa4,
									workspace->f_coa );
							Add_dBO( system, lists, j, pk, CEcoa2 - CEcoa5,
									workspace->f_coa );
							Add_dDelta( system, lists, j, CEcoa3, workspace->f_coa );
							Add_dDelta( system, lists, i, CEcoa4, workspace->f_coa );
							Add_dDelta( system, lists, k, CEcoa5, workspace->f_coa );
							/* end coalition forces */
#endif
						}
					}
				}
			}
		}
		Cuda_Set_End_Index( pi, num_thb_intrs, thb_intrs );
	}
}


CUDA_GLOBAL void Cuda_Valence_Angles_PostProcess( reax_atom *atoms,
		control_params *control, storage p_workspace,
		reax_list p_bonds, int N )
{
	int i, pj;
	bond_data *pbond;
	bond_data *sym_index_bond;
	reax_list *bonds;
	storage *workspace;

	i = blockIdx.x * blockDim.x + threadIdx.x;

	if ( i >= N )
	{
		return;
	}

	bonds = &p_bonds;
	workspace = &p_workspace;

	for( pj = Cuda_Start_Index(i, bonds); pj < Cuda_End_Index(i, bonds); ++pj )
	{
		pbond = &bonds->select.bond_list[pj];

		sym_index_bond = &bonds->select.bond_list[ pbond->sym_index ];
		workspace->CdDelta[i] += sym_index_bond->va_CdDelta;

		//rvec_Add( atoms[i].f, sym_index_bond->va_f );
		rvec_Add( workspace->f[i], sym_index_bond->va_f);
	}
}


/* Estimate the num. of three-body interactions */
CUDA_GLOBAL void Estimate_Cuda_Valence_Angles( reax_atom *my_atoms,
		control_params *control, reax_list p_bonds, int n, int N, int *count )
{
	int j, pi, pk;
	int start_j, end_j;
	int num_thb_intrs;
	real BOA_ij, BOA_jk;
	bond_data *pbond_ij, *pbond_jk;
	bond_order_data *bo_ij, *bo_jk;
	reax_list *bonds;

	j = blockIdx.x * blockDim.x + threadIdx.x;

	if ( j >= N )
	{
		return;
	}

	bonds = &p_bonds;
	start_j = Cuda_Start_Index( j, bonds );
	end_j = Cuda_End_Index( j, bonds );

	for ( pi = start_j; pi < end_j; ++pi )
	{
		num_thb_intrs = 0;
		count[ pi ] = 0;

		pbond_ij = &bonds->select.bond_list[pi];
		bo_ij = &pbond_ij->bo_data;
		BOA_ij = bo_ij->BO - control->thb_cut;





		if ( BOA_ij > 0.0 && ( j < n || pbond_ij->nbr < n ) )
		{
			for ( pk = start_j; pk < end_j; ++pk )
			{
				if ( pk == pi )
				{
					continue;
				}

				pbond_jk = &bonds->select.bond_list[pk];
				bo_jk = &pbond_jk->bo_data;
				BOA_jk = bo_jk->BO - control->thb_cut;



				//CHANGE ORIGINAL
				//if ((BOA_jk <= 0) || ((j >= n) && (k >= n))) continue;
				if ( BOA_jk <= 0.0 )
				{
					continue;
				}
				//CHANGE ORIGINAL

				//printf("start %d, end %d,%f,%f \n", start_j,end_j,bo_jk->BO,BOA_jk);


				++num_thb_intrs;
			}

		}

		count[ pi ] = num_thb_intrs;

	}
}
