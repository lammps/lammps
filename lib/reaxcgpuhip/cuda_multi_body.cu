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

#include "cuda_multi_body.h"

#include "cuda_helpers.h"
#include "cuda_list.h"

#include "index_utils.h"


CUDA_GLOBAL void Cuda_Atom_Energy( reax_atom *my_atoms, global_parameters gp, 
		single_body_parameters *sbp, two_body_parameters *tbp,
		storage p_workspace, reax_list p_bonds, int n, int num_atom_types,
		real *data_elp, real *data_eov, real *data_eun )
{
	int i, j, pj, type_i, type_j;
	real Delta_lpcorr, dfvl;
	real e_lp, expvd2, inv_expvd2, dElp, CElp, DlpVi;
	real e_lph, Di, vov3, deahu2dbo, deahu2dsbo;
	real e_ov, CEover1, CEover2, CEover3, CEover4;
	real exp_ovun1, exp_ovun2, sum_ovun1, sum_ovun2;
	real exp_ovun2n, exp_ovun6, exp_ovun8;
	real inv_exp_ovun1, inv_exp_ovun2, inv_exp_ovun2n, inv_exp_ovun8;
	real e_un, CEunder1, CEunder2, CEunder3, CEunder4;
	real p_lp2, p_lp3;
	real p_ovun2, p_ovun3, p_ovun4, p_ovun5, p_ovun6, p_ovun7, p_ovun8;
	single_body_parameters *sbp_i;
	two_body_parameters *twbp;
	bond_data *pbond;
	bond_order_data *bo_ij;

	i = blockIdx.x * blockDim.x + threadIdx.x;

	if ( i >= n )
	{
		return;
	}

	reax_list *bonds = &p_bonds;
	storage *workspace = &p_workspace;

	/* Initialize parameters */
	p_lp3 = gp.l[5];
	p_ovun3 = gp.l[32];
	p_ovun4 = gp.l[31];
	p_ovun6 = gp.l[6];
	p_ovun7 = gp.l[8];
	p_ovun8 = gp.l[9];

	//for( i = 0; i < system->n; ++i ) {
	/* set the parameter pointer */
	type_i = my_atoms[i].type;
	sbp_i = &sbp[ type_i ];

	/* lone-pair Energy */
	p_lp2 = sbp_i->p_lp2;
	expvd2 = EXP( -75 * workspace->Delta_lp[i] );
	inv_expvd2 = 1. / (1. + expvd2 );

	/* calculate the energy */
	e_lp = p_lp2 * workspace->Delta_lp[i] * inv_expvd2;
	data_elp[i] += e_lp;

	dElp = p_lp2 * inv_expvd2 +
			75 * p_lp2 * workspace->Delta_lp[i] * expvd2 * SQR(inv_expvd2);
	CElp = dElp * workspace->dDelta_lp[i];

	workspace->CdDelta[i] += CElp;  // lp - 1st term

#ifdef TEST_ENERGY
	//  fprintf( out_control->elp, "%24.15e%24.15e%24.15e%24.15e\n",
	//         p_lp2, workspace->Delta_lp_temp[i], expvd2, dElp );
	//  fprintf( out_control->elp, "%6d%24.15e%24.15e%24.15e\n",
	fprintf( out_control->elp, "%6d%12.4f%12.4f%12.4f\n",
			system->my_atoms[i].orig_id, workspace->nlp[i],
			e_lp, data->my_en.e_lp );
#endif

#ifdef TEST_FORCES
	Add_dDelta( system, lists, i, CElp, workspace->f_lp );  // lp - 1st term
#endif

	/* correction for C2 */
	if ( gp.l[5] > 0.001
			&& !cuda_strcmp( sbp[type_i].name, "C", 1 ) )
	{
		for ( pj = Cuda_Start_Index(i, bonds); pj < Cuda_End_Index(i, bonds); ++pj )
		{
			if ( my_atoms[i].orig_id <
					my_atoms[bonds->select.bond_list[pj].nbr].orig_id )
			{
				j = bonds->select.bond_list[pj].nbr;
				type_j = my_atoms[j].type;

				if ( !cuda_strcmp( sbp[type_j].name, "C", 1 ) )
				{
					twbp = &tbp[index_tbp (type_i,type_j, num_atom_types) ];
					bo_ij = &bonds->select.bond_list[pj].bo_data;

					Di = workspace->Delta[i];
					vov3 = bo_ij->BO - Di - 0.040 * POW(Di, 4.);

					if ( vov3 > 3. )
					{
						e_lph = p_lp3 * SQR( vov3 - 3.0 );
						data_elp[i] += e_lph;



						deahu2dbo = 2. * p_lp3 * (vov3 - 3.);
						deahu2dsbo = 2. * p_lp3 * (vov3 - 3.) *
								(-1. - 0.16 * POW(Di, 3.));

						bo_ij->Cdbo += deahu2dbo;

						//printf("cdbo %f\n",bo_ij->Cdbo);

						workspace->CdDelta[i] += deahu2dsbo;

#ifdef TEST_ENERGY
						fprintf(out_control->elp,"C2cor%6d%6d%12.6f%12.6f%12.6f\n",
								system->my_atoms[i].orig_id, system->my_atoms[j].orig_id,
								e_lph, deahu2dbo, deahu2dsbo );
#endif

#ifdef TEST_FORCES
						Add_dBO(system, lists, i, pj, deahu2dbo, workspace->f_lp);
						Add_dDelta(system, lists, i, deahu2dsbo, workspace->f_lp);
#endif
					}
				}
			}
		}
	}
	//}

	//for( i = 0; i < system->n; ++i ) {
	type_i = my_atoms[i].type;
	sbp_i = &sbp[ type_i ];

	/* over-coordination energy */
	if( sbp_i->mass > 21.0 )
	{
		dfvl = 0.0;
	}
	else
	{
		dfvl = 1.0; // only for 1st-row elements
	}

	p_ovun2 = sbp_i->p_ovun2;
	sum_ovun1 = sum_ovun2 = 0;
	for( pj = Cuda_Start_Index(i, bonds); pj < Cuda_End_Index(i, bonds); ++pj )
	{
		j = bonds->select.bond_list[pj].nbr;
		type_j = my_atoms[j].type;
		bo_ij = &bonds->select.bond_list[pj].bo_data;


		twbp = &tbp[ index_tbp(type_i, type_j, num_atom_types )];

		sum_ovun1 += twbp->p_ovun1 * twbp->De_s * bo_ij->BO;
		sum_ovun2 += (workspace->Delta[j] - dfvl*workspace->Delta_lp_temp[j])*
				( bo_ij->BO_pi + bo_ij->BO_pi2 );
	}

	exp_ovun1 = p_ovun3 * EXP( p_ovun4 * sum_ovun2 );
	inv_exp_ovun1 = 1.0 / (1 + exp_ovun1);
	Delta_lpcorr  = workspace->Delta[i] -
			(dfvl * workspace->Delta_lp_temp[i]) * inv_exp_ovun1;

	exp_ovun2 = EXP( p_ovun2 * Delta_lpcorr );
	inv_exp_ovun2 = 1.0 / (1.0 + exp_ovun2);

	DlpVi = 1.0 / (Delta_lpcorr + sbp_i->valency + 1e-8);
	CEover1 = Delta_lpcorr * DlpVi * inv_exp_ovun2;

	e_ov = sum_ovun1 * CEover1;
	data_eov[i] += e_ov;

	CEover2 = sum_ovun1 * DlpVi * inv_exp_ovun2 *
			(1.0 - Delta_lpcorr * ( DlpVi + p_ovun2 * exp_ovun2 * inv_exp_ovun2 ));

	CEover3 = CEover2 * (1.0 - dfvl * workspace->dDelta_lp[i] * inv_exp_ovun1 );

	CEover4 = CEover2 * (dfvl * workspace->Delta_lp_temp[i]) *
			p_ovun4 * exp_ovun1 * SQR(inv_exp_ovun1);

	/* under-coordination potential */
	p_ovun2 = sbp_i->p_ovun2;
	p_ovun5 = sbp_i->p_ovun5;

	exp_ovun2n = 1.0 / exp_ovun2;
	exp_ovun6 = EXP( p_ovun6 * Delta_lpcorr );
	exp_ovun8 = p_ovun7 * EXP(p_ovun8 * sum_ovun2);
	inv_exp_ovun2n = 1.0 / (1.0 + exp_ovun2n);
	inv_exp_ovun8 = 1.0 / (1.0 + exp_ovun8);

	e_un = -p_ovun5 * (1.0 - exp_ovun6) * inv_exp_ovun2n * inv_exp_ovun8;
	data_eun[i] += e_un;

	CEunder1 = inv_exp_ovun2n *
			( p_ovun5 * p_ovun6 * exp_ovun6 * inv_exp_ovun8 +
					p_ovun2 * e_un * exp_ovun2n );
	CEunder2 = -e_un * p_ovun8 * exp_ovun8 * inv_exp_ovun8;
	CEunder3 = CEunder1 * (1.0 - dfvl*workspace->dDelta_lp[i]*inv_exp_ovun1);
	CEunder4 = CEunder1 * (dfvl*workspace->Delta_lp_temp[i]) *
			p_ovun4 * exp_ovun1 * SQR(inv_exp_ovun1) + CEunder2;

	/* forces */
	workspace->CdDelta[i] += CEover3;   // OvCoor - 2nd term
	workspace->CdDelta[i] += CEunder3;  // UnCoor - 1st term

#ifdef TEST_FORCES
	Add_dDelta( system, lists, i, CEover3, workspace->f_ov ); // OvCoor 2nd
	Add_dDelta( system, lists, i, CEunder3, workspace->f_un ); // UnCoor 1st
#endif

	for( pj = Cuda_Start_Index(i, bonds); pj < Cuda_End_Index(i, bonds); ++pj )
	{
		pbond = &bonds->select.bond_list[pj];
		j = pbond->nbr;
		bo_ij = &pbond->bo_data;
		twbp  = &tbp[ index_tbp(my_atoms[i].type, my_atoms[pbond->nbr].type,
				num_atom_types) ];

		bo_ij->Cdbo += CEover1 * twbp->p_ovun1 * twbp->De_s;// OvCoor-1st


		//workspace->CdDelta[j] += CEover4 * (1.0 - dfvl*workspace->dDelta_lp[j]) *
		pbond->ae_CdDelta += CEover4 * (1.0 - dfvl*workspace->dDelta_lp[j]) *
				(bo_ij->BO_pi + bo_ij->BO_pi2); // OvCoor-3a
		bo_ij->Cdbopi += CEover4 *
				(workspace->Delta[j] - dfvl*workspace->Delta_lp_temp[j]); // OvCoor-3b
		bo_ij->Cdbopi2 += CEover4 *
				(workspace->Delta[j] - dfvl*workspace->Delta_lp_temp[j]);  // OvCoor-3b

		//workspace->CdDelta[j] += CEunder4 * (1.0 - dfvl*workspace->dDelta_lp[j]) *
		pbond->ae_CdDelta += CEunder4 * (1.0 - dfvl*workspace->dDelta_lp[j]) *
				(bo_ij->BO_pi + bo_ij->BO_pi2);   // UnCoor - 2a
		bo_ij->Cdbopi += CEunder4 *
				(workspace->Delta[j] - dfvl*workspace->Delta_lp_temp[j]);  // UnCoor-2b
		bo_ij->Cdbopi2 += CEunder4 *
				(workspace->Delta[j] - dfvl*workspace->Delta_lp_temp[j]);  // UnCoor-2b

#ifdef TEST_ENERGY
		/*      fprintf( out_control->eov, "%6d%12.6f\n",
              workspace->reverse_map[j],
        // CEover1 * twbp->p_ovun1 * twbp->De_s, CEover3,
        CEover4 * (1.0 - workspace->dDelta_lp[j]) *
        (bo_ij->BO_pi + bo_ij->BO_pi2)
		 *///           /*CEover4 * (workspace->Delta[j]-workspace->Delta_lp[j])*/);
		//      fprintf( out_control->eov, "%6d%12.6f\n",
		//      fprintf( out_control->eov, "%6d%24.15e\n",
		//           system->my_atoms[j].orig_id,
		// CEover1 * twbp->p_ovun1 * twbp->De_s, CEover3,
		//           CEover4 * (1.0 - workspace->dDelta_lp[j]) *
		//           (bo_ij->BO_pi + bo_ij->BO_pi2)
		//           /*CEover4 * (workspace->Delta[j]-workspace->Delta_lp[j])*/);

		// CEunder4 * (1.0 - workspace->dDelta_lp[j]) *
		// (bo_ij->BO_pi + bo_ij->BO_pi2),
		// CEunder4 * (workspace->Delta[j] - workspace->Delta_lp[j]) );
#endif

#ifdef TEST_FORCES
		Add_dBO( system, lists, i, pj, CEover1 * twbp->p_ovun1 * twbp->De_s,
				workspace->f_ov ); // OvCoor - 1st term

		Add_dDelta( system, lists, j,
				CEover4 * (1.0 - dfvl*workspace->dDelta_lp[j]) *
				(bo_ij->BO_pi + bo_ij->BO_pi2),
				workspace->f_ov );   // OvCoor - 3a

		Add_dBOpinpi2( system, lists, i, pj,
				CEover4 * (workspace->Delta[j] -
						dfvl * workspace->Delta_lp_temp[j]),
						CEover4 * (workspace->Delta[j] -
								dfvl * workspace->Delta_lp_temp[j]),
								workspace->f_ov, workspace->f_ov ); // OvCoor - 3b

		Add_dDelta( system, lists, j,
				CEunder4 * (1.0 - dfvl*workspace->dDelta_lp[j]) *
				(bo_ij->BO_pi + bo_ij->BO_pi2),
				workspace->f_un ); // UnCoor - 2a

		Add_dBOpinpi2( system, lists, i, pj,
				CEunder4 * (workspace->Delta[j] -
						dfvl * workspace->Delta_lp_temp[j]),
						CEunder4 * (workspace->Delta[j] -
								dfvl * workspace->Delta_lp_temp[j]),
								workspace->f_un, workspace->f_un ); // UnCoor - 2b
#endif
	}

#ifdef TEST_ENERGY
	//fprintf( out_control->elp, "%6d%24.15e%24.15e%24.15e\n",
	//fprintf( out_control->elp, "%6d%12.4f%12.4f%12.4f\n",
	//     system->my_atoms[i].orig_id, workspace->nlp[i],
	//     e_lp, data->my_en.e_lp );

	//fprintf( out_control->eov, "%6d%24.15e%24.15e\n",
	fprintf( out_control->eov, "%6d%12.4f%12.4f\n",
			system->my_atoms[i].orig_id,
			e_ov, data->my_en.e_ov + data->my_en.e_un );

	//fprintf( out_control->eun, "%6d%24.15e%24.15e\n",
	fprintf( out_control->eun, "%6d%12.4f%12.4f\n",
			system->my_atoms[i].orig_id,
			e_un, data->my_en.e_ov + data->my_en.e_un );
#endif
	//}



}


CUDA_GLOBAL void Cuda_Atom_Energy_PostProcess( reax_list p_bonds, 
		storage p_workspace, int n )
{
	int i, pj;
	bond_data *sbond;
	//    bond_data *pbond;
	bond_data *sym_index_bond;
	reax_list *bonds;
	storage *workspace;

	i = blockIdx.x * blockDim.x + threadIdx.x;

	if ( i >= n )
	{
		return;
	}

	bonds = &p_bonds;
	workspace = &p_workspace;

	for ( pj = Cuda_Start_Index(i, bonds); pj < Cuda_End_Index(i, bonds); ++pj )
	{
		//        pbond = &bonds->bond_list[pj];
		//        dbond_index_bond = &bonds->bond_list[ pbond->dbond_index ];
		//        workspace->CdDelta[i] += dbond_index_bond->ae_CdDelta;

		sbond = &bonds->select.bond_list[pj];
		sym_index_bond = &bonds->select.bond_list[ sbond->sym_index ];
		workspace->CdDelta[i] += sym_index_bond->ae_CdDelta;
	}



}
