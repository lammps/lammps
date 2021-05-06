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

#include "cuda_bonds.h"

#include "cuda_list.h"

#include "index_utils.h"


CUDA_GLOBAL void Cuda_Bonds( reax_atom *my_atoms, global_parameters gp, 
        single_body_parameters *sbp, two_body_parameters *tbp, 
        storage p_workspace, reax_list p_bonds, int n, int num_atom_types, 
        real *e_bond )
{
    int i, j, pj;
    int start_i, end_i;
    int type_i, type_j;
    real ebond, pow_BOs_be2, exp_be12, CEbo;
    real gp3, gp4, gp7, gp10, gp37;
    real exphu, exphua1, exphub1, exphuov, hulpov, estriph;
    real decobdbo, decobdboua, decobdboub;
    single_body_parameters *sbp_i, *sbp_j;
    two_body_parameters *twbp;
    bond_order_data *bo_ij;
    reax_list *bonds;
    storage *workspace;

    i = blockIdx.x * blockDim.x + threadIdx.x;

    if ( i >= n )
    {
        return;
    }

    bonds = &p_bonds;
    workspace = &p_workspace;
    gp3 = gp.l[3];
    gp4 = gp.l[4];
    gp7 = gp.l[7];
    gp10 = gp.l[10];
    gp37 = (int) gp.l[37];

    start_i = Cuda_Start_Index( i, bonds );
    end_i = Cuda_End_Index( i, bonds );

    for ( pj = start_i; pj < end_i; ++pj )
    {
        j = bonds->select.bond_list[pj].nbr;

        if ( my_atoms[i].orig_id <= my_atoms[j].orig_id )
        {
            /* set the pointers */
            type_i = my_atoms[i].type;
            type_j = my_atoms[j].type;
            sbp_i = &sbp[type_i];
            sbp_j = &sbp[type_j];

            twbp = &tbp[ index_tbp(type_i,type_j, num_atom_types) ];
            bo_ij = &bonds->select.bond_list[pj].bo_data;

            /* calculate the constants */
            pow_BOs_be2 = POW( bo_ij->BO_s, twbp->p_be2 );
            exp_be12 = EXP( twbp->p_be1 * ( 1.0 - pow_BOs_be2 ) );
            CEbo = -twbp->De_s * exp_be12 * 
                ( 1.0 - twbp->p_be1 * twbp->p_be2 * pow_BOs_be2 );

            /* calculate the Bond Energy */
            ebond = -twbp->De_s * bo_ij->BO_s * exp_be12 
                -twbp->De_p * bo_ij->BO_pi 
                -twbp->De_pp * bo_ij->BO_pi2;
            e_bond[ i ] += ebond;

            /* calculate derivatives of Bond Orders */
            bo_ij->Cdbo += CEbo;

           // printf("cebo %f\n", bo_ij->Cdbo);



            bo_ij->Cdbopi -= (CEbo + twbp->De_p);
            bo_ij->Cdbopi2 -= (CEbo + twbp->De_pp);

#ifdef TEST_ENERGY
            //fprintf( out_control->ebond, "%6d%6d%24.15e%24.15e%24.15e\n",
            fprintf( out_control->ebond, "%6d%6d%12.4f%12.4f%12.4f\n",
                    system->my_atoms[i].orig_id, 
                    system->my_atoms[j].orig_id, 
                    bo_ij->BO, ebond, data->my_en.e_bond );
#endif

#ifdef TEST_FORCES
            Add_dBO( system, lists, i, pj, CEbo, workspace->f_be );
            Add_dBOpinpi2( system, lists, i, pj, 
                    -(CEbo + twbp->De_p), -(CEbo + twbp->De_pp), 
                    workspace->f_be, workspace->f_be );
#endif

            /* Stabilisation terminal triple bond */
            if ( bo_ij->BO >= 1.00 )
            {
                if ( gp37 == 2 ||
                        (sbp_i->mass == 12.0000 && sbp_j->mass == 15.9990) || 
                        (sbp_j->mass == 12.0000 && sbp_i->mass == 15.9990) )
                {
                    exphu = EXP( -gp7 * SQR(bo_ij->BO - 2.50) );
                    exphua1 = EXP(-gp3 * (workspace->total_bond_order[i]-bo_ij->BO));
                    exphub1 = EXP(-gp3 * (workspace->total_bond_order[j]-bo_ij->BO));
                    exphuov = EXP(gp4 * (workspace->Delta[i] + workspace->Delta[j]));
                    hulpov = 1.0 / (1.0 + 25.0 * exphuov);

                    estriph = gp10 * exphu * hulpov * (exphua1 + exphub1);
                    e_bond[i] += estriph;

                    decobdbo = gp10 * exphu * hulpov * (exphua1 + exphub1) *
                        ( gp3 - 2.0 * gp7 * (bo_ij->BO-2.50) );
                    decobdboua = -gp10 * exphu * hulpov *
                        (gp3*exphua1 + 25.0*gp4*exphuov*hulpov*(exphua1+exphub1));
                    decobdboub = -gp10 * exphu * hulpov *
                        (gp3*exphub1 + 25.0*gp4*exphuov*hulpov*(exphua1+exphub1));

                   // printf(" %f,%f,%f,%f,%f,%f,%f\n",gp10,exphu,hulpov,exphua1,exphub1,gp3,gp7);
                    //printf("%d,%d,%f\n",i,j,workspace->total_bond_order[j]);

                    bo_ij->Cdbo += decobdbo;

                   // printf("%f,%f\n",workspace->total_bond_order[j], bo_ij->Cdbo);

                    workspace->CdDelta[i] += decobdboua;
                    workspace->CdDelta[j] += decobdboub;

#ifdef TEST_ENERGY
                    //fprintf( out_control->ebond, 
                    //  "%6d%6d%24.15e%24.15e%24.15e%24.15e\n",
                    //  system->my_atoms[i].orig_id, system->my_atoms[j].orig_id,
                    //  estriph, decobdbo, decobdboua, decobdboub );
#endif

#ifdef TEST_FORCES
                    Add_dBO( system, lists, i, pj, decobdbo, workspace->f_be );
                    Add_dDelta( system, lists, i, decobdboua, workspace->f_be );
                    Add_dDelta( system, lists, j, decobdboub, workspace->f_be );
#endif
                }
            }
        }
    }
}
