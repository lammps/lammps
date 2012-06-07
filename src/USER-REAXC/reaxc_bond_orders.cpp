/*----------------------------------------------------------------------
  PuReMD - Purdue ReaxFF Molecular Dynamics Program

  Copyright (2010) Purdue University
  Hasan Metin Aktulga, hmaktulga@lbl.gov
  Joseph Fogarty, jcfogart@mail.usf.edu
  Sagar Pandit, pandit@usf.edu
  Ananth Y Grama, ayg@cs.purdue.edu

  Please cite the related publication:
  H. M. Aktulga, J. C. Fogarty, S. A. Pandit, A. Y. Grama,
  "Parallel Reactive Molecular Dynamics: Numerical Methods and
  Algorithmic Techniques", Parallel Computing, in press.

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

#include "pair_reax_c.h"
#include "reaxc_types.h"
#if defined(PURE_REAX)
#include "bond_orders.h"
#include "list.h"
#include "vector.h"
#include "io_tools.h"
#elif defined(LAMMPS_REAX)
#include "reaxc_bond_orders.h"
#include "reaxc_list.h"
#include "reaxc_vector.h"
#endif

#ifdef TEST_FORCES
void Get_dBO( reax_system *system, reax_list **lists,
              int i, int pj, real C, rvec *v )
{
  reax_list *bonds = (*lists) + BONDS;
  reax_list *dBOs = (*lists) + DBOS;
  int start_pj, end_pj, k;

  pj = bonds->select.bond_list[pj].dbond_index;
  start_pj = Start_Index(pj, dBOs);
  end_pj = End_Index(pj, dBOs);

  for( k = start_pj; k < end_pj; ++k )
    rvec_Scale( v[dBOs->select.dbo_list[k].wrt],
                C, dBOs->select.dbo_list[k].dBO );
}


void Get_dBOpinpi2( reax_system *system, reax_list **lists,
                    int i, int pj, real Cpi, real Cpi2, rvec *vpi, rvec *vpi2 )
{
  reax_list *bonds = (*lists) + BONDS;
  reax_list *dBOs = (*lists) + DBOS;
  dbond_data *dbo_k;
  int start_pj, end_pj, k;

  pj = bonds->select.bond_list[pj].dbond_index;
  start_pj = Start_Index(pj, dBOs);
  end_pj = End_Index(pj, dBOs);

  for( k = start_pj; k < end_pj; ++k )
    {
      dbo_k = &(dBOs->select.dbo_list[k]);
      rvec_Scale( vpi[dbo_k->wrt], Cpi, dbo_k->dBOpi );
      rvec_Scale( vpi2[dbo_k->wrt], Cpi2, dbo_k->dBOpi2 );
    }
}


void Add_dBO( reax_system *system, reax_list **lists,
              int i, int pj, real C, rvec *v )
{
  reax_list *bonds = (*lists) + BONDS;
  reax_list *dBOs = (*lists) + DBOS;
  int start_pj, end_pj, k;

  pj = bonds->select.bond_list[pj].dbond_index;
  start_pj = Start_Index(pj, dBOs);
  end_pj = End_Index(pj, dBOs);
  //fprintf( stderr, "i=%d j=%d start=%d end=%d\n", i, pj, start_pj, end_pj );

  for( k = start_pj; k < end_pj; ++k )
    rvec_ScaledAdd( v[dBOs->select.dbo_list[k].wrt],
                    C, dBOs->select.dbo_list[k].dBO );

}


void Add_dBOpinpi2( reax_system *system, reax_list **lists,
                    int i, int pj, real Cpi, real Cpi2, rvec *vpi, rvec *vpi2 )
{
  reax_list *bonds = (*lists) + BONDS;
  reax_list *dBOs = (*lists) + DBOS;
  dbond_data *dbo_k;
  int start_pj, end_pj, k;

  pj = bonds->select.bond_list[pj].dbond_index;
  start_pj = Start_Index(pj, dBOs);
  end_pj = End_Index(pj, dBOs);

  for( k = start_pj; k < end_pj; ++k )
    {
      dbo_k = &(dBOs->select.dbo_list[k]);
      rvec_ScaledAdd( vpi[dbo_k->wrt], Cpi, dbo_k->dBOpi );
      rvec_ScaledAdd( vpi2[dbo_k->wrt], Cpi2, dbo_k->dBOpi2 );
    }
}


void Add_dBO_to_Forces( reax_system *system, reax_list **lists,
                        int i, int pj, real C )
{
  reax_list *bonds = (*lists) + BONDS;
  reax_list *dBOs = (*lists) + DBOS;
  int start_pj, end_pj, k;

  pj = bonds->select.bond_list[pj].dbond_index;
  start_pj = Start_Index(pj, dBOs);
  end_pj = End_Index(pj, dBOs);

  for( k = start_pj; k < end_pj; ++k )
    rvec_ScaledAdd( system->my_atoms[dBOs->select.dbo_list[k].wrt].f,
                    C, dBOs->select.dbo_list[k].dBO );
}


void Add_dBOpinpi2_to_Forces( reax_system *system, reax_list **lists,
                              int i, int pj, real Cpi, real Cpi2 )
{
  reax_list *bonds = (*lists) + BONDS;
  reax_list *dBOs = (*lists) + DBOS;
  dbond_data *dbo_k;
  int start_pj, end_pj, k;

  pj = bonds->select.bond_list[pj].dbond_index;
  start_pj = Start_Index(pj, dBOs);
  end_pj = End_Index(pj, dBOs);

  for( k = start_pj; k < end_pj; ++k )
    {
      dbo_k = &(dBOs->select.dbo_list[k]);
      rvec_ScaledAdd( system->my_atoms[dbo_k->wrt].f, Cpi, dbo_k->dBOpi );
      rvec_ScaledAdd( system->my_atoms[dbo_k->wrt].f, Cpi2, dbo_k->dBOpi2 );
    }
}


void Add_dDelta( reax_system *system, reax_list **lists, int i, real C, rvec *v )
{
  reax_list *dDeltas = &((*lists)[DDELTAS]);
  int start = Start_Index(i, dDeltas);
  int end = End_Index(i, dDeltas);
  int k;

  for( k = start; k < end; ++k )
    rvec_ScaledAdd( v[dDeltas->select.dDelta_list[k].wrt],
                    C, dDeltas->select.dDelta_list[k].dVal );
}


void Add_dDelta_to_Forces( reax_system *system, reax_list **lists,
                           int i, real C )
{
  reax_list *dDeltas = &((*lists)[DDELTAS]);
  int start = Start_Index(i, dDeltas);
  int end = End_Index(i, dDeltas);
  int k;

  for( k = start; k < end; ++k )
    rvec_ScaledAdd( system->my_atoms[dDeltas->select.dDelta_list[k].wrt].f,
                    C, dDeltas->select.dDelta_list[k].dVal );
}



void Calculate_dBO( int i, int pj,
                    storage *workspace, reax_list **lists, int *top )
{
  /* Initializations */
  reax_list *bonds, *dBOs;
  int j, k, l, start_i, end_i, end_j;
  bond_data *nbr_l, *nbr_k;
  bond_order_data *bo_ij;
  dbond_data *top_dbo;

  //  rvec due_j[1000], due_i[1000];
  //  rvec due_j_pi[1000], due_i_pi[1000];

  //  memset(due_j, 0, sizeof(rvec)*1000 );
  //  memset(due_i, 0, sizeof(rvec)*1000 );
  //  memset(due_j_pi, 0, sizeof(rvec)*1000 );
  //  memset(due_i_pi, 0, sizeof(rvec)*1000 );

  bonds = *lists + BONDS;
  dBOs  = *lists + DBOS;

  start_i = Start_Index(i, bonds);
  end_i = End_Index(i, bonds);

  j = bonds->select.bond_list[pj].nbr;
  l = Start_Index(j, bonds);
  end_j = End_Index(j, bonds);

  bo_ij = &(bonds->select.bond_list[pj].bo_data);
  top_dbo = &(dBOs->select.dbo_list[ (*top) ]);

  for( k = start_i; k < end_i; ++k ) {
    nbr_k = &(bonds->select.bond_list[k]);

    for( ; l < end_j && bonds->select.bond_list[l].nbr < nbr_k->nbr; ++l ) {
      /* These are the neighbors of j which are not in the nbr_list of i
         Note that they might also include i! */
      nbr_l = &(bonds->select.bond_list[l]);
      top_dbo->wrt = nbr_l->nbr;

      //rvec_ScaledAdd( due_j[top_dbo->wrt],
      //                -bo_ij->BO * bo_ij->A2_ji, nbr_l->bo_data.dBOp );

      /*3rd, dBO*/
      rvec_Scale( top_dbo->dBO, -bo_ij->C3dbo, nbr_l->bo_data.dBOp );
      /*4th, dBOpi*/
      rvec_Scale( top_dbo->dBOpi, -bo_ij->C4dbopi, nbr_l->bo_data.dBOp );
      /*4th, dBOpp*/
      rvec_Scale( top_dbo->dBOpi2, -bo_ij->C4dbopi2, nbr_l->bo_data.dBOp );

      if( nbr_l->nbr == i ) {
        /* do the adjustments on i */
        //rvec_ScaledAdd( due_i[i], bo_ij->A0_ij+bo_ij->BO*bo_ij->A1_ij,
        //                bo_ij->dBOp );  /*1st, dBO*/
        //rvec_ScaledAdd( due_i[i], bo_ij->BO * bo_ij->A2_ij,
        //                workspace->dDeltap_self[i] ); /*2nd, dBO*/

        /*1st, dBO*/
        rvec_ScaledAdd( top_dbo->dBO, bo_ij->C1dbo, bo_ij->dBOp );
        /*2nd, dBO*/
        rvec_ScaledAdd( top_dbo->dBO,
                        bo_ij->C2dbo, workspace->dDeltap_self[i] );

        /*1st, dBOpi*/
        rvec_ScaledAdd(top_dbo->dBOpi,bo_ij->C1dbopi,bo_ij->dln_BOp_pi );
        /*2nd, dBOpi*/
        rvec_ScaledAdd( top_dbo->dBOpi, bo_ij->C2dbopi, bo_ij->dBOp );
        /*3rd, dBOpi*/
        rvec_ScaledAdd( top_dbo->dBOpi,
                        bo_ij->C3dbopi, workspace->dDeltap_self[i] );

        /*1st, dBO_p*/
        rvec_ScaledAdd( top_dbo->dBOpi2,
                        bo_ij->C1dbopi2, bo_ij->dln_BOp_pi2 );
        /*2nd, dBO_p*/
        rvec_ScaledAdd( top_dbo->dBOpi2, bo_ij->C2dbopi2, bo_ij->dBOp );
        /*3rd, dBO_p*/
        rvec_ScaledAdd( top_dbo->dBOpi2,
                        bo_ij->C3dbopi2, workspace->dDeltap_self[i] );
      }

      //rvec_Add( workspace->dDelta[nbr_l->nbr], top_dbo->dBO );
      ++(*top), ++top_dbo;
    }

    /* Now we are processing neighbor k of i. */
    top_dbo->wrt = nbr_k->nbr;

    //2nd, dBO
    //rvec_ScaledAdd( due_i[top_dbo->wrt],
    //                -bo_ij->BO * bo_ij->A2_ij, nbr_k->bo_data.dBOp );

    /*2nd, dBO*/
    rvec_Scale( top_dbo->dBO, -bo_ij->C2dbo, nbr_k->bo_data.dBOp );
    /*3rd, dBOpi*/
    rvec_Scale( top_dbo->dBOpi, -bo_ij->C3dbopi, nbr_k->bo_data.dBOp );
    /*3rd, dBOpp*/
    rvec_Scale( top_dbo->dBOpi2, -bo_ij->C3dbopi2, nbr_k->bo_data.dBOp );

    if( l < end_j && bonds->select.bond_list[l].nbr == nbr_k->nbr ) {
      /* This is a common neighbor of i and j. */
      nbr_l = &(bonds->select.bond_list[l]);

      /*3rd, dBO*/
      //rvec_ScaledAdd( due_j[top_dbo->wrt],
      //                -bo_ij->BO * bo_ij->A2_ji, nbr_l->bo_data.dBOp );

      /*3rd, dBO*/
      rvec_ScaledAdd( top_dbo->dBO, -bo_ij->C3dbo, nbr_l->bo_data.dBOp );
      /*4th, dBOpi*/
      rvec_ScaledAdd(top_dbo->dBOpi, -bo_ij->C4dbopi, nbr_l->bo_data.dBOp);
      /*4th, dBOpp*/
      rvec_ScaledAdd(top_dbo->dBOpi2,-bo_ij->C4dbopi2,nbr_l->bo_data.dBOp);
      ++l;
    }
    else if( k == pj ) {
      /* This negihbor is j. */
      //rvec_ScaledAdd( due_j[j], -(bo_ij->A0_ij+bo_ij->BO*bo_ij->A1_ij),
      //                bo_ij->dBOp );  /*1st, dBO*/
      //rvec_ScaledAdd( due_j[j], bo_ij->BO * bo_ij->A2_ji,
      //                workspace->dDeltap_self[j] );  /*3rd, dBO*/

      /*1st, dBO*/
      rvec_ScaledAdd( top_dbo->dBO, -bo_ij->C1dbo, bo_ij->dBOp );
      /*3rd, dBO*/
      rvec_ScaledAdd(top_dbo->dBO,bo_ij->C3dbo,workspace->dDeltap_self[j]);

      /*1st, dBOpi*/
      rvec_ScaledAdd( top_dbo->dBOpi, -bo_ij->C1dbopi, bo_ij->dln_BOp_pi );
      /*2nd, dBOpi*/
      rvec_ScaledAdd( top_dbo->dBOpi, -bo_ij->C2dbopi, bo_ij->dBOp );
      /*4th, dBOpi*/
      rvec_ScaledAdd( top_dbo->dBOpi,
                      bo_ij->C4dbopi, workspace->dDeltap_self[j] );

      /*1st, dBOpi2*/
      rvec_ScaledAdd(top_dbo->dBOpi2,-bo_ij->C1dbopi2,bo_ij->dln_BOp_pi2);
      /*2nd, dBOpi2*/
      rvec_ScaledAdd(top_dbo->dBOpi2,-bo_ij->C2dbopi2,bo_ij->dBOp );
      /*4th, dBOpi2*/
      rvec_ScaledAdd(top_dbo->dBOpi2, bo_ij->C4dbopi2,
                     workspace->dDeltap_self[j] );
    }

    // rvec_Add( workspace->dDelta[nbr_k->nbr], top_dbo->dBO );
    ++(*top), ++top_dbo;
  }

  for( ; l < end_j; ++l ) {
    /* These are the remaining neighbors of j which are not in the
       neighbor_list of i. Note that they might also include i!*/
    nbr_l = &(bonds->select.bond_list[l]);
    top_dbo->wrt = nbr_l->nbr;

    // fprintf( stdout, "\tl: %d nbr:%d\n", l, nbr_l->nbr+1 );

    // rvec_ScaledAdd( due_j[top_dbo->wrt],
    //                 -bo_ij->BO * bo_ij->A2_ji, nbr_l->bo_data.dBOp );

    /*3rd, dBO*/
    rvec_Scale( top_dbo->dBO, -bo_ij->C3dbo, nbr_l->bo_data.dBOp );
    /*4th, dBOpi*/
    rvec_Scale( top_dbo->dBOpi, -bo_ij->C4dbopi, nbr_l->bo_data.dBOp );
    /*4th, dBOpp*/
    rvec_Scale( top_dbo->dBOpi2, -bo_ij->C4dbopi2, nbr_l->bo_data.dBOp );

    if( nbr_l->nbr == i )
      {
        /* do the adjustments on i */
        //rvec_ScaledAdd( due_i[i], bo_ij->A0_ij + bo_ij->BO * bo_ij->A1_ij,
        //                bo_ij->dBOp );  /*1st, dBO*/
        //rvec_ScaledAdd( due_i[i], bo_ij->BO * bo_ij->A2_ij,
        //                workspace->dDeltap_self[i] ); /*2nd, dBO*/

        /*1st, dBO*/
        rvec_ScaledAdd( top_dbo->dBO, bo_ij->C1dbo, bo_ij->dBOp );
        /*2nd, dBO*/
        rvec_ScaledAdd(top_dbo->dBO,bo_ij->C2dbo,workspace->dDeltap_self[i]);

        /*1st, dBO_p*/
        rvec_ScaledAdd( top_dbo->dBOpi, bo_ij->C1dbopi, bo_ij->dln_BOp_pi );
        /*2nd, dBOpi*/
        rvec_ScaledAdd( top_dbo->dBOpi, bo_ij->C2dbopi, bo_ij->dBOp );
        /*3rd,dBOpi*/
        rvec_ScaledAdd( top_dbo->dBOpi,
                        bo_ij->C3dbopi, workspace->dDeltap_self[i] );

        /*1st, dBO_p*/
        rvec_ScaledAdd(top_dbo->dBOpi2, bo_ij->C1dbopi2, bo_ij->dln_BOp_pi2);
        /*2nd, dBO_p*/
        rvec_ScaledAdd(top_dbo->dBOpi2, bo_ij->C2dbopi2, bo_ij->dBOp);
        /*3rd,dBO_p*/
        rvec_ScaledAdd(top_dbo->dBOpi2,
                       bo_ij->C3dbopi2, workspace->dDeltap_self[i]);
      }

    // rvec_Add( workspace->dDelta[nbr_l->nbr], top_dbo->dBO );
    ++(*top), ++top_dbo;
  }

  //  for( k = 0; k < 21; ++k ){
  //  fprintf( stderr, "%d %d %d, due_i:[%g %g %g]\n",
  //           i+1, j+1, k+1, due_i[k][0], due_i[k][1], due_i[k][2] );
  //  fprintf( stderr, "%d %d %d, due_j:[%g %g %g]\n",
  //           i+1, j+1, k+1, due_j[k][0], due_j[k][1], due_j[k][2] );
}
#endif


void Add_dBond_to_Forces_NPT( int i, int pj, simulation_data *data,
                              storage *workspace, reax_list **lists )
{
  reax_list *bonds = (*lists) + BONDS;
  bond_data *nbr_j, *nbr_k;
  bond_order_data *bo_ij, *bo_ji;
  dbond_coefficients coef;
  rvec temp, ext_press;
  ivec rel_box;
  int pk, k, j;

  /* Initializations */
  nbr_j = &(bonds->select.bond_list[pj]);
  j = nbr_j->nbr;
  bo_ij = &(nbr_j->bo_data);
  bo_ji = &(bonds->select.bond_list[ nbr_j->sym_index ].bo_data);

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
  for( pk = Start_Index(i, bonds); pk < End_Index(i, bonds); ++pk ) {
    nbr_k = &(bonds->select.bond_list[pk]);
    k = nbr_k->nbr;

    rvec_Scale(temp, -coef.C2dbo, nbr_k->bo_data.dBOp);       /*2nd, dBO*/
    rvec_ScaledAdd(temp, -coef.C2dDelta, nbr_k->bo_data.dBOp);/*dDelta*/
    rvec_ScaledAdd(temp, -coef.C3dbopi, nbr_k->bo_data.dBOp); /*3rd, dBOpi*/
    rvec_ScaledAdd(temp, -coef.C3dbopi2, nbr_k->bo_data.dBOp);/*3rd, dBOpi2*/

    /* force */
    rvec_Add( workspace->f[k], temp );
    /* pressure */
    rvec_iMultiply( ext_press, nbr_k->rel_box, temp );
    rvec_Add( data->my_ext_press, ext_press );

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


  /******************************************************
   * forces and pressure related to atom j               *
   * first neighbors of atom j                           *
   ******************************************************/
  for( pk = Start_Index(j, bonds); pk < End_Index(j, bonds); ++pk ) {
    nbr_k = &(bonds->select.bond_list[pk]);
    k = nbr_k->nbr;

    rvec_Scale( temp, -coef.C3dbo, nbr_k->bo_data.dBOp );      /*3rd,dBO*/
    rvec_ScaledAdd( temp, -coef.C3dDelta, nbr_k->bo_data.dBOp);/*dDelta*/
    rvec_ScaledAdd( temp, -coef.C4dbopi, nbr_k->bo_data.dBOp); /*4th,dBOpi*/
    rvec_ScaledAdd( temp, -coef.C4dbopi2, nbr_k->bo_data.dBOp);/*4th,dBOpi2*/

    /* force */
    rvec_Add( workspace->f[k], temp );
    /* pressure */
    if( k != i ) {
      ivec_Sum( rel_box, nbr_k->rel_box, nbr_j->rel_box ); //rel_box(k, i)
      rvec_iMultiply( ext_press, rel_box, temp );
      rvec_Add( data->my_ext_press, ext_press );

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



void Add_dBond_to_Forces( reax_system *system, int i, int pj,
                          storage *workspace, reax_list **lists )
{
  reax_list *bonds = (*lists) + BONDS;
  bond_data *nbr_j, *nbr_k;
  bond_order_data *bo_ij, *bo_ji;
  dbond_coefficients coef;
  int pk, k, j;

  /* Virial Tallying variables */
  real f_scaler;
  rvec fi_tmp, fj_tmp, fk_tmp, delij, delji, delki, delkj, temp;

  /* Initializations */
  nbr_j = &(bonds->select.bond_list[pj]);
  j = nbr_j->nbr;
  bo_ij = &(nbr_j->bo_data);
  bo_ji = &(bonds->select.bond_list[ nbr_j->sym_index ].bo_data);

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

  // forces on i
  rvec_Scale(           temp, coef.C1dbo,    bo_ij->dBOp );
  rvec_ScaledAdd( temp, coef.C2dbo,    workspace->dDeltap_self[i] );
  rvec_ScaledAdd( temp, coef.C1dDelta, bo_ij->dBOp );
  rvec_ScaledAdd( temp, coef.C2dDelta, workspace->dDeltap_self[i] );
  rvec_ScaledAdd( temp, coef.C1dbopi,  bo_ij->dln_BOp_pi );
  rvec_ScaledAdd( temp, coef.C2dbopi,  bo_ij->dBOp );
  rvec_ScaledAdd( temp, coef.C3dbopi,  workspace->dDeltap_self[i]);
  rvec_ScaledAdd( temp, coef.C1dbopi2, bo_ij->dln_BOp_pi2 );
  rvec_ScaledAdd( temp, coef.C2dbopi2, bo_ij->dBOp );
  rvec_ScaledAdd( temp, coef.C3dbopi2, workspace->dDeltap_self[i] );
  rvec_Add( workspace->f[i], temp );

  if( system->pair_ptr->vflag_atom) {
    rvec_Scale(fi_tmp, -1.0, temp);
    rvec_ScaledSum( delij, 1., system->my_atoms[i].x,-1., system->my_atoms[j].x );
    system->pair_ptr->v_tally(i,fi_tmp,delij);
  }

  // forces on j
  rvec_Scale(           temp, -coef.C1dbo,    bo_ij->dBOp );
  rvec_ScaledAdd( temp,  coef.C3dbo,    workspace->dDeltap_self[j] );
  rvec_ScaledAdd( temp, -coef.C1dDelta, bo_ij->dBOp );
  rvec_ScaledAdd( temp,  coef.C3dDelta, workspace->dDeltap_self[j]);
  rvec_ScaledAdd( temp, -coef.C1dbopi,  bo_ij->dln_BOp_pi );
  rvec_ScaledAdd( temp, -coef.C2dbopi,  bo_ij->dBOp );
  rvec_ScaledAdd( temp,  coef.C4dbopi,  workspace->dDeltap_self[j]);
  rvec_ScaledAdd( temp, -coef.C1dbopi2, bo_ij->dln_BOp_pi2 );
  rvec_ScaledAdd( temp, -coef.C2dbopi2, bo_ij->dBOp );
  rvec_ScaledAdd( temp,  coef.C4dbopi2, workspace->dDeltap_self[j]);
  rvec_Add( workspace->f[j], temp );

  if( system->pair_ptr->vflag_atom) {
    rvec_Scale(fj_tmp, -1.0, temp);
    rvec_ScaledSum( delji, 1., system->my_atoms[j].x,-1., system->my_atoms[i].x );
    system->pair_ptr->v_tally(j,fj_tmp,delji);
  }

  // forces on k: i neighbor
  for( pk = Start_Index(i, bonds); pk < End_Index(i, bonds); ++pk ) {
    nbr_k = &(bonds->select.bond_list[pk]);
    k = nbr_k->nbr;

    rvec_Scale(     temp, -coef.C2dbo,    nbr_k->bo_data.dBOp);
    rvec_ScaledAdd( temp, -coef.C2dDelta, nbr_k->bo_data.dBOp);
    rvec_ScaledAdd( temp, -coef.C3dbopi,  nbr_k->bo_data.dBOp);
    rvec_ScaledAdd( temp, -coef.C3dbopi2, nbr_k->bo_data.dBOp);
    rvec_Add( workspace->f[k], temp );

    if( system->pair_ptr->vflag_atom ) {
      rvec_Scale(fk_tmp, -1.0, temp);
      rvec_ScaledSum(delki,1.,system->my_atoms[k].x,-1.,system->my_atoms[i].x);
      system->pair_ptr->v_tally(k,fk_tmp,delki);
      rvec_ScaledSum(delkj,1.,system->my_atoms[k].x,-1.,system->my_atoms[j].x);
      system->pair_ptr->v_tally(k,fk_tmp,delkj);
    }
  }

  // forces on k: j neighbor
  for( pk = Start_Index(j, bonds); pk < End_Index(j, bonds); ++pk ) {
    nbr_k = &(bonds->select.bond_list[pk]);
    k = nbr_k->nbr;

    rvec_Scale(     temp, -coef.C3dbo,    nbr_k->bo_data.dBOp );
    rvec_ScaledAdd( temp, -coef.C3dDelta, nbr_k->bo_data.dBOp);
    rvec_ScaledAdd( temp, -coef.C4dbopi,  nbr_k->bo_data.dBOp);
    rvec_ScaledAdd( temp, -coef.C4dbopi2, nbr_k->bo_data.dBOp);
    rvec_Add( workspace->f[k], temp );

    if( system->pair_ptr->vflag_atom ) {
      rvec_Scale(fk_tmp, -1.0, temp);
      rvec_ScaledSum(delki,1.,system->my_atoms[k].x,-1.,system->my_atoms[i].x);
      system->pair_ptr->v_tally(k,fk_tmp,delki);
      rvec_ScaledSum(delkj,1.,system->my_atoms[k].x,-1.,system->my_atoms[j].x);
      system->pair_ptr->v_tally(k,fk_tmp,delkj);
    }
  }

}


int BOp( storage *workspace, reax_list *bonds, real bo_cut,
         int i, int btop_i, far_neighbor_data *nbr_pj,
         single_body_parameters *sbp_i, single_body_parameters *sbp_j,
         two_body_parameters *twbp ) {
  int j, btop_j;
  real r2, C12, C34, C56;
  real Cln_BOp_s, Cln_BOp_pi, Cln_BOp_pi2;
  real BO, BO_s, BO_pi, BO_pi2;
  bond_data *ibond, *jbond;
  bond_order_data *bo_ij, *bo_ji;

  j = nbr_pj->nbr;
  r2 = SQR(nbr_pj->d);

  if( sbp_i->r_s > 0.0 && sbp_j->r_s > 0.0 ) {
    C12 = twbp->p_bo1 * pow( nbr_pj->d / twbp->r_s, twbp->p_bo2 );
    BO_s = (1.0 + bo_cut) * exp( C12 );
  }
  else BO_s = C12 = 0.0;

  if( sbp_i->r_pi > 0.0 && sbp_j->r_pi > 0.0 ) {
    C34 = twbp->p_bo3 * pow( nbr_pj->d / twbp->r_p, twbp->p_bo4 );
    BO_pi = exp( C34 );
  }
  else BO_pi = C34 = 0.0;

  if( sbp_i->r_pi_pi > 0.0 && sbp_j->r_pi_pi > 0.0 ) {
    C56 = twbp->p_bo5 * pow( nbr_pj->d / twbp->r_pp, twbp->p_bo6 );
    BO_pi2= exp( C56 );
  }
  else BO_pi2 = C56 = 0.0;

  /* Initially BO values are the uncorrected ones, page 1 */
  BO = BO_s + BO_pi + BO_pi2;

  if( BO >= bo_cut ) {
    /****** bonds i-j and j-i ******/
    ibond = &( bonds->select.bond_list[btop_i] );
    btop_j = End_Index( j, bonds );
    jbond = &(bonds->select.bond_list[btop_j]);

    ibond->nbr = j;
    jbond->nbr = i;
    ibond->d = nbr_pj->d;
    jbond->d = nbr_pj->d;
    rvec_Copy( ibond->dvec, nbr_pj->dvec );
    rvec_Scale( jbond->dvec, -1, nbr_pj->dvec );
    ivec_Copy( ibond->rel_box, nbr_pj->rel_box );
    ivec_Scale( jbond->rel_box, -1, nbr_pj->rel_box );
    ibond->dbond_index = btop_i;
    jbond->dbond_index = btop_i;
    ibond->sym_index = btop_j;
    jbond->sym_index = btop_i;
    Set_End_Index( j, btop_j+1, bonds );

    bo_ij = &( ibond->bo_data );
    bo_ji = &( jbond->bo_data );
    bo_ji->BO = bo_ij->BO = BO;
    bo_ji->BO_s = bo_ij->BO_s = BO_s;
    bo_ji->BO_pi = bo_ij->BO_pi = BO_pi;
    bo_ji->BO_pi2 = bo_ij->BO_pi2 = BO_pi2;

    /* Bond Order page2-3, derivative of total bond order prime */
    Cln_BOp_s = twbp->p_bo2 * C12 / r2;
    Cln_BOp_pi = twbp->p_bo4 * C34 / r2;
    Cln_BOp_pi2 = twbp->p_bo6 * C56 / r2;

    /* Only dln_BOp_xx wrt. dr_i is stored here, note that
       dln_BOp_xx/dr_i = -dln_BOp_xx/dr_j and all others are 0 */
    rvec_Scale(bo_ij->dln_BOp_s,-bo_ij->BO_s*Cln_BOp_s,ibond->dvec);
    rvec_Scale(bo_ij->dln_BOp_pi,-bo_ij->BO_pi*Cln_BOp_pi,ibond->dvec);
    rvec_Scale(bo_ij->dln_BOp_pi2,
               -bo_ij->BO_pi2*Cln_BOp_pi2,ibond->dvec);
    rvec_Scale(bo_ji->dln_BOp_s, -1., bo_ij->dln_BOp_s);
    rvec_Scale(bo_ji->dln_BOp_pi, -1., bo_ij->dln_BOp_pi );
    rvec_Scale(bo_ji->dln_BOp_pi2, -1., bo_ij->dln_BOp_pi2 );

    /* Only dBOp wrt. dr_i is stored here, note that
       dBOp/dr_i = -dBOp/dr_j and all others are 0 */
    rvec_Scale( bo_ij->dBOp,
                -(bo_ij->BO_s * Cln_BOp_s +
                  bo_ij->BO_pi * Cln_BOp_pi +
                  bo_ij->BO_pi2 * Cln_BOp_pi2), ibond->dvec );
    rvec_Scale( bo_ji->dBOp, -1., bo_ij->dBOp );

    rvec_Add( workspace->dDeltap_self[i], bo_ij->dBOp );
    rvec_Add( workspace->dDeltap_self[j], bo_ji->dBOp );

    bo_ij->BO_s -= bo_cut;
    bo_ij->BO -= bo_cut;
    bo_ji->BO_s -= bo_cut;
    bo_ji->BO -= bo_cut;
    workspace->total_bond_order[i] += bo_ij->BO; //currently total_BOp
    workspace->total_bond_order[j] += bo_ji->BO; //currently total_BOp
    bo_ij->Cdbo = bo_ij->Cdbopi = bo_ij->Cdbopi2 = 0.0;
    bo_ji->Cdbo = bo_ji->Cdbopi = bo_ji->Cdbopi2 = 0.0;

    /*fprintf( stderr, "%d %d %g %g %g\n",
      i+1, j+1, bo_ij->BO, bo_ij->BO_pi, bo_ij->BO_pi2 );*/

    /*fprintf( stderr, "Cln_BOp_s: %f, pbo2: %f, C12:%f\n",
      Cln_BOp_s, twbp->p_bo2, C12 );
      fprintf( stderr, "Cln_BOp_pi: %f, pbo4: %f, C34:%f\n",
      Cln_BOp_pi, twbp->p_bo4, C34 );
      fprintf( stderr, "Cln_BOp_pi2: %f, pbo6: %f, C56:%f\n",
      Cln_BOp_pi2, twbp->p_bo6, C56 );*/
    /*fprintf(stderr, "pbo1: %f, pbo2:%f\n", twbp->p_bo1, twbp->p_bo2);
      fprintf(stderr, "pbo3: %f, pbo4:%f\n", twbp->p_bo3, twbp->p_bo4);
      fprintf(stderr, "pbo5: %f, pbo6:%f\n", twbp->p_bo5, twbp->p_bo6);
      fprintf( stderr, "r_s: %f, r_p: %f, r_pp: %f\n",
      twbp->r_s, twbp->r_p, twbp->r_pp );
      fprintf( stderr, "C12: %g, C34:%g, C56:%g\n", C12, C34, C56 );*/

    /*fprintf( stderr, "\tfactors: %g %g %g\n",
      -(bo_ij->BO_s * Cln_BOp_s + bo_ij->BO_pi * Cln_BOp_pi +
      bo_ij->BO_pi2 * Cln_BOp_pp),
      -bo_ij->BO_pi * Cln_BOp_pi, -bo_ij->BO_pi2 * Cln_BOp_pi2 );*/
    /*fprintf( stderr, "dBOpi:\t[%g, %g, %g]\n",
      bo_ij->dBOp[0], bo_ij->dBOp[1], bo_ij->dBOp[2] );
      fprintf( stderr, "dBOpi:\t[%g, %g, %g]\n",
      bo_ij->dln_BOp_pi[0], bo_ij->dln_BOp_pi[1],
      bo_ij->dln_BOp_pi[2] );
      fprintf( stderr, "dBOpi2:\t[%g, %g, %g]\n\n",
      bo_ij->dln_BOp_pi2[0], bo_ij->dln_BOp_pi2[1],
      bo_ij->dln_BOp_pi2[2] );*/

    return 1;
  }

  return 0;
}


int compare_bonds( const void *p1, const void *p2 )
{
  return ((bond_data *)p1)->nbr - ((bond_data *)p2)->nbr;
}


void BO( reax_system *system, control_params *control, simulation_data *data,
         storage *workspace, reax_list **lists, output_controls *out_control )
{
  int i, j, pj, type_i, type_j;
  int start_i, end_i, sym_index, num_bonds;
  real val_i, Deltap_i, Deltap_boc_i;
  real val_j, Deltap_j, Deltap_boc_j;
  real f1, f2, f3, f4, f5, f4f5, exp_f4, exp_f5;
  real exp_p1i,        exp_p2i, exp_p1j, exp_p2j;
  real temp, u1_ij, u1_ji, Cf1A_ij, Cf1B_ij, Cf1_ij, Cf1_ji;
  real Cf45_ij, Cf45_ji, p_lp1; //u_ij, u_ji
  real A0_ij, A1_ij, A2_ij, A2_ji, A3_ij, A3_ji;
  real explp1, p_boc1, p_boc2;
  single_body_parameters *sbp_i, *sbp_j;
  two_body_parameters *twbp;
  bond_order_data *bo_ij, *bo_ji;
  reax_list *bonds = (*lists) + BONDS;
#ifdef TEST_FORCES
  int  k, pk, start_j, end_j;
  int  top_dbo, top_dDelta;
  dbond_data *pdbo;
  dDelta_data *ptop_dDelta;
  reax_list *dDeltas, *dBOs;
  top_dbo=0;
  top_dDelta=0;
  dDeltas = (*lists) + DDELTAS;
  dBOs = (*lists) + DBOS;

  //for( i = 0; i < system->N; ++i )
  //  qsort( &(bonds->select.bond_list[Start_Index(i, bonds)]),
  //   Num_Entries(i, bonds), sizeof(bond_data), compare_bonds );
#endif

  num_bonds = 0;
  p_boc1 = system->reax_param.gp.l[0];
  p_boc2 = system->reax_param.gp.l[1];

 /* Calculate Deltaprime, Deltaprime_boc values */
  for( i = 0; i < system->N; ++i ) {
    type_i = system->my_atoms[i].type;
    sbp_i = &(system->reax_param.sbp[type_i]);
    workspace->Deltap[i] = workspace->total_bond_order[i] - sbp_i->valency;
    workspace->Deltap_boc[i] =
      workspace->total_bond_order[i] - sbp_i->valency_val;

    //fprintf( stdout, "%d(%d) %24.15f\n",
    //     i, workspace->bond_mark[i], workspace->total_bond_order[i] );
    workspace->total_bond_order[i] = 0;
  }
  // fprintf( stderr, "done with uncorrected bond orders\n" );

  /* Corrected Bond Order calculations */
  for( i = 0; i < system->N; ++i ) {
    type_i = system->my_atoms[i].type;
    sbp_i = &(system->reax_param.sbp[type_i]);
    val_i = sbp_i->valency;
    Deltap_i = workspace->Deltap[i];
    Deltap_boc_i = workspace->Deltap_boc[i];
    start_i = Start_Index(i, bonds);
    end_i = End_Index(i, bonds);
    // fprintf( stderr, "i:%d Dp:%g Dbocp:%g s:%d e:%d\n",
    //       i+1, Deltap_i, Deltap_boc_i, start_i, end_i );

    for( pj = start_i; pj < end_i; ++pj ) {
      j = bonds->select.bond_list[pj].nbr;
      type_j = system->my_atoms[j].type;
      bo_ij = &( bonds->select.bond_list[pj].bo_data );
      // fprintf( stderr, "\tj:%d - ubo: %8.3f\n", j+1, bo_ij->BO );

      if( i < j || workspace->bond_mark[j] > 3 ) {
        twbp = &( system->reax_param.tbp[type_i][type_j] );
#ifdef TEST_FORCES
        Set_Start_Index( pj, top_dbo, dBOs );
        /* fprintf( stderr, "%6d%6d%12.6f%12.6f%12.6f\n",
           workspace->reverse_map[i], workspace->reverse_map[j],
           twbp->ovc, twbp->v13cor, bo_ij->BO ); */
#endif

        if( twbp->ovc < 0.001 && twbp->v13cor < 0.001 ) {
          /* There is no correction to bond orders nor to derivatives
             of bond order prime! So we leave bond orders unchanged and
             set derivative of bond order coefficients such that
             dBO = dBOp & dBOxx = dBOxxp in Add_dBO_to_Forces */
          bo_ij->C1dbo = 1.000000;
          bo_ij->C2dbo = 0.000000;
          bo_ij->C3dbo = 0.000000;

          bo_ij->C1dbopi = bo_ij->BO_pi;
          bo_ij->C2dbopi = 0.000000;
          bo_ij->C3dbopi = 0.000000;
          bo_ij->C4dbopi = 0.000000;

          bo_ij->C1dbopi2 = bo_ij->BO_pi2;
          bo_ij->C2dbopi2 = 0.000000;
          bo_ij->C3dbopi2 = 0.000000;
          bo_ij->C4dbopi2 = 0.000000;

#ifdef TEST_FORCES
          pdbo = &(dBOs->select.dbo_list[ top_dbo ]);

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
        else {
          val_j = system->reax_param.sbp[type_j].valency;
          Deltap_j = workspace->Deltap[j];
          Deltap_boc_j = workspace->Deltap_boc[j];

          /* on page 1 */
          if( twbp->ovc >= 0.001 ) {
            /* Correction for overcoordination */
            exp_p1i = exp( -p_boc1 * Deltap_i );
            exp_p2i = exp( -p_boc2 * Deltap_i );
            exp_p1j = exp( -p_boc1 * Deltap_j );
            exp_p2j = exp( -p_boc2 * Deltap_j );

            f2 = exp_p1i + exp_p1j;
            f3 = -1.0 / p_boc2 * log( 0.5 * ( exp_p2i  + exp_p2j ) );
            f1 = 0.5 * ( ( val_i + f2 )/( val_i + f2 + f3 ) +
                         ( val_j + f2 )/( val_j + f2 + f3 ) );

            /*fprintf( stderr,"%d %d\t%g %g   j:%g %g  p_boc:%g %g\n"
              "\tf:%g  %g  %g, exp:%g %g %g %g\n",
              i+1, j+1,
              val_i, Deltap_i, val_j, Deltap_j, p_boc1, p_boc2,
              f1, f2, f3, exp_p1i, exp_p2i, exp_p1j, exp_p2j );*/

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
          else {
            /* No overcoordination correction! */
            f1 = 1.0;
            Cf1_ij = Cf1_ji = 0.0;
          }

          if( twbp->v13cor >= 0.001 ) {
            /* Correction for 1-3 bond orders */
            exp_f4 =exp(-(twbp->p_boc4 * SQR( bo_ij->BO ) -
                          Deltap_boc_i) * twbp->p_boc3 + twbp->p_boc5);
            exp_f5 =exp(-(twbp->p_boc4 * SQR( bo_ij->BO ) -
                          Deltap_boc_j) * twbp->p_boc3 + twbp->p_boc5);

            f4 = 1. / (1. + exp_f4);
            f5 = 1. / (1. + exp_f5);
            f4f5 = f4 * f5;

            /* Bond Order pages 8-9, derivative of f4 and f5 */
            /*temp = twbp->p_boc5 -
              twbp->p_boc3 * twbp->p_boc4 * SQR( bo_ij->BO );
              u_ij = temp + twbp->p_boc3 * Deltap_boc_i;
              u_ji = temp + twbp->p_boc3 * Deltap_boc_j;
              Cf45_ij = Cf45( u_ij, u_ji ) / f4f5;
              Cf45_ji = Cf45( u_ji, u_ij ) / f4f5;*/
            Cf45_ij = -f4 * exp_f4;
            Cf45_ji = -f5 * exp_f5;
          }
          else {
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
          bo_ij->BO    = bo_ij->BO    * A0_ij;
          bo_ij->BO_pi = bo_ij->BO_pi * A0_ij *f1;
          bo_ij->BO_pi2= bo_ij->BO_pi2* A0_ij *f1;
          bo_ij->BO_s  = bo_ij->BO - ( bo_ij->BO_pi + bo_ij->BO_pi2 );

          bo_ij->C1dbo = A0_ij + bo_ij->BO * A1_ij;
          bo_ij->C2dbo = bo_ij->BO * A2_ij;
          bo_ij->C3dbo = bo_ij->BO * A2_ji;

          bo_ij->C1dbopi = f1*f1*f4*f5;
          bo_ij->C2dbopi = bo_ij->BO_pi * A1_ij;
          bo_ij->C3dbopi = bo_ij->BO_pi * A3_ij;
          bo_ij->C4dbopi = bo_ij->BO_pi * A3_ji;

          bo_ij->C1dbopi2 = f1*f1*f4*f5;
          bo_ij->C2dbopi2 = bo_ij->BO_pi2 * A1_ij;
          bo_ij->C3dbopi2 = bo_ij->BO_pi2 * A3_ij;
          bo_ij->C4dbopi2 = bo_ij->BO_pi2 * A3_ji;

#ifdef TEST_FORCES
          Calculate_dBO( i, pj, workspace, lists, &top_dbo );
#endif
        }

        /* neglect bonds that are < 1e-10 */
        if( bo_ij->BO < 1e-10 )
          bo_ij->BO = 0.0;
        if( bo_ij->BO_s < 1e-10 )
          bo_ij->BO_s = 0.0;
        if( bo_ij->BO_pi < 1e-10 )
          bo_ij->BO_pi = 0.0;
        if( bo_ij->BO_pi2 < 1e-10 )
          bo_ij->BO_pi2 = 0.0;

        workspace->total_bond_order[i] += bo_ij->BO; //now keeps total_BO


        /* fprintf( stderr, "%d %d\t%g %g %g %g\n"
           "Cdbo:\t%g %g %g\n"
           "Cdbopi:\t%g %g %g %g\n"
           "Cdbopi2:%g %g %g %g\n\n",
           i+1, j+1,
           bonds->select.bond_list[ pj ].d,
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
      }
      else {
        /* We only need to update bond orders from bo_ji
           everything else is set in uncorrected_bo calculations */
        sym_index = bonds->select.bond_list[pj].sym_index;
        bo_ji = &(bonds->select.bond_list[ sym_index ].bo_data);
        bo_ij->BO = bo_ji->BO;
        bo_ij->BO_s = bo_ji->BO_s;
        bo_ij->BO_pi = bo_ji->BO_pi;
        bo_ij->BO_pi2 = bo_ji->BO_pi2;

        workspace->total_bond_order[i] += bo_ij->BO;// now keeps total_BO
#ifdef TEST_FORCES
        Add_dBO( system, lists, j, sym_index, 1.0, workspace->dDelta );
#endif
      }
    }

#ifdef TEST_FORCES
    // fprintf( stderr, "dDelta computations\nj:" );
    Set_Start_Index( i, top_dDelta, dDeltas );
    ptop_dDelta = &( dDeltas->select.dDelta_list[top_dDelta] );

    for( pj = start_i; pj < end_i; ++pj ) {
      j = bonds->select.bond_list[pj].nbr;
      // fprintf( stderr, "%d  ", j );

      if( !rvec_isZero( workspace->dDelta[j] ) ) {
        ptop_dDelta->wrt = j;
        rvec_Copy( ptop_dDelta->dVal, workspace->dDelta[j] );
        rvec_MakeZero( workspace->dDelta[j] );
        ++top_dDelta, ++ptop_dDelta;
      }

      start_j = Start_Index(j, bonds);
      end_j = End_Index(j, bonds);
      for( pk = start_j; pk < end_j; ++pk ) {
        k = bonds->select.bond_list[pk].nbr;
        if( !rvec_isZero( workspace->dDelta[k] ) ) {
          ptop_dDelta->wrt = k;
          rvec_Copy( ptop_dDelta->dVal, workspace->dDelta[k] );
          rvec_MakeZero( workspace->dDelta[k] );
          ++top_dDelta, ++ptop_dDelta;
        }
      }
    }

    Set_End_Index( i, top_dDelta, dDeltas );

    /*for(pj = Start_Index(i,dDeltas); pj < End_Index(i,dDeltas); ++pj)
      fprintf( stdout, "dDel: %d %d [%g %g %g]\n",
      i+1, dDeltas->select.dDelta_list[pj].wrt+1,
      dDeltas->select.dDelta_list[pj].dVal[0],
      dDeltas->select.dDelta_list[pj].dVal[1],
      dDeltas->select.dDelta_list[pj].dVal[2] );*/
#endif
  }

  /* fprintf( stderr, "\tCalculated actual bond orders ...\n" );
     fprintf( stderr, "j\tDelta\tDelta_e\tDelta_boc\tnlp"
     "\tDelta_lp\tClp\tdDelta_lp\n" );
     fprintf( stderr, "Atom\tDelta\t\tDelta_e\t\tDelta_boc\tnlp"
     "\t\tnlp_opt\t\tDelta_lp\tClp\t\tdDelta_lp\n" );*/

  p_lp1 = system->reax_param.gp.l[15];
  /* Calculate some helper variables that are  used at many places
     throughout force calculations */
  for( j = 0; j < system->N; ++j ){
    type_j = system->my_atoms[j].type;
    sbp_j = &(system->reax_param.sbp[ type_j ]);

    workspace->Delta[j] = workspace->total_bond_order[j] - sbp_j->valency;
    workspace->Delta_e[j] = workspace->total_bond_order[j] - sbp_j->valency_e;
    workspace->Delta_boc[j] = workspace->total_bond_order[j] -
      sbp_j->valency_boc;
    workspace->Delta_val[j] = workspace->total_bond_order[j] -
      sbp_j->valency_val;

    workspace->vlpex[j] = workspace->Delta_e[j] -
      2.0 * (int)(workspace->Delta_e[j]/2.0);
    explp1 = exp(-p_lp1 * SQR(2.0 + workspace->vlpex[j]));
    workspace->nlp[j] = explp1 - (int)(workspace->Delta_e[j] / 2.0);
    workspace->Delta_lp[j] = sbp_j->nlp_opt - workspace->nlp[j];
    workspace->Clp[j] = 2.0 * p_lp1 * explp1 * (2.0 + workspace->vlpex[j]);
    /* Adri uses different dDelta_lp values than the ones in notes... */
    workspace->dDelta_lp[j] = workspace->Clp[j];
    //workspace->dDelta_lp[j] = workspace->Clp[j] + (0.5-workspace->Clp[j]) *
    //((fabs(workspace->Delta_e[j]/2.0 -
    //       (int)(workspace->Delta_e[j]/2.0)) < 0.1) ? 1 : 0 );

    if( sbp_j->mass > 21.0 ) {
      workspace->nlp_temp[j] = 0.5 * (sbp_j->valency_e - sbp_j->valency);
      workspace->Delta_lp_temp[j] = sbp_j->nlp_opt - workspace->nlp_temp[j];
      workspace->dDelta_lp_temp[j] = 0.;
    }
    else {
      workspace->nlp_temp[j] = workspace->nlp[j];
      workspace->Delta_lp_temp[j] = sbp_j->nlp_opt - workspace->nlp_temp[j];
      workspace->dDelta_lp_temp[j] = workspace->Clp[j];
    }

    //fprintf( stderr, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
    //j, workspace->Delta[j], workspace->Delta_e[j], workspace->Delta_boc[j],
    //workspace->nlp[j], system->reaxprm.sbp[type_j].nlp_opt,
    //workspace->Delta_lp[j], workspace->Clp[j], workspace->dDelta_lp[j] );

    //fprintf( stdout, "%d(%d) %24.15f\n",
    //     j, workspace->bond_mark[j], workspace->total_bond_order[j] );
  }

  //Print_Bonds( system, bonds, "pbonds.out" );

#if defined(TEST_ENERGIES) || defined(TEST_FORCES)
  fprintf( stderr, "Number of bonds: %d\n", num_bonds );
  Print_Bond_List( system, control, data, lists, out_control);
#endif
}


#if defined(DEPRECATED)
/* Locate j on i's list.
   This function assumes that j is there for sure!
   And this is the case given our method of neighbor generation*/
int Locate_Symmetric_Bond( reax_list *bonds, int i, int j )
{
  int start = Start_Index(i, bonds);
  int end = End_Index(i, bonds);
  int mid = (start + end) / 2;
  int mid_nbr;

  while( (mid_nbr = bonds->select.bond_list[mid].nbr) != j ) {
    /*fprintf( stderr, "\tstart: %d   end: %d   mid: %d\n",
      start, end, mid );*/
    if( mid_nbr < j )
      start = mid+1;
    else end = mid - 1;

    mid = (start + end) / 2;
  }

  return mid;
}


inline void Copy_Bond_Order_Data( bond_order_data *dest, bond_order_data *src )
{
  dest->BO = src->BO;
  dest->BO_s = src->BO_s;
  dest->BO_pi = src->BO_pi;
  dest->BO_pi2 = src->BO_pi2;

  rvec_Scale( dest->dBOp, -1.0, src->dBOp );
  rvec_Scale( dest->dln_BOp_s, -1.0, src->dln_BOp_s );
  rvec_Scale( dest->dln_BOp_pi, -1.0, src->dln_BOp_pi );
  rvec_Scale( dest->dln_BOp_pi2, -1.0, src->dln_BOp_pi2 );
}
#endif
