/*----------------------------------------------------------------------
  PuReMD - Purdue ReaxFF Molecular Dynamics Program
  Website: https://www.cs.purdue.edu/puremd

  Copyright (2010) Purdue University

  Contributing authors:
  H. M. Aktulga, J. Fogarty, S. Pandit, A. Grama
  Corresponding author:
  Hasan Metin Aktulga, Michigan State University, hma@cse.msu.edu

  Please cite the related publication:
  H. M. Aktulga, J. C. Fogarty, S. A. Pandit, A. Y. Grama,
  "Parallel Reactive Molecular Dynamics: Numerical Methods and
  Algorithmic Techniques", Parallel Computing, 38 (4-5), 245-259

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

#include "reaxc_bond_orders_omp.h"
#include <mpi.h>
#include <cmath>
#include "fix_omp.h"
#include "reaxc_defs.h"
#include "pair_reaxc_omp.h"
#include "reaxc_types.h"
#include "reaxc_list.h"
#include "reaxc_vector.h"

#if defined(_OPENMP)
#include  <omp.h>
#endif

using namespace LAMMPS_NS;

void Add_dBond_to_ForcesOMP( reax_system *system, int i, int pj,
                             storage *workspace, reax_list **lists ) {
  reax_list *bonds = (*lists) + BONDS;
  bond_data *nbr_j, *nbr_k;
  bond_order_data *bo_ij, *bo_ji;
  dbond_coefficients coef;
  int pk, k, j;

  PairReaxCOMP *pair_reax_ptr = static_cast<class PairReaxCOMP*>(system->pair_ptr);

#if defined(_OPENMP)
  int tid = omp_get_thread_num();
#else
  int tid = 0;
#endif
  ThrData *thr = pair_reax_ptr->getFixOMP()->get_thr(tid);
  long reductionOffset = (system->N * tid);

  /* Virial Tallying variables */
  rvec fi_tmp, fj_tmp, fk_tmp, delij, delji, delki, delkj, temp;

  /* Initializations */
  nbr_j = &(bonds->select.bond_list[pj]);
  j = nbr_j->nbr;
  bo_ij = &(nbr_j->bo_data);
  bo_ji = &(bonds->select.bond_list[ nbr_j->sym_index ].bo_data);

  double c = bo_ij->Cdbo + bo_ji->Cdbo;
  coef.C1dbo = bo_ij->C1dbo * c;
  coef.C2dbo = bo_ij->C2dbo * c;
  coef.C3dbo = bo_ij->C3dbo * c;

  c = bo_ij->Cdbopi + bo_ji->Cdbopi;
  coef.C1dbopi = bo_ij->C1dbopi * c;
  coef.C2dbopi = bo_ij->C2dbopi * c;
  coef.C3dbopi = bo_ij->C3dbopi * c;
  coef.C4dbopi = bo_ij->C4dbopi * c;

  c = bo_ij->Cdbopi2 + bo_ji->Cdbopi2;
  coef.C1dbopi2 = bo_ij->C1dbopi2 * c;
  coef.C2dbopi2 = bo_ij->C2dbopi2 * c;
  coef.C3dbopi2 = bo_ij->C3dbopi2 * c;
  coef.C4dbopi2 = bo_ij->C4dbopi2 * c;

  c = workspace->CdDelta[i] + workspace->CdDelta[j];
  coef.C1dDelta = bo_ij->C1dbo * c;
  coef.C2dDelta = bo_ij->C2dbo * c;
  coef.C3dDelta = bo_ij->C3dbo * c;

  // The same "c" refactoring here can be replicated below in Add_dBond_to_Forces_NPTOMP(), but
  // I'd prefer to wait for a test to verify changes before doing so (just to be safe).

  // forces on i
  // rvec_Scale(     temp, coef.C1dbo,    bo_ij->dBOp );
  // rvec_ScaledAdd( temp, coef.C2dbo,    workspace->dDeltap_self[i] );
  // rvec_ScaledAdd( temp, coef.C1dDelta, bo_ij->dBOp );
  // rvec_ScaledAdd( temp, coef.C2dDelta, workspace->dDeltap_self[i] );
  // rvec_ScaledAdd( temp, coef.C1dbopi,  bo_ij->dln_BOp_pi );
  // rvec_ScaledAdd( temp, coef.C2dbopi,  bo_ij->dBOp );
  // rvec_ScaledAdd( temp, coef.C3dbopi,  workspace->dDeltap_self[i]);
  // rvec_ScaledAdd( temp, coef.C1dbopi2, bo_ij->dln_BOp_pi2 );
  // rvec_ScaledAdd( temp, coef.C2dbopi2, bo_ij->dBOp );
  // rvec_ScaledAdd( temp, coef.C3dbopi2, workspace->dDeltap_self[i] );

  c = (coef.C1dbo + coef.C1dDelta + coef.C2dbopi + coef.C2dbopi2);
  rvec_Scale(     temp, c,    bo_ij->dBOp );

  c = (coef.C2dbo + coef.C2dDelta + coef.C3dbopi + coef.C3dbopi2);
  rvec_ScaledAdd( temp, c,    workspace->dDeltap_self[i] );

  rvec_ScaledAdd( temp, coef.C1dbopi,  bo_ij->dln_BOp_pi );
  rvec_ScaledAdd( temp, coef.C1dbopi2, bo_ij->dln_BOp_pi2 );

  rvec_Add(workspace->forceReduction[reductionOffset+i],temp );

  if (system->pair_ptr->vflag_atom) {
    rvec_Scale(fi_tmp, -1.0, temp);
    rvec_ScaledSum( delij, 1., system->my_atoms[i].x,-1., system->my_atoms[j].x );

    pair_reax_ptr->ev_tally_xyz_thr_proxy(system->pair_ptr,i,j,system->N,0,0,0,
                                          fi_tmp[0],fi_tmp[1],fi_tmp[2],
                                          delij[0],delij[1],delij[2],thr);
  }

  // forces on j
  // rvec_Scale(     temp, -coef.C1dbo,    bo_ij->dBOp );
  // rvec_ScaledAdd( temp,  coef.C3dbo,    workspace->dDeltap_self[j] );
  // rvec_ScaledAdd( temp, -coef.C1dDelta, bo_ij->dBOp );
  // rvec_ScaledAdd( temp,  coef.C3dDelta, workspace->dDeltap_self[j]);
  // rvec_ScaledAdd( temp, -coef.C1dbopi,  bo_ij->dln_BOp_pi );
  // rvec_ScaledAdd( temp, -coef.C2dbopi,  bo_ij->dBOp );
  // rvec_ScaledAdd( temp,  coef.C4dbopi,  workspace->dDeltap_self[j]);
  // rvec_ScaledAdd( temp, -coef.C1dbopi2, bo_ij->dln_BOp_pi2 );
  // rvec_ScaledAdd( temp, -coef.C2dbopi2, bo_ij->dBOp );
  // rvec_ScaledAdd( temp,  coef.C4dbopi2, workspace->dDeltap_self[j]);


  c = -(coef.C1dbo + coef.C1dDelta + coef.C2dbopi + coef.C2dbopi2);
  rvec_Scale(     temp, c,    bo_ij->dBOp );

  c = (coef.C3dbo + coef.C3dDelta + coef.C4dbopi + coef.C4dbopi2);
  rvec_ScaledAdd( temp,  c,    workspace->dDeltap_self[j] );

  rvec_ScaledAdd( temp, -coef.C1dbopi,  bo_ij->dln_BOp_pi );
  rvec_ScaledAdd( temp, -coef.C1dbopi2, bo_ij->dln_BOp_pi2 );


  rvec_Add(workspace->forceReduction[reductionOffset+j],temp );

  if (system->pair_ptr->vflag_atom) {
    rvec_Scale(fj_tmp, -1.0, temp);
    rvec_ScaledSum( delji, 1., system->my_atoms[j].x,-1., system->my_atoms[i].x );

    pair_reax_ptr->ev_tally_xyz_thr_proxy(system->pair_ptr,j,i,system->N,0,0,0,
                                          fj_tmp[0],fj_tmp[1],fj_tmp[2],
                                          delji[0],delji[1],delji[2],thr);
  }

  // forces on k: i neighbor
  for( pk = Start_Index(i, bonds); pk < End_Index(i, bonds); ++pk ) {
    nbr_k = &(bonds->select.bond_list[pk]);
    k = nbr_k->nbr;

    // rvec_Scale(     temp, -coef.C2dbo,    nbr_k->bo_data.dBOp);
    // rvec_ScaledAdd( temp, -coef.C2dDelta, nbr_k->bo_data.dBOp);
    // rvec_ScaledAdd( temp, -coef.C3dbopi,  nbr_k->bo_data.dBOp);
    // rvec_ScaledAdd( temp, -coef.C3dbopi2, nbr_k->bo_data.dBOp);

    const double c = -(coef.C2dbo + coef.C2dDelta + coef.C3dbopi + coef.C3dbopi2);
    rvec_Scale(temp, c, nbr_k->bo_data.dBOp);

    rvec_Add(workspace->forceReduction[reductionOffset+k],temp );

    if (system->pair_ptr->vflag_atom) {
      rvec_Scale(fk_tmp, -1.0, temp);
      rvec_ScaledSum(delki,1.,system->my_atoms[k].x,-1.,system->my_atoms[i].x);

      pair_reax_ptr->ev_tally_xyz_thr_proxy(system->pair_ptr,k,i,system->N,0,0,0,
                                            fk_tmp[0],fk_tmp[1],fk_tmp[2],
                                            delki[0],delki[1],delki[2],thr);
      rvec_ScaledSum(delkj,1.,system->my_atoms[k].x,-1.,system->my_atoms[j].x);

      pair_reax_ptr->ev_tally_xyz_thr_proxy(system->pair_ptr,k,j,system->N,0,0,0,
                                            fk_tmp[0],fk_tmp[1],fk_tmp[2],
                                            delkj[0],delkj[1],delkj[2],thr);
    }
  }

  // forces on k: j neighbor
  for( pk = Start_Index(j, bonds); pk < End_Index(j, bonds); ++pk ) {
    nbr_k = &(bonds->select.bond_list[pk]);
    k = nbr_k->nbr;

    // rvec_Scale(     temp, -coef.C3dbo,    nbr_k->bo_data.dBOp );
    // rvec_ScaledAdd( temp, -coef.C3dDelta, nbr_k->bo_data.dBOp);
    // rvec_ScaledAdd( temp, -coef.C4dbopi,  nbr_k->bo_data.dBOp);
    // rvec_ScaledAdd( temp, -coef.C4dbopi2, nbr_k->bo_data.dBOp);

    const double c = -(coef.C3dbo + coef.C3dDelta + coef.C4dbopi + coef.C4dbopi2);
    rvec_Scale(temp, c, nbr_k->bo_data.dBOp);

    rvec_Add(workspace->forceReduction[reductionOffset+k],temp );

    if (system->pair_ptr->vflag_atom) {
      rvec_Scale(fk_tmp, -1.0, temp);
      rvec_ScaledSum(delki,1.,system->my_atoms[k].x,-1.,system->my_atoms[i].x);

      pair_reax_ptr->ev_tally_xyz_thr_proxy(system->pair_ptr,k,i,system->N,0,0,0,
                                            fk_tmp[0],fk_tmp[1],fk_tmp[2],
                                            delki[0],delki[1],delki[2],thr);

      rvec_ScaledSum(delkj,1.,system->my_atoms[k].x,-1.,system->my_atoms[j].x);

      pair_reax_ptr->ev_tally_xyz_thr_proxy(system->pair_ptr,k,j,system->N,0,0,0,
                                            fk_tmp[0],fk_tmp[1],fk_tmp[2],
                                            delkj[0],delkj[1],delkj[2],thr);
    }
  }
}

/* ---------------------------------------------------------------------- */

void Add_dBond_to_Forces_NPTOMP( reax_system *system, int i, int pj,
                                 simulation_data * /* data */,
                                 storage *workspace, reax_list **lists ) {
  reax_list *bonds = (*lists) + BONDS;
  bond_data *nbr_j, *nbr_k;
  bond_order_data *bo_ij, *bo_ji;
  dbond_coefficients coef;
  rvec temp, ext_press;
  ivec rel_box;
  int pk, k, j;

#if defined(_OPENMP)
  int tid = omp_get_thread_num();
#else
  int tid = 0;
#endif
  long reductionOffset = (system->N * tid);

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
    rvec_Add(workspace->forceReduction[reductionOffset+k],temp );

    /* pressure */
    rvec_iMultiply( ext_press, nbr_k->rel_box, temp );
    rvec_Add( workspace->my_ext_pressReduction[tid], ext_press );
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
  rvec_Add(workspace->forceReduction[reductionOffset+i],temp );

  for( pk = Start_Index(j, bonds); pk < End_Index(j, bonds); ++pk ) {
    nbr_k = &(bonds->select.bond_list[pk]);
    k = nbr_k->nbr;

    rvec_Scale( temp, -coef.C3dbo, nbr_k->bo_data.dBOp );      /*3rd,dBO*/
    rvec_ScaledAdd( temp, -coef.C3dDelta, nbr_k->bo_data.dBOp);/*dDelta*/
    rvec_ScaledAdd( temp, -coef.C4dbopi, nbr_k->bo_data.dBOp); /*4th,dBOpi*/
    rvec_ScaledAdd( temp, -coef.C4dbopi2, nbr_k->bo_data.dBOp);/*4th,dBOpi2*/

    /* force */
    rvec_Add(workspace->forceReduction[reductionOffset+k],temp );

    /* pressure */
    if (k != i) {
      ivec_Sum( rel_box, nbr_k->rel_box, nbr_j->rel_box ); //rel_box(k, i)
      rvec_iMultiply( ext_press, rel_box, temp );
      rvec_Add( workspace->my_ext_pressReduction[tid], ext_press );
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
  rvec_Add(workspace->forceReduction[reductionOffset+j],temp );

  /* pressure */
  rvec_iMultiply( ext_press, nbr_j->rel_box, temp );
    rvec_Add( workspace->my_ext_pressReduction[tid], ext_press );
}

/* ---------------------------------------------------------------------- */

int BOp_OMP( storage * /* workspace */, reax_list *bonds, double bo_cut,
             int i, int btop_i, far_neighbor_data *nbr_pj,
             single_body_parameters * /* sbp_i */, single_body_parameters * /* sbp_j */,
             two_body_parameters *twbp,
             int btop_j, double C12, double C34, double C56, double BO, double BO_s, double BO_pi, double BO_pi2) {
  int j;
  double rr2;
  double Cln_BOp_s, Cln_BOp_pi, Cln_BOp_pi2;
  bond_data *ibond, *jbond;
  bond_order_data *bo_ij, *bo_ji;

  j = nbr_pj->nbr;
  rr2 = 1.0 / SQR(nbr_pj->d);

  // Top portion of BOp() moved to reaxc_forces_omp.cpp::Init_Forces_noQEq_OMP()

  /* Initially BO values are the uncorrected ones, page 1 */

  /****** bonds i-j and j-i ******/
  ibond = &( bonds->select.bond_list[btop_i] );
  jbond = &( bonds->select.bond_list[btop_j] );

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

  bo_ij = &( ibond->bo_data );
  bo_ji = &( jbond->bo_data );
  bo_ji->BO     = bo_ij->BO     = BO;
  bo_ji->BO_s   = bo_ij->BO_s   = BO_s;
  bo_ji->BO_pi  = bo_ij->BO_pi  = BO_pi;
  bo_ji->BO_pi2 = bo_ij->BO_pi2 = BO_pi2;

  /* Bond Order page2-3, derivative of total bond order prime */
  Cln_BOp_s   = twbp->p_bo2 * C12 * rr2;
  Cln_BOp_pi  = twbp->p_bo4 * C34 * rr2;
  Cln_BOp_pi2 = twbp->p_bo6 * C56 * rr2;

  /* Only dln_BOp_xx wrt. dr_i is stored here, note that
     dln_BOp_xx/dr_i = -dln_BOp_xx/dr_j and all others are 0 */
  rvec_Scale(bo_ij->dln_BOp_s,-bo_ij->BO_s*Cln_BOp_s,ibond->dvec);
  rvec_Scale(bo_ij->dln_BOp_pi,-bo_ij->BO_pi*Cln_BOp_pi,ibond->dvec);
  rvec_Scale(bo_ij->dln_BOp_pi2,
             -bo_ij->BO_pi2*Cln_BOp_pi2,ibond->dvec);
  rvec_Scale(bo_ji->dln_BOp_s,   -1., bo_ij->dln_BOp_s);
  rvec_Scale(bo_ji->dln_BOp_pi,  -1., bo_ij->dln_BOp_pi );
  rvec_Scale(bo_ji->dln_BOp_pi2, -1., bo_ij->dln_BOp_pi2 );

  rvec_Scale( bo_ij->dBOp,
              -(bo_ij->BO_s * Cln_BOp_s +
                bo_ij->BO_pi * Cln_BOp_pi +
                bo_ij->BO_pi2 * Cln_BOp_pi2), ibond->dvec );
  rvec_Scale( bo_ji->dBOp, -1., bo_ij->dBOp );

  bo_ij->BO_s -= bo_cut;
  bo_ij->BO   -= bo_cut;
  bo_ji->BO_s -= bo_cut;
  bo_ji->BO   -= bo_cut;

  bo_ij->Cdbo = bo_ij->Cdbopi = bo_ij->Cdbopi2 = 0.0;
  bo_ji->Cdbo = bo_ji->Cdbopi = bo_ji->Cdbopi2 = 0.0;

  return 1;
}

/* ---------------------------------------------------------------------- */

void BOOMP( reax_system *system, control_params * /* control */, simulation_data * /* data */,
            storage *workspace, reax_list **lists, output_controls * /* out_control */)
{
#ifdef OMP_TIMING
  double endTimeBase, startTimeBase;
  startTimeBase = MPI_Wtime();
#endif

  double p_lp1 = system->reax_param.gp.l[15];
  double p_boc1 = system->reax_param.gp.l[0];
  double p_boc2 = system->reax_param.gp.l[1];
  reax_list *bonds = (*lists) + BONDS;

#if defined(_OPENMP)
#pragma omp parallel default(shared)
#endif
  {
    int  i, j, pj, type_i, type_j;
    int  start_i, end_i, sym_index;
    double val_i, Deltap_i, Deltap_boc_i;
    double val_j, Deltap_j, Deltap_boc_j;
    double f1, f2, f3, f4, f5, f4f5, exp_f4, exp_f5;
    double exp_p1i, exp_p2i, exp_p1j, exp_p2j, explp1;
    double temp, u1_ij, u1_ji, Cf1A_ij, Cf1B_ij, Cf1_ij, Cf1_ji;
    double Cf45_ij, Cf45_ji; //u_ij, u_ji
    double A0_ij, A1_ij, A2_ij, A2_ji, A3_ij, A3_ji;
    single_body_parameters *sbp_i, *sbp_j;
    two_body_parameters *twbp;
    bond_order_data *bo_ij, *bo_ji;

    /* Calculate Deltaprime, Deltaprime_boc values */
#if defined(_OPENMP)
#pragma omp for schedule(static)
#endif
    for (i = 0; i < system->N; ++i) {
      type_i = system->my_atoms[i].type;
      if(type_i < 0) continue;
      sbp_i = &(system->reax_param.sbp[type_i]);
      workspace->Deltap[i] = workspace->total_bond_order[i] - sbp_i->valency;
      workspace->Deltap_boc[i] =
        workspace->total_bond_order[i] - sbp_i->valency_val;

      workspace->total_bond_order[i] = 0;
    }

    // Wait till initialization complete
#if defined(_OPENMP)
#pragma omp barrier
#endif

    /* Corrected Bond Order calculations */
#if defined(_OPENMP)
#pragma omp for schedule(guided)
#endif
    for (i = 0; i < system->N; ++i) {
      type_i = system->my_atoms[i].type;
      if(type_i < 0) continue;
      sbp_i = &(system->reax_param.sbp[type_i]);
      val_i = sbp_i->valency;
      Deltap_i = workspace->Deltap[i];
      Deltap_boc_i = workspace->Deltap_boc[i];
      start_i = Start_Index(i, bonds);
      end_i = End_Index(i, bonds);

      for (pj = start_i; pj < end_i; ++pj) {
        j = bonds->select.bond_list[pj].nbr;
        type_j = system->my_atoms[j].type;
        if(type_j < 0) continue;
        bo_ij = &( bonds->select.bond_list[pj].bo_data );

        if (i < j || workspace->bond_mark[j] > 3) {
          twbp = &( system->reax_param.tbp[type_i][type_j] );

          if (twbp->ovc < 0.001 && twbp->v13cor < 0.001) {
            bo_ij->C1dbo = 1.000000;
            bo_ij->C2dbo = 0.000000;
            bo_ij->C3dbo = 0.000000;

            bo_ij->C1dbopi = 1.000000;
            bo_ij->C2dbopi = 0.000000;
            bo_ij->C3dbopi = 0.000000;
            bo_ij->C4dbopi = 0.000000;

            bo_ij->C1dbopi2 = 1.000000;
            bo_ij->C2dbopi2 = 0.000000;
            bo_ij->C3dbopi2 = 0.000000;
            bo_ij->C4dbopi2 = 0.000000;

          }
          else {
            val_j = system->reax_param.sbp[type_j].valency;
            Deltap_j = workspace->Deltap[j];
            Deltap_boc_j = workspace->Deltap_boc[j];

            /* on page 1 */
            if (twbp->ovc >= 0.001) {
              /* Correction for overcoordination */
              exp_p1i = exp( -p_boc1 * Deltap_i );
              exp_p2i = exp( -p_boc2 * Deltap_i );
              exp_p1j = exp( -p_boc1 * Deltap_j );
              exp_p2j = exp( -p_boc2 * Deltap_j );

              f2 = exp_p1i + exp_p1j;
              f3 = -1.0 / p_boc2 * log( 0.5 * ( exp_p2i  + exp_p2j ) );
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
            }
            else {
              /* No overcoordination correction! */
              f1 = 1.0;
              Cf1_ij = Cf1_ji = 0.0;
            }

            if (twbp->v13cor >= 0.001) {
              /* Correction for 1-3 bond orders */
              exp_f4 =exp(-(twbp->p_boc4 * SQR( bo_ij->BO ) -
                            Deltap_boc_i) * twbp->p_boc3 + twbp->p_boc5);
              exp_f5 =exp(-(twbp->p_boc4 * SQR( bo_ij->BO ) -
                            Deltap_boc_j) * twbp->p_boc3 + twbp->p_boc5);

              f4 = 1. / (1. + exp_f4);
              f5 = 1. / (1. + exp_f5);
              f4f5 = f4 * f5;

              /* Bond Order pages 8-9, derivative of f4 and f5 */
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
          }

          /* neglect bonds that are < 1e-10 */
          if (bo_ij->BO < 1e-10)
            bo_ij->BO = 0.0;
          if (bo_ij->BO_s < 1e-10)
            bo_ij->BO_s = 0.0;
          if (bo_ij->BO_pi < 1e-10)
            bo_ij->BO_pi = 0.0;
          if (bo_ij->BO_pi2 < 1e-10)
            bo_ij->BO_pi2 = 0.0;

          workspace->total_bond_order[i] += bo_ij->BO; //now keeps total_BO
        }
        // else {
        //   /* We only need to update bond orders from bo_ji
        //      everything else is set in uncorrected_bo calculations */
        //   sym_index = bonds->select.bond_list[pj].sym_index;
        //   bo_ji = &(bonds->select.bond_list[ sym_index ].bo_data);
        //   bo_ij->BO = bo_ji->BO;
        //   bo_ij->BO_s = bo_ji->BO_s;
        //   bo_ij->BO_pi = bo_ji->BO_pi;
        //   bo_ij->BO_pi2 = bo_ji->BO_pi2;

        //   workspace->total_bond_order[i] += bo_ij->BO;// now keeps total_BO
        // }
      }

    }

    // Wait for bo_ij to be updated
#if defined(_OPENMP)
#pragma omp barrier
#endif
    // Try to combine the following for-loop back into the for-loop above
    /*-------------------------*/
#if defined(_OPENMP)
#pragma omp for schedule(guided)
#endif
    for (i = 0; i < system->N; ++i) {
      type_i = system->my_atoms[i].type;
      if(type_i < 0) continue;
      start_i = Start_Index(i, bonds);
      end_i = End_Index(i, bonds);

      for (pj = start_i; pj < end_i; ++pj) {
        j = bonds->select.bond_list[pj].nbr;
        type_j = system->my_atoms[j].type;
        if(type_j < 0) continue;

        if (i < j || workspace->bond_mark[j] > 3) {
          // Computed in previous for-loop
        } else {
          /* We only need to update bond orders from bo_ji
             everything else is set in uncorrected_bo calculations */
          sym_index = bonds->select.bond_list[pj].sym_index;

          bo_ij = &( bonds->select.bond_list[pj].bo_data );
          bo_ji = &(bonds->select.bond_list[ sym_index ].bo_data);
          bo_ij->BO = bo_ji->BO;
          bo_ij->BO_s = bo_ji->BO_s;
          bo_ij->BO_pi = bo_ji->BO_pi;
          bo_ij->BO_pi2 = bo_ji->BO_pi2;

          workspace->total_bond_order[i] += bo_ij->BO;// now keeps total_BO
        }
      }

    }

    /*-------------------------*/

    // Need to wait for total_bond_order to be accumulated.
#if defined(_OPENMP)
#pragma omp barrier
#endif
    /* Calculate some helper variables that are  used at many places
       throughout force calculations */
#if defined(_OPENMP)
#pragma omp for schedule(guided)
#endif
    for(j = 0; j < system->N; ++j ) {
      type_j = system->my_atoms[j].type;
      if(type_j < 0) continue;
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
      workspace->dDelta_lp[j] = workspace->Clp[j];

      if (sbp_j->mass > 21.0) {
        workspace->nlp_temp[j] = 0.5 * (sbp_j->valency_e - sbp_j->valency);
        workspace->Delta_lp_temp[j] = sbp_j->nlp_opt - workspace->nlp_temp[j];
        workspace->dDelta_lp_temp[j] = 0.;
      }
      else {
        workspace->nlp_temp[j] = workspace->nlp[j];
        workspace->Delta_lp_temp[j] = sbp_j->nlp_opt - workspace->nlp_temp[j];
        workspace->dDelta_lp_temp[j] = workspace->Clp[j];
      }
    }

  } // parallel region

#ifdef OMP_TIMING
        endTimeBase = MPI_Wtime();
        ompTimingData[COMPUTEBOINDEX] += (endTimeBase-startTimeBase);

#endif
}
