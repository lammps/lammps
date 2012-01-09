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
#if defined(PURE_REAX)
#include "torsion_angles.h"
#include "bond_orders.h"
#include "list.h"
#include "tool_box.h"
#include "vector.h"
#elif defined(LAMMPS_REAX)
#include "reaxc_torsion_angles.h"
#include "reaxc_bond_orders.h"
#include "reaxc_list.h"
#include "reaxc_tool_box.h"
#include "reaxc_vector.h"
#endif

#define MIN_SINE 1e-10

real Calculate_Omega( rvec dvec_ij, real r_ij,
		      rvec dvec_jk, real r_jk,
		      rvec dvec_kl, real r_kl,
		      rvec dvec_li, real r_li,
		      three_body_interaction_data *p_ijk, 
		      three_body_interaction_data *p_jkl, 
		      rvec dcos_omega_di, rvec dcos_omega_dj, 
		      rvec dcos_omega_dk, rvec dcos_omega_dl, 
		      output_controls *out_control )
{
  real unnorm_cos_omega, unnorm_sin_omega, omega;
  real sin_ijk, cos_ijk, sin_jkl, cos_jkl;
  real htra, htrb, htrc, hthd, hthe, hnra, hnrc, hnhd, hnhe;
  real arg, poem, tel;
  rvec cross_jk_kl;

  sin_ijk = sin( p_ijk->theta );
  cos_ijk = cos( p_ijk->theta );
  sin_jkl = sin( p_jkl->theta );
  cos_jkl = cos( p_jkl->theta );

  /* omega */
  unnorm_cos_omega = -rvec_Dot(dvec_ij, dvec_jk) * rvec_Dot(dvec_jk, dvec_kl) + 
    SQR( r_jk ) *  rvec_Dot( dvec_ij, dvec_kl );

  rvec_Cross( cross_jk_kl, dvec_jk, dvec_kl );
  unnorm_sin_omega = -r_jk * rvec_Dot( dvec_ij, cross_jk_kl );

  omega = atan2( unnorm_sin_omega, unnorm_cos_omega ); 


  /* derivatives */
  /* coef for adjusments to cos_theta's */
  /* rla = r_ij, rlb = r_jk, rlc = r_kl, r4 = r_li;
     coshd = cos_ijk, coshe = cos_jkl;
     sinhd = sin_ijk, sinhe = sin_jkl; */
  htra = r_ij + cos_ijk * ( r_kl * cos_jkl - r_jk );
  htrb = r_jk - r_ij * cos_ijk - r_kl * cos_jkl;
  htrc = r_kl + cos_jkl * ( r_ij * cos_ijk - r_jk );
  hthd = r_ij * sin_ijk * ( r_jk - r_kl * cos_jkl );
  hthe = r_kl * sin_jkl * ( r_jk - r_ij * cos_ijk );
  hnra = r_kl * sin_ijk * sin_jkl;
  hnrc = r_ij * sin_ijk * sin_jkl;
  hnhd = r_ij * r_kl * cos_ijk * sin_jkl;
  hnhe = r_ij * r_kl * sin_ijk * cos_jkl;


  poem = 2.0 * r_ij * r_kl * sin_ijk * sin_jkl;
  if( poem < 1e-20 ) poem = 1e-20;

  tel  = SQR( r_ij ) + SQR( r_jk ) + SQR( r_kl ) - SQR( r_li ) - 
    2.0 * ( r_ij * r_jk * cos_ijk - r_ij * r_kl * cos_ijk * cos_jkl + 
	    r_jk * r_kl * cos_jkl );
				  
  arg  = tel / poem;
  if( arg >  1.0 ) arg =  1.0;
  if( arg < -1.0 ) arg = -1.0;


  /* fprintf( out_control->etor, 
     "%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f\n",
     htra, htrb, htrc, hthd, hthe, hnra, hnrc, hnhd, hnhe );
     fprintf( out_control->etor, "%12.6f%12.6f%12.6f\n",
     dvec_ij[0]/r_ij, dvec_ij[1]/r_ij, dvec_ij[2]/r_ij );
     fprintf( out_control->etor, "%12.6f%12.6f%12.6f\n",
     -dvec_jk[0]/r_jk, -dvec_jk[1]/r_jk, -dvec_jk[2]/r_jk );
     fprintf( out_control->etor, "%12.6f%12.6f%12.6f\n",
     -dvec_kl[0]/r_kl, -dvec_kl[1]/r_kl, -dvec_kl[2]/r_kl );
     fprintf( out_control->etor, "%12.6f%12.6f%12.6f%12.6f\n",
     r_li, dvec_li[0], dvec_li[1], dvec_li[2] );
     fprintf( out_control->etor, "%12.6f%12.6f%12.6f\n",
     poem, tel, arg ); */
  /* fprintf( out_control->etor, "%12.6f%12.6f%12.6f\n",
     -p_ijk->dcos_dk[0]/sin_ijk, -p_ijk->dcos_dk[1]/sin_ijk, 
     -p_ijk->dcos_dk[2]/sin_ijk );
     fprintf( out_control->etor, "%12.6f%12.6f%12.6f\n",
     -p_jkl->dcos_dk[0]/sin_jkl, -p_jkl->dcos_dk[1]/sin_jkl, 
     -p_jkl->dcos_dk[2]/sin_jkl );*/

  if( sin_ijk >= 0 && sin_ijk <= MIN_SINE ) sin_ijk = MIN_SINE;
  else if( sin_ijk <= 0 && sin_ijk >= -MIN_SINE ) sin_ijk = -MIN_SINE;
  if( sin_jkl >= 0 && sin_jkl <= MIN_SINE ) sin_jkl = MIN_SINE;
  else if( sin_jkl <= 0 && sin_jkl >= -MIN_SINE ) sin_jkl = -MIN_SINE;

  // dcos_omega_di
  rvec_ScaledSum( dcos_omega_di, (htra-arg*hnra)/r_ij, dvec_ij, -1., dvec_li );
  rvec_ScaledAdd( dcos_omega_di,-(hthd-arg*hnhd)/sin_ijk, p_ijk->dcos_dk );
  rvec_Scale( dcos_omega_di, 2.0 / poem, dcos_omega_di );

  // dcos_omega_dj
  rvec_ScaledSum( dcos_omega_dj,-(htra-arg*hnra)/r_ij, dvec_ij, 
		  -htrb / r_jk, dvec_jk );
  rvec_ScaledAdd( dcos_omega_dj,-(hthd-arg*hnhd)/sin_ijk, p_ijk->dcos_dj );
  rvec_ScaledAdd( dcos_omega_dj,-(hthe-arg*hnhe)/sin_jkl, p_jkl->dcos_di );
  rvec_Scale( dcos_omega_dj, 2.0 / poem, dcos_omega_dj );

  // dcos_omega_dk
  rvec_ScaledSum( dcos_omega_dk,-(htrc-arg*hnrc)/r_kl, dvec_kl,  
		  htrb / r_jk, dvec_jk );
  rvec_ScaledAdd( dcos_omega_dk,-(hthd-arg*hnhd)/sin_ijk, p_ijk->dcos_di );
  rvec_ScaledAdd( dcos_omega_dk,-(hthe-arg*hnhe)/sin_jkl, p_jkl->dcos_dj );
  rvec_Scale( dcos_omega_dk, 2.0 / poem, dcos_omega_dk );

  // dcos_omega_dl
  rvec_ScaledSum( dcos_omega_dl, (htrc-arg*hnrc)/r_kl, dvec_kl, 1., dvec_li );
  rvec_ScaledAdd( dcos_omega_dl,-(hthe-arg*hnhe)/sin_jkl, p_jkl->dcos_dk );
  rvec_Scale( dcos_omega_dl, 2.0 / poem, dcos_omega_dl );

  return omega;  
}



void Torsion_Angles( reax_system *system, control_params *control, 
		     simulation_data *data, storage *workspace, 
		     reax_list **lists, output_controls *out_control )
{
  int i, j, k, l, pi, pj, pk, pl, pij, plk, natoms;
  int type_i, type_j, type_k, type_l;
  int start_j, end_j, start_k, end_k;
  int start_pj, end_pj, start_pk, end_pk;
  int num_frb_intrs = 0;

  real Delta_j, Delta_k;
  real r_ij, r_jk, r_kl, r_li;
  real BOA_ij, BOA_jk, BOA_kl;

  real exp_tor2_ij, exp_tor2_jk, exp_tor2_kl;
  real exp_tor1, exp_tor3_DjDk, exp_tor4_DjDk, exp_tor34_inv;
  real exp_cot2_jk, exp_cot2_ij, exp_cot2_kl;
  real fn10, f11_DjDk, dfn11, fn12;
  real theta_ijk, theta_jkl;
  real sin_ijk, sin_jkl;
  real cos_ijk, cos_jkl;
  real tan_ijk_i, tan_jkl_i;
  real omega, cos_omega, cos2omega, cos3omega;
  rvec dcos_omega_di, dcos_omega_dj, dcos_omega_dk, dcos_omega_dl;
  real CV, cmn, CEtors1, CEtors2, CEtors3, CEtors4;
  real CEtors5, CEtors6, CEtors7, CEtors8, CEtors9;
  real Cconj, CEconj1, CEconj2, CEconj3;
  real CEconj4, CEconj5, CEconj6;
  real e_tor, e_con;
  rvec dvec_li;
  rvec force, ext_press;
  ivec rel_box_jl;
  // rtensor total_rtensor, temp_rtensor;
  four_body_header *fbh;
  four_body_parameters *fbp;
  bond_data *pbond_ij, *pbond_jk, *pbond_kl;
  bond_order_data *bo_ij, *bo_jk, *bo_kl;
  three_body_interaction_data *p_ijk, *p_jkl;
  real p_tor2 = system->reax_param.gp.l[23];
  real p_tor3 = system->reax_param.gp.l[24];
  real p_tor4 = system->reax_param.gp.l[25];
  real p_cot2 = system->reax_param.gp.l[27];
  reax_list *bonds = (*lists) + BONDS;
  reax_list *thb_intrs = (*lists) + THREE_BODIES;
  // char  fname[100];
  // FILE *ftor;

  // Virial tallying variables
  real delil[3], deljl[3], delkl[3];
  real eng_tmp, f_scaler, fi_tmp[3], fj_tmp[3], fk_tmp[3];

  // sprintf( fname, "tor%d.out", system->my_rank );
  // ftor = fopen( fname, "w" );

  natoms = system->n;

  for( j = 0; j < natoms; ++j ) {
    type_j = system->my_atoms[j].type;
    Delta_j = workspace->Delta_boc[j];
    start_j = Start_Index(j, bonds);
    end_j = End_Index(j, bonds);
        
    for( pk = start_j; pk < end_j; ++pk ) {
      pbond_jk = &( bonds->select.bond_list[pk] );
      k = pbond_jk->nbr;
      bo_jk = &( pbond_jk->bo_data );
      BOA_jk = bo_jk->BO - control->thb_cut;
	
      /* see if there are any 3-body interactions involving j&k
	 where j is the central atom. Otherwise there is no point in
	 trying to form a 4-body interaction out of this neighborhood */
      if( system->my_atoms[j].orig_id < system->my_atoms[k].orig_id && 
	  bo_jk->BO > control->thb_cut/*0*/ && Num_Entries(pk, thb_intrs) ) {
	start_k = Start_Index(k, bonds);
	end_k = End_Index(k, bonds);	    	       
	pj = pbond_jk->sym_index; // pj points to j on k's list
	  
	/* do the same check as above: 
	   are there any 3-body interactions involving k&j 
	   where k is the central atom */
	if( Num_Entries(pj, thb_intrs) ) {
	  type_k = system->my_atoms[k].type;
	  Delta_k = workspace->Delta_boc[k];
	  r_jk = pbond_jk->d;
	    
	  start_pk = Start_Index(pk, thb_intrs );
	  end_pk = End_Index(pk, thb_intrs );
	  start_pj = Start_Index(pj, thb_intrs );
	  end_pj = End_Index(pj, thb_intrs );		
	    
	  exp_tor2_jk = exp( -p_tor2 * BOA_jk );
	  exp_cot2_jk = exp( -p_cot2 * SQR(BOA_jk - 1.5) );
	  exp_tor3_DjDk = exp( -p_tor3 * (Delta_j + Delta_k) );
	  exp_tor4_DjDk = exp( p_tor4  * (Delta_j + Delta_k) );
	  exp_tor34_inv = 1.0 / (1.0 + exp_tor3_DjDk + exp_tor4_DjDk);
	  f11_DjDk = (2.0 + exp_tor3_DjDk) * exp_tor34_inv;
	    
	    
	  /* pick i up from j-k interaction where j is the central atom */
	  for( pi = start_pk; pi < end_pk; ++pi ) {
	    p_ijk = &( thb_intrs->select.three_body_list[pi] );
	    pij = p_ijk->pthb; // pij is pointer to i on j's bond_list
	    pbond_ij = &( bonds->select.bond_list[pij] );
	    bo_ij = &( pbond_ij->bo_data );
	      
	      
	    if( bo_ij->BO > control->thb_cut/*0*/ ) {
	      i = p_ijk->thb;
	      type_i = system->my_atoms[i].type;
	      r_ij = pbond_ij->d;
	      BOA_ij = bo_ij->BO - control->thb_cut;
		
	      theta_ijk = p_ijk->theta;
	      sin_ijk = sin( theta_ijk );
	      cos_ijk = cos( theta_ijk );
	      //tan_ijk_i = 1. / tan( theta_ijk );
	      if( sin_ijk >= 0 && sin_ijk <= MIN_SINE ) 
		tan_ijk_i = cos_ijk / MIN_SINE;
	      else if( sin_ijk <= 0 && sin_ijk >= -MIN_SINE ) 
		tan_ijk_i = cos_ijk / -MIN_SINE;
	      else tan_ijk_i = cos_ijk / sin_ijk;
		
	      exp_tor2_ij = exp( -p_tor2 * BOA_ij );
	      exp_cot2_ij = exp( -p_cot2 * SQR(BOA_ij -1.5) );
		
		
	      /* pick l up from j-k interaction where k is the central atom */
	      for( pl = start_pj; pl < end_pj; ++pl ) {
		p_jkl = &( thb_intrs->select.three_body_list[pl] );
		l = p_jkl->thb;
		plk = p_jkl->pthb; //pointer to l on k's bond_list!
		pbond_kl = &( bonds->select.bond_list[plk] );
		bo_kl = &( pbond_kl->bo_data );
		type_l = system->my_atoms[l].type;
		fbh = &(system->reax_param.fbp[type_i][type_j]
			[type_k][type_l]);
		fbp = &(system->reax_param.fbp[type_i][type_j]
			[type_k][type_l].prm[0]);

		  
		if( i != l && fbh->cnt && 
		    bo_kl->BO > control->thb_cut/*0*/ &&
		    bo_ij->BO * bo_jk->BO * bo_kl->BO > control->thb_cut/*0*/ ){
		  ++num_frb_intrs;
		  r_kl = pbond_kl->d;
		  BOA_kl = bo_kl->BO - control->thb_cut;
		    
		  theta_jkl = p_jkl->theta;
		  sin_jkl = sin( theta_jkl );
		  cos_jkl = cos( theta_jkl );
		  //tan_jkl_i = 1. / tan( theta_jkl );
		  if( sin_jkl >= 0 && sin_jkl <= MIN_SINE ) 
		    tan_jkl_i = cos_jkl / MIN_SINE;
		  else if( sin_jkl <= 0 && sin_jkl >= -MIN_SINE ) 
		    tan_jkl_i = cos_jkl / -MIN_SINE;
		  else tan_jkl_i = cos_jkl /sin_jkl;
		    
		  rvec_ScaledSum( dvec_li, 1., system->my_atoms[i].x, 
				  -1., system->my_atoms[l].x );
		  r_li = rvec_Norm( dvec_li );				 
		    
		    
		  /* omega and its derivative */
		  omega = Calculate_Omega( pbond_ij->dvec, r_ij, 
					   pbond_jk->dvec, r_jk, 
					   pbond_kl->dvec, r_kl,
					   dvec_li, r_li,
					   p_ijk, p_jkl,
					   dcos_omega_di, dcos_omega_dj,
					   dcos_omega_dk, dcos_omega_dl,
					   out_control );
		    
		  cos_omega = cos( omega );
		  cos2omega = cos( 2. * omega );
		  cos3omega = cos( 3. * omega );
		  /* end omega calculations */
		    
		  /* torsion energy */
		  exp_tor1 = exp( fbp->p_tor1 * 
				  SQR(2.0 - bo_jk->BO_pi - f11_DjDk) );
		  exp_tor2_kl = exp( -p_tor2 * BOA_kl );
		  exp_cot2_kl = exp( -p_cot2 * SQR(BOA_kl - 1.5) );
		  fn10 = (1.0 - exp_tor2_ij) * (1.0 - exp_tor2_jk) * 
		    (1.0 - exp_tor2_kl);
		    
		  CV = 0.5 * ( fbp->V1 * (1.0 + cos_omega) + 
			       fbp->V2 * exp_tor1 * (1.0 - cos2omega) +
			       fbp->V3 * (1.0 + cos3omega) );
		    
		  data->my_en.e_tor += e_tor = fn10 * sin_ijk * sin_jkl * CV;
		    
		  dfn11 = (-p_tor3 * exp_tor3_DjDk +
			   (p_tor3 * exp_tor3_DjDk - p_tor4 * exp_tor4_DjDk) *
			   (2.0 + exp_tor3_DjDk) * exp_tor34_inv) * 
		    exp_tor34_inv;
		    
		  CEtors1 = sin_ijk * sin_jkl * CV;
		    
		  CEtors2 = -fn10 * 2.0 * fbp->p_tor1 * fbp->V2 * exp_tor1 *
		    (2.0 - bo_jk->BO_pi - f11_DjDk) * (1.0 - SQR(cos_omega)) *
		    sin_ijk * sin_jkl; 
		  CEtors3 = CEtors2 * dfn11;
		    
		  CEtors4 = CEtors1 * p_tor2 * exp_tor2_ij * 
		    (1.0 - exp_tor2_jk) * (1.0 - exp_tor2_kl);
		  CEtors5 = CEtors1 * p_tor2 * 
		    (1.0 - exp_tor2_ij) * exp_tor2_jk * (1.0 - exp_tor2_kl);
		  CEtors6 = CEtors1 * p_tor2 * 
		    (1.0 - exp_tor2_ij) * (1.0 - exp_tor2_jk) * exp_tor2_kl;
		    
		  cmn = -fn10 * CV;
		  CEtors7 = cmn * sin_jkl * tan_ijk_i;
		  CEtors8 = cmn * sin_ijk * tan_jkl_i;
		    
		  CEtors9 = fn10 * sin_ijk * sin_jkl * 
		    (0.5 * fbp->V1 - 2.0 * fbp->V2 * exp_tor1 * cos_omega +
		     1.5 * fbp->V3 * (cos2omega + 2.0 * SQR(cos_omega)));
		  /* end  of torsion energy */
				  

		  /* 4-body conjugation energy */
		  fn12 = exp_cot2_ij * exp_cot2_jk * exp_cot2_kl;
		  data->my_en.e_con += e_con =
		    fbp->p_cot1 * fn12 * 
		    (1.0 + (SQR(cos_omega) - 1.0) * sin_ijk * sin_jkl);
		    
		  Cconj = -2.0 * fn12 * fbp->p_cot1 * p_cot2 * 
		    (1.0 + (SQR(cos_omega) - 1.0) * sin_ijk * sin_jkl);
		    
		  CEconj1 = Cconj * (BOA_ij - 1.5e0);
		  CEconj2 = Cconj * (BOA_jk - 1.5e0);
		  CEconj3 = Cconj * (BOA_kl - 1.5e0);
		    
		  CEconj4 = -fbp->p_cot1 * fn12 * 
		    (SQR(cos_omega) - 1.0) * sin_jkl * tan_ijk_i;
		  CEconj5 = -fbp->p_cot1 * fn12 * 
		    (SQR(cos_omega) - 1.0) * sin_ijk * tan_jkl_i;
		  CEconj6 = 2.0 * fbp->p_cot1 * fn12 * 
		    cos_omega * sin_ijk * sin_jkl;
		  /* end 4-body conjugation energy */
		  
		  /* forces */
		  bo_jk->Cdbopi += CEtors2;
		  workspace->CdDelta[j] += CEtors3;
		  workspace->CdDelta[k] += CEtors3;
		  bo_ij->Cdbo += (CEtors4 + CEconj1);
		  bo_jk->Cdbo += (CEtors5 + CEconj2);
		  bo_kl->Cdbo += (CEtors6 + CEconj3);
		    
		  if( control->virial == 0 ) {
		    /* dcos_theta_ijk */
		    rvec_ScaledAdd( workspace->f[i], 
				    CEtors7 + CEconj4, p_ijk->dcos_dk );
		    rvec_ScaledAdd( workspace->f[j], 
				    CEtors7 + CEconj4, p_ijk->dcos_dj );
		    rvec_ScaledAdd( workspace->f[k], 
				    CEtors7 + CEconj4, p_ijk->dcos_di );
		      
		    /* dcos_theta_jkl */
		    rvec_ScaledAdd( workspace->f[j], 
				    CEtors8 + CEconj5, p_jkl->dcos_di );
		    rvec_ScaledAdd( workspace->f[k], 
				    CEtors8 + CEconj5, p_jkl->dcos_dj );
		    rvec_ScaledAdd( workspace->f[l], 
				    CEtors8 + CEconj5, p_jkl->dcos_dk );
		      
		    /* dcos_omega */
		    rvec_ScaledAdd( workspace->f[i], 
				    CEtors9 + CEconj6, dcos_omega_di );
		    rvec_ScaledAdd( workspace->f[j], 
				    CEtors9 + CEconj6, dcos_omega_dj );
		    rvec_ScaledAdd( workspace->f[k], 
				    CEtors9 + CEconj6, dcos_omega_dk );
		    rvec_ScaledAdd( workspace->f[l], 
				    CEtors9 + CEconj6, dcos_omega_dl );
		  }
		  else {
		    ivec_Sum(rel_box_jl, pbond_jk->rel_box, pbond_kl->rel_box);
		      
		    /* dcos_theta_ijk */
		    rvec_Scale( force, CEtors7 + CEconj4, p_ijk->dcos_dk );
		    rvec_Add( workspace->f[i], force );
		    rvec_iMultiply( ext_press, pbond_ij->rel_box, force );
		    rvec_Add( data->my_ext_press, ext_press );
			
		    rvec_ScaledAdd( workspace->f[j], 
				    CEtors7 + CEconj4, p_ijk->dcos_dj );
			
		    rvec_Scale( force, CEtors7 + CEconj4, p_ijk->dcos_di );
		    rvec_Add( workspace->f[k], force );
		    rvec_iMultiply( ext_press, pbond_jk->rel_box, force );
		    rvec_Add( data->my_ext_press, ext_press );
		      
		      
		    /* dcos_theta_jkl */
		    rvec_ScaledAdd( workspace->f[j], 
				    CEtors8 + CEconj5, p_jkl->dcos_di );
			
		    rvec_Scale( force, CEtors8 + CEconj5, p_jkl->dcos_dj );
		    rvec_Add( workspace->f[k], force );
		    rvec_iMultiply( ext_press, pbond_jk->rel_box, force );
		    rvec_Add( data->my_ext_press, ext_press );
		      
		    rvec_Scale( force, CEtors8 + CEconj5, p_jkl->dcos_dk );
		    rvec_Add( workspace->f[l], force );
		    rvec_iMultiply( ext_press, rel_box_jl, force );
		    rvec_Add( data->my_ext_press, ext_press );
		      
		      
		    /* dcos_omega */				      
		    rvec_Scale( force, CEtors9 + CEconj6, dcos_omega_di );
		    rvec_Add( workspace->f[i], force );
		    rvec_iMultiply( ext_press, pbond_ij->rel_box, force );
		    rvec_Add( data->my_ext_press, ext_press );
		      
		    rvec_ScaledAdd( workspace->f[j], 
				    CEtors9 + CEconj6, dcos_omega_dj );
		      
		    rvec_Scale( force, CEtors9 + CEconj6, dcos_omega_dk );
		    rvec_Add( workspace->f[k], force );
		    rvec_iMultiply( ext_press, pbond_jk->rel_box, force );
		    rvec_Add( data->my_ext_press, ext_press );
		      
		    rvec_Scale( force, CEtors9 + CEconj6, dcos_omega_dl );
		    rvec_Add( workspace->f[l], force );
		    rvec_iMultiply( ext_press, rel_box_jl, force );
		    rvec_Add( data->my_ext_press, ext_press );
		  }

		  /* tally into per-atom virials */
		  if( system->pair_ptr->vflag_atom || system->pair_ptr->evflag) {

		    // acquire vectors
		    rvec_ScaledSum( delil, 1., system->my_atoms[l].x,
					  -1., system->my_atoms[i].x );
		    rvec_ScaledSum( deljl, 1., system->my_atoms[l].x,
					  -1., system->my_atoms[j].x );
		    rvec_ScaledSum( delkl, 1., system->my_atoms[l].x,
					  -1., system->my_atoms[k].x );
		    // dcos_theta_ijk
		    rvec_Scale( fi_tmp, CEtors7 + CEconj4, p_ijk->dcos_dk );
		    rvec_Scale( fj_tmp, CEtors7 + CEconj4, p_ijk->dcos_dj );
		    rvec_Scale( fk_tmp, CEtors7 + CEconj4, p_ijk->dcos_di );

		    // dcos_theta_jkl
		    rvec_ScaledAdd( fj_tmp, CEtors8 + CEconj5, p_jkl->dcos_di );
		    rvec_ScaledAdd( fk_tmp, CEtors8 + CEconj5, p_jkl->dcos_dj );

		    // dcos_omega
		    rvec_ScaledAdd( fi_tmp, CEtors9 + CEconj6, dcos_omega_di );
		    rvec_ScaledAdd( fj_tmp, CEtors9 + CEconj6, dcos_omega_dj );
		    rvec_ScaledAdd( fk_tmp, CEtors9 + CEconj6, dcos_omega_dk );

		    // tally
		    eng_tmp = e_tor + e_con;
		    if( system->pair_ptr->evflag)
			    system->pair_ptr->ev_tally(j,k,natoms,1,eng_tmp,0.0,0.0,0.0,0.0,0.0);
		    if( system->pair_ptr->vflag_atom)
			    system->pair_ptr->v_tally4(i,j,k,l,fi_tmp,fj_tmp,fk_tmp,delil,deljl,delkl);
		  }
		  

#ifdef TEST_ENERGY
		  /* fprintf( out_control->etor, 
		     "%12.8f%12.8f%12.8f%12.8f%12.8f%12.8f%12.8f\n",
		     r_ij, r_jk, r_kl, cos_ijk, cos_jkl, sin_ijk, sin_jkl );
		     fprintf( out_control->etor, "%12.8f\n", dfn11 ); */
		  /* fprintf( out_control->etor, 
		     "%12.8f%12.8f%12.8f%12.8f%12.8f%12.8f%12.8f%12.8f\n",
		     CEtors2, CEtors3, CEtors4, CEtors5, CEtors6, 
		     CEtors7, CEtors8, CEtors9 ); */
		  /* fprintf( out_control->etor, 
		     "%12.8f%12.8f%12.8f%12.8f%12.8f%12.8f%12.8f%12.8f\n",
		     htra, htrb, htrc, hthd, hthe, hnra, hnrc, hnhd, hnhe ); */
		  /* fprintf( out_control->etor, 
		     "%12.8f%12.8f%12.8f%12.8f%12.8f%12.8f\n",
		     CEconj1, CEconj2, CEconj3, CEconj4, CEconj5, CEconj6 ); */
				  
		  /* fprintf( out_control->etor, "%12.6f%12.6f%12.6f%12.6f\n",
		     fbp->V1, fbp->V2, fbp->V3, fbp->p_tor1 );*/
		  
		  fprintf(out_control->etor, 
			  //"%6d%6d%6d%6d%24.15e%24.15e%24.15e%24.15e\n", 
			  "%6d%6d%6d%6d%12.4f%12.4f%12.4f%12.4f\n", 
			  system->my_atoms[i].orig_id,system->my_atoms[j].orig_id, 
			  system->my_atoms[k].orig_id,system->my_atoms[l].orig_id, 
			   RAD2DEG(omega), BOA_jk, e_tor, data->my_en.e_tor );
		  
		  fprintf(out_control->econ, 
			  //"%6d%6d%6d%6d%24.15e%24.15e%24.15e%24.15e%24.15e%24.15e\n", 
			  "%6d%6d%6d%6d%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f\n", 
			  system->my_atoms[i].orig_id,system->my_atoms[j].orig_id, 
			  system->my_atoms[k].orig_id,system->my_atoms[l].orig_id, 
			  RAD2DEG(omega), BOA_ij, BOA_jk, BOA_kl, 
			  e_con, data->my_en.e_con );
#endif
		    
#ifdef TEST_FORCES
		  /* Torsion Forces */
		  Add_dBOpinpi2( system, lists, j, pk, CEtors2, 0.0, 
				 workspace->f_tor, workspace->f_tor );
		  Add_dDelta( system, lists, j, CEtors3, workspace->f_tor );
		  Add_dDelta( system, lists, k, CEtors3, workspace->f_tor );
		  Add_dBO( system, lists, j, pij, CEtors4, workspace->f_tor );
		  Add_dBO( system, lists, j, pk, CEtors5, workspace->f_tor );
		  Add_dBO( system, lists, k, plk, CEtors6, workspace->f_tor );
		    
		  rvec_ScaledAdd( workspace->f_tor[i], 
				  CEtors7, p_ijk->dcos_dk );
		  rvec_ScaledAdd( workspace->f_tor[j], 
				  CEtors7, p_ijk->dcos_dj );
		  rvec_ScaledAdd( workspace->f_tor[k], 
				  CEtors7, p_ijk->dcos_di );
			  
		  rvec_ScaledAdd( workspace->f_tor[j], 
				  CEtors8, p_jkl->dcos_di );
		  rvec_ScaledAdd( workspace->f_tor[k], 
				  CEtors8, p_jkl->dcos_dj );
		  rvec_ScaledAdd( workspace->f_tor[l], 
				  CEtors8, p_jkl->dcos_dk );
			  
		  rvec_ScaledAdd( workspace->f_tor[i], 
				  CEtors9, dcos_omega_di );
		  rvec_ScaledAdd( workspace->f_tor[j], 
				  CEtors9, dcos_omega_dj );
		  rvec_ScaledAdd( workspace->f_tor[k], 
				  CEtors9, dcos_omega_dk );
		  rvec_ScaledAdd( workspace->f_tor[l], 
				  CEtors9, dcos_omega_dl );
		    
		  /* Conjugation Forces */
		  Add_dBO( system, lists, j, pij, CEconj1, workspace->f_con );
		  Add_dBO( system, lists, j, pk, CEconj2, workspace->f_con );
		  Add_dBO( system, lists, k, plk, CEconj3, workspace->f_con );
		    
		  rvec_ScaledAdd( workspace->f_con[i], 
				  CEconj4, p_ijk->dcos_dk );
		  rvec_ScaledAdd( workspace->f_con[j], 
				  CEconj4, p_ijk->dcos_dj );
		  rvec_ScaledAdd( workspace->f_con[k], 
				  CEconj4, p_ijk->dcos_di );
		    
		  rvec_ScaledAdd( workspace->f_con[j], 
				  CEconj5, p_jkl->dcos_di );
		  rvec_ScaledAdd( workspace->f_con[k], 
				  CEconj5, p_jkl->dcos_dj );
		  rvec_ScaledAdd( workspace->f_con[l], 
				  CEconj5, p_jkl->dcos_dk );
		    
		  rvec_ScaledAdd( workspace->f_con[i], 
				  CEconj6, dcos_omega_di );
		  rvec_ScaledAdd( workspace->f_con[j], 
				  CEconj6, dcos_omega_dj );
		  rvec_ScaledAdd( workspace->f_con[k], 
				  CEconj6, dcos_omega_dk );
		  rvec_ScaledAdd( workspace->f_con[l], 
				  CEconj6, dcos_omega_dl );
#endif
		} // pl check ends
	      } // pl loop ends
	    } // pi check ends
	  } // pi loop ends
	} // k-j neighbor check ends
      } // j<k && j-k neighbor check ends
    } // pk loop ends
  } // j loop
  
#if defined(DEBUG)
  fprintf( stderr, "Number of torsion angles: %d\n", num_frb_intrs );
  fprintf( stderr, "Torsion Energy: %g\t Conjugation Energy: %g\n", 
	   data->my_en.e_tor, data->my_en.e_con );

  fprintf( stderr, "4body: ext_press (%12.6f %12.6f %12.6f)\n",
	   data->ext_press[0], data->ext_press[1], data->ext_press[2] );
#endif
}
