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

#include "reaxc_hydrogen_bonds_omp.h"
#include "pair_reaxc_omp.h"
#include "reaxc_defs.h"
#include "reaxc_bond_orders_omp.h"
#include "reaxc_list.h"
#include "reaxc_valence_angles.h"     // To access Calculate_Theta()
#include "reaxc_valence_angles_omp.h" // To access Calculate_dCos_ThetaOMP()
#include "reaxc_vector.h"

#if defined(_OPENMP)
#include  <omp.h>
#endif

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

void Hydrogen_BondsOMP( reax_system *system, control_params *control,
                     simulation_data *data, storage *workspace,
                     reax_list **lists, output_controls * /* out_control */)
{
#ifdef OMP_TIMING
  double endTimeBase, startTimeBase;
  startTimeBase = MPI_Wtime();
#endif

  const int nthreads = control->nthreads;

#if defined(_OPENMP)
#pragma omp parallel default(shared) //default(none)
#endif
  {
  int  i, j, k, pi, pk;
  int  type_i, type_j, type_k;
  int  start_j, end_j, hb_start_j, hb_end_j;
  int  hblist[MAX_BONDS];
  int  itr, top;
  int  num_hb_intrs = 0;
  ivec rel_jk;
  double r_jk, theta, cos_theta, sin_xhz4, cos_xhz1, sin_theta2;
  double e_hb, e_hb_thr = 0.0, exp_hb2, exp_hb3, CEhb1, CEhb2, CEhb3;
  rvec dcos_theta_di, dcos_theta_dj, dcos_theta_dk;
  rvec dvec_jk, force, ext_press;
  hbond_parameters *hbp;
  bond_order_data *bo_ij;
  bond_data *pbond_ij;
  far_neighbor_data *nbr_jk;
  reax_list *bonds, *hbonds;
  bond_data *bond_list;
  hbond_data *hbond_list;

  // tally variables
  double fi_tmp[3], fk_tmp[3], delij[3], delkj[3];

  bonds = (*lists) + BONDS;
  bond_list = bonds->select.bond_list;
  hbonds = (*lists) + HBONDS;
  hbond_list = hbonds->select.hbond_list;

  int natoms = system->n;
#if defined(_OPENMP)
  int tid = omp_get_thread_num();
#else
  int tid = 0;
#endif
  const int idelta = 1 + natoms/nthreads;
  int ifrom = tid*idelta;
  int ito   = ((ifrom + idelta) > natoms) ? natoms : ifrom + idelta;

  long reductionOffset = (system->N * tid);

  class PairReaxCOMP *pair_reax_ptr;
  pair_reax_ptr = static_cast<class PairReaxCOMP*>(system->pair_ptr);

  class ThrData *thr = pair_reax_ptr->getFixOMP()->get_thr(tid);

  /* loops below discover the Hydrogen bonds between i-j-k triplets.
     here j is H atom and there has to be some bond between i and j.
     Hydrogen bond is between j and k.
     so in this function i->X, j->H, k->Z when we map
     variables onto the ones in the handout.*/
  //  for( j = 0; j < system->n; ++j )
  for( j = ifrom; j < ito; ++j ) {
    /* j has to be of type H */
    if (system->reax_param.sbp[system->my_atoms[j].type].p_hbond == 1) {
      /*set j's variables */
      type_j     = system->my_atoms[j].type;
      start_j    = Start_Index(j, bonds);
      end_j      = End_Index(j, bonds);
      hb_start_j = Start_Index( system->my_atoms[j].Hindex, hbonds );
      hb_end_j   = End_Index( system->my_atoms[j].Hindex, hbonds );
      if(type_j < 0) continue;

      top = 0;
      for( pi = start_j; pi < end_j; ++pi )  {
        pbond_ij = &( bond_list[pi] );
        i = pbond_ij->nbr;
        type_i = system->my_atoms[i].type;
        if(type_i < 0) continue;
        bo_ij = &(pbond_ij->bo_data);

        if( system->reax_param.sbp[type_i].p_hbond == 2 &&
            bo_ij->BO >= HB_THRESHOLD )
          hblist[top++] = pi;
      }

      for( pk = hb_start_j; pk < hb_end_j; ++pk ) {
        /* set k's varibles */
        k = hbond_list[pk].nbr;
        type_k = system->my_atoms[k].type;
        if(type_k < 0) continue;
        nbr_jk = hbond_list[pk].ptr;
        r_jk = nbr_jk->d;
        rvec_Scale( dvec_jk, hbond_list[pk].scl, nbr_jk->dvec );

        for( itr = 0; itr < top; ++itr ) {
          pi = hblist[itr];
          pbond_ij = &( bonds->select.bond_list[pi] );
          i = pbond_ij->nbr;

          if (system->my_atoms[i].orig_id != system->my_atoms[k].orig_id) {
            bo_ij = &(pbond_ij->bo_data);
            type_i = system->my_atoms[i].type;
            if(type_i < 0) continue;
            hbp = &(system->reax_param.hbp[ type_i ][ type_j ][ type_k ]);
            ++num_hb_intrs;

            Calculate_Theta( pbond_ij->dvec, pbond_ij->d, dvec_jk, r_jk,
                             &theta, &cos_theta );
            /* the derivative of cos(theta) */
            Calculate_dCos_ThetaOMP( pbond_ij->dvec, pbond_ij->d, dvec_jk, r_jk,
                                     &dcos_theta_di, &dcos_theta_dj,
                                     &dcos_theta_dk );

            /* hydrogen bond energy*/
            sin_theta2 = sin( theta/2.0 );
            sin_xhz4 = SQR(sin_theta2);
            sin_xhz4 *= sin_xhz4;
            cos_xhz1 = ( 1.0 - cos_theta );
            exp_hb2 = exp( -hbp->p_hb2 * bo_ij->BO );
            exp_hb3 = exp( -hbp->p_hb3 * ( hbp->r0_hb / r_jk +
                                           r_jk / hbp->r0_hb - 2.0 ) );

            e_hb_thr += e_hb = hbp->p_hb1 * (1.0 - exp_hb2) * exp_hb3 * sin_xhz4;

            CEhb1 = hbp->p_hb1 * hbp->p_hb2 * exp_hb2 * exp_hb3 * sin_xhz4;
            CEhb2 = -hbp->p_hb1/2.0 * (1.0 - exp_hb2) * exp_hb3 * cos_xhz1;
            CEhb3 = -hbp->p_hb3 *
              (-hbp->r0_hb / SQR(r_jk) + 1.0 / hbp->r0_hb) * e_hb;

            /* hydrogen bond forces */
            bo_ij->Cdbo += CEhb1; // dbo term

            if (control->virial == 0) {
              // dcos terms
              rvec_ScaledAdd(workspace->forceReduction[reductionOffset+i], +CEhb2, dcos_theta_di );
              rvec_ScaledAdd(workspace->forceReduction[reductionOffset+j], +CEhb2, dcos_theta_dj );
              rvec_ScaledAdd(workspace->forceReduction[reductionOffset+k], +CEhb2, dcos_theta_dk );
              // dr terms
              rvec_ScaledAdd(workspace->forceReduction[reductionOffset+j], -CEhb3/r_jk, dvec_jk );
              rvec_ScaledAdd(workspace->forceReduction[reductionOffset+k], +CEhb3/r_jk, dvec_jk );
            }
            else {
              /* for pressure coupling, terms that are not related to bond order
                 derivatives are added directly into pressure vector/tensor */
              rvec_Scale( force, +CEhb2, dcos_theta_di ); // dcos terms
              rvec_Add(workspace->forceReduction[reductionOffset+i], force );
              rvec_iMultiply( ext_press, pbond_ij->rel_box, force );
              rvec_ScaledAdd( workspace->my_ext_pressReduction[tid],1.0, ext_press );

              rvec_ScaledAdd(workspace->forceReduction[reductionOffset+j], +CEhb2, dcos_theta_dj );

              ivec_Scale( rel_jk, hbond_list[pk].scl, nbr_jk->rel_box );
              rvec_Scale( force, +CEhb2, dcos_theta_dk );
              rvec_Add(workspace->forceReduction[reductionOffset+k], force );
              rvec_iMultiply( ext_press, rel_jk, force );
              rvec_ScaledAdd( workspace->my_ext_pressReduction[tid],1.0, ext_press );
              // dr terms
              rvec_ScaledAdd(workspace->forceReduction[reductionOffset+j],-CEhb3/r_jk, dvec_jk );

              rvec_Scale( force, CEhb3/r_jk, dvec_jk );
              rvec_Add(workspace->forceReduction[reductionOffset+k], force );
              rvec_iMultiply( ext_press, rel_jk, force );
              rvec_ScaledAdd( workspace->my_ext_pressReduction[tid],1.0, ext_press );
            }

            /* tally into per-atom virials */
            if (system->pair_ptr->vflag_atom || system->pair_ptr->evflag) {
              rvec_ScaledSum( delij, 1., system->my_atoms[j].x,
                                    -1., system->my_atoms[i].x );
              rvec_ScaledSum( delkj, 1., system->my_atoms[j].x,
                                     -1., system->my_atoms[k].x );

              rvec_Scale(fi_tmp, CEhb2, dcos_theta_di);
              rvec_Scale(fk_tmp, CEhb2, dcos_theta_dk);
              rvec_ScaledAdd(fk_tmp, CEhb3/r_jk, dvec_jk);

              pair_reax_ptr->ev_tally3_thr_proxy(system->pair_ptr,i,j,k,e_hb,0.0,fi_tmp,fk_tmp,delij,delkj,thr);
            }
          }
        }
      }

    }
  }
#if defined(_OPENMP)
#pragma omp critical
#endif
  {
    data->my_en.e_hb += e_hb_thr;
  }
}

#ifdef OMP_TIMING
  endTimeBase = MPI_Wtime();
  ompTimingData[COMPUTEHBONDSINDEX] += (endTimeBase-startTimeBase);
#endif
}
