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

#ifndef __CUDA_VALENCE_ANGLES_H_
#define __CUDA_VALENCE_ANGLES_H_

#include "reaxc_types.h"

#include "vector.h"


CUDA_GLOBAL void Cuda_Valence_Angles( reax_atom *, global_parameters,
        single_body_parameters *, three_body_header *, control_params *,
        storage, reax_list, reax_list, int, int, int, real *,
        real *, real *, rvec *);

CUDA_GLOBAL void Cuda_Valence_Angles_PostProcess ( reax_atom *, control_params *,
        storage , reax_list, int );

CUDA_GLOBAL void Estimate_Cuda_Valence_Angles( reax_atom *, control_params *,
        reax_list , int , int, int *);


/* calculates the angle (theta) between atom triplet i-j-k */
CUDA_DEVICE static inline void Calculate_Theta( rvec dvec_ji, real d_ji, rvec dvec_jk,
        real d_jk, real *theta, real *cos_theta )
{
    (*cos_theta) = Dot( dvec_ji, dvec_jk, 3 ) / ( d_ji * d_jk );

    if ( *cos_theta > 1. )
    {
        *cos_theta  = 1.0;
    }
    if ( *cos_theta < -1. )
    {
        *cos_theta  = -1.0;
    }

    (*theta) = ACOS( *cos_theta );
}


/* calculates the derivative of the cosine of the angle between atom triplet i-j-k */
CUDA_DEVICE static inline void Calculate_dCos_Theta( rvec dvec_ji, real d_ji,
        rvec dvec_jk, real d_jk, rvec* dcos_theta_di, rvec* dcos_theta_dj,
        rvec* dcos_theta_dk )
{
    int t;
    real sqr_d_ji = SQR( d_ji );
    real sqr_d_jk = SQR( d_jk );
    real inv_dists = 1.0 / (d_ji * d_jk);
    real inv_dists3 = POW( inv_dists, 3 );
    real dot_dvecs = Dot( dvec_ji, dvec_jk, 3 );
    real Cdot_inv3 = dot_dvecs * inv_dists3;

    for ( t = 0; t < 3; ++t )
    {
        (*dcos_theta_di)[t] = dvec_jk[t] * inv_dists -
            Cdot_inv3 * sqr_d_jk * dvec_ji[t];
        (*dcos_theta_dj)[t] = -(dvec_jk[t] + dvec_ji[t]) * inv_dists +
            Cdot_inv3 * ( sqr_d_jk * dvec_ji[t] + sqr_d_ji * dvec_jk[t] );
        (*dcos_theta_dk)[t] = dvec_ji[t] * inv_dists -
            Cdot_inv3 * sqr_d_ji * dvec_jk[t];
    }
}


#endif
