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

#ifndef __VECTOR_H_
#define __VECTOR_H_

#include "pair.h"
#include "reaxc_types.h"
#include "reaxc_defs.h"

void rvec_Copy( rvec, rvec );
void rvec_Scale( rvec, real, rvec );
void rvec_Add( rvec, rvec );
void rvec_ScaledAdd( rvec, real, rvec );
void rvec_Sum( rvec, rvec, rvec );
void rvec_ScaledSum( rvec, real, rvec, real, rvec );
real rvec_Dot( rvec, rvec );
real rvec_ScaledDot( real, rvec, real, rvec );
void rvec_Multiply( rvec, rvec, rvec );
void rvec_iMultiply( rvec, ivec, rvec );
void rvec_Divide( rvec, rvec, rvec );
void rvec_iDivide( rvec, rvec, ivec );
void rvec_Invert( rvec, rvec );
void rvec_Cross( rvec, rvec, rvec );
void rvec_OuterProduct( rtensor, rvec, rvec );
real rvec_Norm_Sqr( rvec );
real rvec_Norm( rvec );
int  rvec_isZero( rvec );
void rvec_MakeZero( rvec );
void rvec_Random( rvec );

void rtensor_MakeZero( rtensor );
void rtensor_MatVec( rvec, rtensor, rvec );
void rtensor_Scale( rtensor, real, rtensor );

void ivec_MakeZero( ivec );
void ivec_Copy( ivec, ivec );
void ivec_Scale( ivec, real, ivec );
void ivec_Sum( ivec, ivec, ivec );

#endif
