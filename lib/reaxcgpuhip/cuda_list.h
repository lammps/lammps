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

#ifndef __CUDA_LIST_H_
#define __CUDA_LIST_H_

#include "reaxc_types.h"


#ifdef __cplusplus
extern "C" {
#endif

void Cuda_Make_List( int, int, int, reax_list* );

void Cuda_Delete_List( reax_list* );

#ifdef __cplusplus
}
#endif


static inline CUDA_HOST_DEVICE int Cuda_Num_Entries( int i, reax_list *l )
{
    return l->end_index[i] - l->index[i];
}

static inline CUDA_HOST_DEVICE int Cuda_Start_Index( int i, reax_list *l )
{
    return l->index[i];
}

static inline CUDA_HOST_DEVICE int Cuda_End_Index( int i, reax_list *l )
{
    return l->end_index[i];
}

static inline CUDA_HOST_DEVICE void Cuda_Set_Start_Index( int i, int val, reax_list *l )
{
    l->index[i] = val;
}

static inline CUDA_HOST_DEVICE void Cuda_Set_End_Index( int i, int val, reax_list *l )
{
    l->end_index[i] = val;
}


#endif
