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
#include "list.h"
#include "tool_box.h"
#elif defined(LAMMPS_REAX)
#include "reaxc_list.h"
#include "reaxc_tool_box.h"
#endif


/************* allocate list space ******************/
int Make_List(int n, int num_intrs, int type, reax_list *l, MPI_Comm comm)
{
  l->allocated = 1;
  
  l->n = n;
  l->num_intrs = num_intrs;
    
  l->index = (int*) smalloc( n * sizeof(int), "list:index", comm );
  l->end_index = (int*) smalloc( n * sizeof(int), "list:end_index", comm );

  l->type = type;
#if defined(DEBUG_FOCUS)
  fprintf( stderr, "list: n=%d num_intrs=%d type=%d\n", n, num_intrs, type );
#endif

  switch(l->type) {
  case TYP_VOID:
    l->select.v = (void*) smalloc(l->num_intrs * sizeof(void*), "list:v", comm);
    break;
    
  case TYP_THREE_BODY:
    l->select.three_body_list = (three_body_interaction_data*) 
      smalloc( l->num_intrs * sizeof(three_body_interaction_data), 
	       "list:three_bodies", comm );
    break;
    
  case TYP_BOND:
    l->select.bond_list = (bond_data*) 
      smalloc( l->num_intrs * sizeof(bond_data), "list:bonds", comm );
    break;
    
  case TYP_DBO:
    l->select.dbo_list = (dbond_data*) 
      smalloc( l->num_intrs * sizeof(dbond_data), "list:dbonds", comm );
    break;
    
  case TYP_DDELTA:
    l->select.dDelta_list = (dDelta_data*) 
      smalloc( l->num_intrs * sizeof(dDelta_data), "list:dDeltas", comm );
    break;
    
  case TYP_FAR_NEIGHBOR:
    l->select.far_nbr_list = (far_neighbor_data*) 
      smalloc(l->num_intrs * sizeof(far_neighbor_data), "list:far_nbrs", comm);
    break;
        
  case TYP_HBOND:
    l->select.hbond_list = (hbond_data*)
      smalloc( l->num_intrs * sizeof(hbond_data), "list:hbonds", comm );
    break;			
    
  default:
    fprintf( stderr, "ERROR: no %d list type defined!\n", l->type );
    MPI_Abort( comm, INVALID_INPUT );
  }

  return SUCCESS;
}


void Delete_List( reax_list *l, MPI_Comm comm )
{
  if( l->allocated == 0 )
    return;
  l->allocated = 0;

  sfree( l->index, "list:index" );
  sfree( l->end_index, "list:end_index" );
  
  switch(l->type) {
  case TYP_VOID:
    sfree( l->select.v, "list:v" );
    break;
  case TYP_HBOND:
    sfree( l->select.hbond_list, "list:hbonds" );
    break;
  case TYP_FAR_NEIGHBOR:
    sfree( l->select.far_nbr_list, "list:far_nbrs" );
    break;
  case TYP_BOND:
    sfree( l->select.bond_list, "list:bonds" );
    break;
  case TYP_DBO:
    sfree( l->select.dbo_list, "list:dbos" );
    break;
  case TYP_DDELTA:
    sfree( l->select.dDelta_list, "list:dDeltas" );
    break;
  case TYP_THREE_BODY:
    sfree( l->select.three_body_list, "list:three_bodies" );
    break;

  default:
    fprintf( stderr, "ERROR: no %d list type defined!\n", l->type );
    MPI_Abort( comm, INVALID_INPUT );
  }
}


#if defined(PURE_REAX)
inline int Num_Entries( int i, reax_list *l )
{
  return l->end_index[i] - l->index[i];
}

inline int Start_Index( int i, reax_list *l )
{
  return l->index[i];
}

inline int End_Index( int i, reax_list *l )
{
  return l->end_index[i];
}

inline void Set_Start_Index( int i, int val, reax_list *l )
{
  l->index[i] = val;
}

inline void Set_End_Index( int i, int val, reax_list *l )
{
  l->end_index[i] = val;
}
#endif
