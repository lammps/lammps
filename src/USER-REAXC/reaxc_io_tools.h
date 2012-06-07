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

#ifndef __IO_TOOLS_H_
#define __IO_TOOLS_H_

#include "reaxc_types.h"

int Init_Output_Files( reax_system*, control_params*,
                       output_controls*, mpi_datatypes*, char* );
int Close_Output_Files( reax_system*, control_params*,
                        output_controls*, mpi_datatypes* );

void  Print_Box( simulation_box*, char*, FILE* );

void  Print_Grid( grid*, FILE* );
void  Print_GCell_Exchange_Bounds( int, neighbor_proc* );
void  Print_Native_GCells( reax_system* );
void  Print_All_GCells( reax_system*);

void  Print_Init_Atoms( reax_system*, storage* );
void  Print_My_Atoms( reax_system* );
void  Print_My_Ext_Atoms( reax_system* );

void  Print_Far_Neighbors( reax_system*, reax_list**, control_params *);
void  Print_Sparse_Matrix( reax_system*, sparse_matrix* );
void  Print_Sparse_Matrix2( reax_system*, sparse_matrix*, char* );
void  Print_Linear_System( reax_system*, control_params*, storage*, int );
void  Print_LinSys_Soln( reax_system*, real*, real*, real* );
void  Print_Charges( reax_system* );
void  Print_Bonds( reax_system*, reax_list*, char* );
void  Print_Bond_List2( reax_system*, reax_list*, char* );
void  Print_Total_Force( reax_system*, simulation_data*, storage* );
void  Output_Results( reax_system*, control_params*, simulation_data*,
                      reax_list**, output_controls*, mpi_datatypes* );
void  fixbond( reax_system*, control_params*, simulation_data*,
                      reax_list**, output_controls*, mpi_datatypes* );
void  fixspecies( reax_system*, control_params*, simulation_data*,
                      reax_list**, output_controls*, mpi_datatypes* );

#if defined(DEBUG_FOCUS) || defined(TEST_FORCES) || defined(TEST_ENERGY)
void Debug_Marker_Bonded( output_controls*, int );
void Debug_Marker_Nonbonded( output_controls*, int );
void  Print_Near_Neighbors_List( reax_system*, reax_list**, control_params*,
                                 simulation_data*, output_controls* );
void  Print_Far_Neighbors_List( reax_system*, reax_list**, control_params*,
                                simulation_data*, output_controls* );
void  Print_Bond_List( reax_system*, control_params*, simulation_data*,
                       reax_list**, output_controls* );
/*void Dummy_Printer( reax_system*, control_params*, simulation_data*,
                    storage*, reax_list**, output_controls* );
void Print_Bond_Orders( reax_system*, control_params*, simulation_data*,
                        storage*, reax_list**, output_controls* );
void Print_Bond_Forces( reax_system*, control_params*, simulation_data*,
                        storage*, reax_list**, output_controls* );
void Print_LonePair_Forces( reax_system*, control_params*, simulation_data*,
                            storage*, reax_list**, output_controls* );
void Print_OverUnderCoor_Forces( reax_system*, control_params*,
                                 simulation_data*, storage*, reax_list**,
                                 output_controls* );
void Print_Three_Body_Forces( reax_system*, control_params*, simulation_data*,
                                 storage*, reax_list**, output_controls* );
void Print_Hydrogen_Bond_Forces( reax_system*, control_params*,
                                 simulation_data*, storage*, reax_list**,
                                 output_controls* );
void Print_Four_Body_Forces( reax_system*, control_params*, simulation_data*,
                                 storage*, reax_list**, output_controls* );
void Print_vdW_Coulomb_Forces( reax_system*, control_params*,
                               simulation_data*, storage*, reax_list**,
                               output_controls* );
void Print_Total_Force( reax_system*, control_params*, simulation_data*,
                        storage*, reax_list**, output_controls* );
void Compare_Total_Forces( reax_system*, control_params*, simulation_data*,
storage*, reax_list**, output_controls* );*/
//void  Print_Total_Force( reax_system*, control_params* );
void Print_Force_Files( reax_system*, control_params*, simulation_data*,
                        storage*, reax_list**, output_controls*,
                        mpi_datatypes * );
//void Init_Force_Test_Functions( );

int fn_qsort_intcmp( const void *, const void * );

void Print_Far_Neighbors_List( reax_system*, reax_list**, control_params*,
                               simulation_data*, output_controls* );

void Print_Near_Neighbors_List( reax_system*, reax_list**, control_params*,
                                simulation_data*, output_controls* );

void Print_Bond_List( reax_system*, control_params*, simulation_data*,
                      reax_list**, output_controls*);

#endif
#endif
