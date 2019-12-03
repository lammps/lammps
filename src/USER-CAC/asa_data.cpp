/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "asa_data.h"
#include "asa_user.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Asa_Data::Asa_Data(LAMMPS *lmp, PairCAC *paircac) : Pointers(lmp)
{
  cgParm=NULL;
  asaParm=NULL;
  Objective=NULL;

  pair_pointer = paircac;
  avec_pointer = NULL;
  class_flag=0;
   
  allocate();
}

/* ---------------------------------------------------------------------- */

Asa_Data::Asa_Data(LAMMPS *lmp, AtomVecCAC *atomveccac) : Pointers(lmp)
{
  cgParm=NULL;
  asaParm=NULL;
  Objective=NULL;

  avec_pointer = atomveccac;
  pair_pointer = NULL;
  class_flag=1;
   
  allocate();
}

/*---------------------------------------------------------------------- */

void Asa_Data::allocate(){

  memory->create(cgParm, 1, "pairCAC:cgParm");

  memory->create(asaParm, 1, "pairCAC:asaParm");

  memory->create(Objective, 1, "pairCAC:asaParm");

  // if you want to change parameter value, initialize strucs with default 
  asa_cg_default(cgParm);
  asa_default(asaParm);

  // if you want to change parameters, change them here: 
  cgParm->PrintParms = FALSE;
  cgParm->PrintLevel = 0;

  asaParm->PrintParms = FALSE;
  asaParm->PrintLevel = 0;
  asaParm->PrintFinal = 0;
  
}

/* ---------------------------------------------------------------------- */

Asa_Data::~Asa_Data()
{
  memory->destroy(cgParm);
  memory->destroy(asaParm);
  memory->destroy(Objective);
}

/* ---------------------------------------------------------------------- */

int Asa_Data::call_asa_cg(double *x,double *lo,double *hi, ASA_INT n,
  double grad_tol, double (*valgrad) (asa_objective *), double *Work, ASA_INT *iWork){
  int flag; //return condition flag
  int max_tries = 100;
  int ntry = 0;
  flag=asa_cg(x, lo, hi, n, NULL, cgParm, asaParm,
    grad_tol, NULL, Work, iWork, this);
  while((flag==4||flag==2)&&ntry!=max_tries){
    ntry++;
    flag=asa_cg(x, lo, hi, n, NULL, cgParm, asaParm,
    grad_tol, NULL, Work, iWork, this);
  }
  if(flag>1){
  if(class_flag==0)
  sprintf(asa_error, "asa_cg iterations failed in finding a solution for quadrature neighboring processes. flag = %d ", flag);
  if(class_flag==1)
  sprintf(asa_error, "asa_cg iterations failed in finding a solution for reneighboring checks. flag = %d ", flag);
  error->one(FLERR,asa_error);
  }
  return flag;

}

/* evaluate the objective function */
/* ---------------------------------------------------------------------- */

double Asa_Data::myvalue(asa_objective *asa) 
{   
  double result;
  if(class_flag==0)
  result = myvalue_surfmin(asa);
  if(class_flag==1)
  result = myvalue_neigh_check(asa);

  return result;
}

/* evaluate the gradient of the objective function */
/* ---------------------------------------------------------------------- */

void Asa_Data::mygrad(asa_objective *asa)
{
  if(class_flag==0)
  mygrad_surfmin(asa);
  if(class_flag==1)
  mygrad_neigh_check(asa);
  
}

/* evaluate the objective function for element to point surface minimization (Q8 element)
/* ---------------------------------------------------------------------- */

double Asa_Data::myvalue_surfmin(asa_objective *asa)
{   
  double f, t, *x;
  double px, py, pz;
  
  double unit_cell_mapped[3];
  double shape_func2;
    double surf_args[3];
  //get necessary data from pair pointer about the current element iteration
  double *quad_r = pair_pointer->quad_r;
  int *neighbor_element_scale = pair_pointer->neighbor_element_scale;
  int *dof_surf_list = pair_pointer->dof_surf_list;
  int *surf_select = pair_pointer->surf_select;
  int neigh_nodes_per_element = pair_pointer->neigh_nodes_per_element;
  int poly_min = pair_pointer->poly_min;
  double **neighbor_element_positions = pair_pointer->neighbor_element_positions[poly_min];
   
  unit_cell_mapped[0] = 2 / double(neighbor_element_scale[0]);
  unit_cell_mapped[1] = 2 / double(neighbor_element_scale[1]);
  unit_cell_mapped[2] = 2 / double(neighbor_element_scale[2]);
  ASA_INT i;
  x = asa->x;
  f = 0;
  
  double r[3] = { quad_r[0]  ,quad_r[1]  , quad_r[2] };

  if (surf_select[0] == 1 && surf_select[1] == -1) {
    surf_args[0] = -1 + unit_cell_mapped[0] / 2;
    surf_args[1] = x[0];
    surf_args[2] = x[1];

  }
  else if (surf_select[0] == 1 && surf_select[1] == 1) {
    surf_args[0] = 1 - unit_cell_mapped[0] / 2;
    surf_args[1] = x[0];
    surf_args[2] = x[1];
  }
  else if (surf_select[0] == 2 && surf_select[1] == -1) {
    surf_args[0] = x[0];
    surf_args[1] = -1 + unit_cell_mapped[1] / 2;
    surf_args[2] = x[1];
  }
  else if (surf_select[0] == 2 && surf_select[1] == 1) {
    surf_args[0] = x[0];
    surf_args[1] = 1 - unit_cell_mapped[1] / 2;
    surf_args[2] = x[1];
  }
  else if (surf_select[0] == 3 && surf_select[1] == -1) {
    surf_args[0] = x[0];
    surf_args[1] = x[1];
    surf_args[2] = -1 + unit_cell_mapped[2] / 2;
  }
  else if (surf_select[0] == 3 && surf_select[1] == 1) {
    surf_args[0] = x[0];
    surf_args[1] = x[1];
    surf_args[2] = 1 - unit_cell_mapped[2] / 2;
  }
  
  px = 0;
  py = 0;
  pz = 0;
  
  for (int kk = 0; kk < neigh_nodes_per_element; kk++) {
    shape_func2 = pair_pointer->shape_function(surf_args[0], surf_args[1], surf_args[2], 2, kk + 1);
    px += neighbor_element_positions[kk][0] * shape_func2;
    py += neighbor_element_positions[kk][1] * shape_func2;
    pz += neighbor_element_positions[kk][2] * shape_func2;
  }


  f = (r[0] - px)*(r[0] - px) + (r[1] - py)*(r[1] - py) + (r[2] - pz)*(r[2] - pz);

  
  return (f);
}

/* evaluate the gradient of the objective function for element to point surface minimization (Q8 element)*/
/* ---------------------------------------------------------------------- */

void Asa_Data::mygrad_surfmin(asa_objective *asa)
{
  double  t, *g, *x;
  double px, py, pz;
  double px1, px2, py1, py2, pz1, pz2;
  double unit_cell_mapped[3];
  double shape_func2, shape_func1;
    double surf_args[3];
  //get necessary data from pair pointer about the current element iteration
  double *quad_r = pair_pointer->quad_r;
  int *neighbor_element_scale = pair_pointer->neighbor_element_scale;
  int *dof_surf_list = pair_pointer->dof_surf_list;
  int *surf_select = pair_pointer->surf_select;
  int neigh_nodes_per_element = pair_pointer->neigh_nodes_per_element;
  int poly_min = pair_pointer->poly_min;
  double **neighbor_element_positions = pair_pointer->neighbor_element_positions[poly_min];
  
  unit_cell_mapped[0] = 2 / double(neighbor_element_scale[0]);
  unit_cell_mapped[1] = 2 / double(neighbor_element_scale[1]);
  unit_cell_mapped[2] = 2 / double(neighbor_element_scale[2]);

  ASA_INT i;
  x = asa->x;
  g = asa->g;
  double r[3] = { quad_r[0]  ,quad_r[1]  , quad_r[2] };
  int deriv_select[2];

  if (surf_select[0] == 1 && surf_select[1] == -1) {
    surf_args[0] = -1 + unit_cell_mapped[0] / 2;
    surf_args[1] = x[0];
    surf_args[2] = x[1];
    deriv_select[0] = 2;
    deriv_select[1] = 3;
  }
  else if (surf_select[0] == 1 && surf_select[1] == 1) {
    surf_args[0] = 1 - unit_cell_mapped[0] / 2;
    surf_args[1] = x[0];
    surf_args[2] = x[1];
    deriv_select[0] = 2;
    deriv_select[1] = 3;
  }
  else if (surf_select[0] == 2 && surf_select[1] == -1) {
    surf_args[0] = x[0];
    surf_args[1] = -1 + unit_cell_mapped[1] / 2;
    surf_args[2] = x[1];
    deriv_select[0] = 1;
    deriv_select[1] = 3;
  }
  else if (surf_select[0] == 2 && surf_select[1] == 1) {
    surf_args[0] = x[0];
    surf_args[1] = 1 - unit_cell_mapped[1] / 2;
    surf_args[2] = x[1];
    deriv_select[0] = 1;
    deriv_select[1] = 3;
  }
  else if (surf_select[0] == 3 && surf_select[1] == -1) {
    surf_args[0] = x[0];
    surf_args[1] = x[1];
    surf_args[2] = -1 + unit_cell_mapped[2] / 2;
    deriv_select[0] = 1;
    deriv_select[1] = 2;
  }
  else if (surf_select[0] == 3 && surf_select[1] == 1) {
    surf_args[0] = x[0];
    surf_args[1] = x[1];
    surf_args[2] = 1 - unit_cell_mapped[2] / 2;
    deriv_select[0] = 1;
    deriv_select[1] = 2;
  }

  px = 0;
  py = 0;
  pz = 0;
  for (int kk = 0; kk < neigh_nodes_per_element; kk++) {
    shape_func2 = pair_pointer->shape_function(surf_args[0], surf_args[1], surf_args[2], 2, kk + 1);

    px += neighbor_element_positions[kk][0] * shape_func2;
    py += neighbor_element_positions[kk][1] * shape_func2;
    pz += neighbor_element_positions[kk][2] * shape_func2;
  }
  
  px1 = 0;
  py1 = 0;
  pz1 = 0;
  px2 = 0;
  py2 = 0;
  pz2 = 0;
  for (int kk = 0; kk < neigh_nodes_per_element; kk++) {
    shape_func1 = pair_pointer->shape_function_derivative(surf_args[0], surf_args[1], surf_args[2], 2, kk + 1, deriv_select[0]);
    shape_func2 = pair_pointer->shape_function_derivative(surf_args[0], surf_args[1], surf_args[2], 2, kk + 1, deriv_select[1]);
    px1 += neighbor_element_positions[kk][0] * shape_func1;
    py1 += neighbor_element_positions[kk][1] * shape_func1;
    pz1 += neighbor_element_positions[kk][2] * shape_func1;
    px2 += neighbor_element_positions[kk][0] * shape_func2;
    py2 += neighbor_element_positions[kk][1] * shape_func2;
    pz2 += neighbor_element_positions[kk][2] * shape_func2;
  }


  g[0] = -2 * (r[0] - px)*px1 - 2 * (r[1] - py)*py1 - 2 * (r[2] - pz)*pz1;
  g[1] = -2 * (r[0] - px)*px2 - 2 * (r[1] - py)*py2 - 2 * (r[2] - pz)*pz2;
  
}

/* evaluate the objective function for the neighbor list decision process
/* ---------------------------------------------------------------------- */
double Asa_Data::myvalue_neigh_check(asa_objective *asa)
{
  double f, t, *x;
  double px, py, pz;
    
  double shape_func2;
    //get class data for current element
  int min_element_index = avec_pointer->min_element_index;
    int *current_element_scale = avec_pointer->check_element_scale[min_element_index];
  int min_nodes_per_element = avec_pointer->min_nodes_per_element;
  int poly_min = avec_pointer->poly_min;
  double **current_hold_positions = avec_pointer->hold_nodal_positions[min_element_index][poly_min];
  double **current_nodal_positions = avec_pointer->check_nodal_positions[min_element_index][poly_min];
  ASA_INT i;
  x = asa->x;
  f = 0;

  px = 0;
  py = 0;
  pz = 0;
  for (int kk = 0; kk < min_nodes_per_element; kk++) {
    shape_func2 = avec_pointer->shape_function(x[0], x[1], x[2], 2, kk + 1);
    px += (current_nodal_positions[kk][0]-current_hold_positions[kk][0]) * shape_func2;
    py += (current_nodal_positions[kk][1]-current_hold_positions[kk][1]) * shape_func2;
    pz += (current_nodal_positions[kk][2]-current_hold_positions[kk][2]) * shape_func2;
  }


  f = -(px*px +  py*py +  pz*pz);


  return (f);
}

/* evaluate the gradient of the objective function for the neighbor list decision process
/* ---------------------------------------------------------------------- */
void Asa_Data::mygrad_neigh_check(asa_objective *asa)
{
  double t, *g, *x;
  double px, py, pz;
  double px1, px2, px3, py1, py2, py3, pz1, pz2, pz3;
  double shape_func3,shape_func2, shape_func1;
    //get class data for current element
  int min_element_index = avec_pointer->min_element_index;
    int *current_element_scale = avec_pointer->check_element_scale[min_element_index];
  int min_nodes_per_element = avec_pointer->min_nodes_per_element;
  int poly_min = avec_pointer->poly_min;
  double **current_hold_positions = avec_pointer->hold_nodal_positions[min_element_index][poly_min];
    double **current_nodal_positions = avec_pointer->check_nodal_positions[min_element_index][poly_min];

  ASA_INT i;
  x = asa->x;
  g = asa->g;


  px = 0;
  py = 0;
  pz = 0;
  for (int kk = 0; kk < min_nodes_per_element; kk++) {
    shape_func2 = avec_pointer->shape_function(x[0], x[1], x[2], 2, kk + 1);
    px += (current_nodal_positions[kk][0]-current_hold_positions[kk][0]) * shape_func2;
    py += (current_nodal_positions[kk][1]-current_hold_positions[kk][1]) * shape_func2;
    pz += (current_nodal_positions[kk][2]-current_hold_positions[kk][2]) * shape_func2;
  }
  
  px1 = 0;
  py1 = 0;
  pz1 = 0;
  px2 = 0;
  py2 = 0;
  pz2 = 0;
  px3 = 0;
  py3 = 0;
  pz3 = 0;
  for (int kk = 0; kk < min_nodes_per_element; kk++) {
    shape_func1 = avec_pointer->shape_function_derivative(x[0], x[1], x[2], 2, kk + 1, 1);
    shape_func2 = avec_pointer->shape_function_derivative(x[0], x[1], x[2], 2, kk + 1, 2);
    shape_func3 = avec_pointer->shape_function_derivative(x[0], x[1], x[2], 2, kk + 1, 3);
    px1 += (current_nodal_positions[kk][0]-current_hold_positions[kk][0]) * shape_func1;
    py1 += (current_nodal_positions[kk][1]-current_hold_positions[kk][1]) * shape_func1;
    pz1 += (current_nodal_positions[kk][2]-current_hold_positions[kk][2]) * shape_func1;
    px2 += (current_nodal_positions[kk][0]-current_hold_positions[kk][0]) * shape_func2;
    py2 += (current_nodal_positions[kk][1]-current_hold_positions[kk][1]) * shape_func2;
    pz2 += (current_nodal_positions[kk][2]-current_hold_positions[kk][2]) * shape_func2;
    px3 += (current_nodal_positions[kk][0]-current_hold_positions[kk][0]) * shape_func3;
    py3 += (current_nodal_positions[kk][1]-current_hold_positions[kk][1]) * shape_func3;
    pz3 += (current_nodal_positions[kk][2]-current_hold_positions[kk][2]) * shape_func3;
  }


  g[0] = -(2 * px*px1 + 2 * py*py1 + 2 *  pz*pz1);
  g[1] = -(2 * px*px2 + 2 * py*py2 + 2 *  pz*pz2);
  g[2] = -(2 * px*px3 + 2 * py*py3 + 2 *  pz*pz3);
}
