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

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include "pair_cac.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "update.h"
#include "neigh_list.h"
#include "integrate.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "domain.h"
#include "asa_user.h"
#include "asa_data.h"

#define MAXESHAPE  10 //maximum number of shape functions per element
#define MAXNEIGHOUT  500
#define MAXNEIGHIN  10
#define MAXNEIGHBUFF 600
#define EXPAND 20
#define MAXLINE 1024
#define DELTA 4
using namespace LAMMPS_NS;


/* ---------------------------------------------------------------------- */

PairCAC::PairCAC(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;
  manybody_flag = 1;
  pre_force_flag = 0;
  nmax = 0;
  nmax_surf = 0;
  cut_global_s = 0;
  quadrature_point_max = 0;
  quadrature_point_data = NULL;
  quadrature_counts = NULL;
  e2quad_index = NULL;
  old_atom_etype = NULL;
  old_all_atom_etype = NULL;
  quad_allocated = 0;
  local_inner_max = 0;
  local_outer_max = 0;
  outer_neighflag = 0;
  sector_flag = 0;
  one_layer_flag = 0;
  shape_quad_result = NULL;
  neighbor->pgsize=10;
  neighbor->oneatom=1;
  old_atom_count=0;
  old_quad_count=0;
  atom->CAC_pair_flag=1;
  //allocate shape function pointers
  shape_functions= (Shape_Functions *) memory->smalloc(sizeof(Shape_Functions)*MAXESHAPE, "Pair CAC:shape_functions");
  set_shape_functions();
}

/* ---------------------------------------------------------------------- */

PairCAC::~PairCAC() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(mass_matrix);
    memory->destroy(quadrature_abcissae);
    memory->destroy(quadrature_weights);
    memory->destroy(force_column);
    memory->destroy(current_force_column);
    memory->destroy(current_nodal_forces);
    memory->destroy(pivot);
    memory->destroy(shape_quad_result);
  }
}

/* ---------------------------------------------------------------------- */

void PairCAC::compute(int eflag, int vflag) {
  int i;
  double delx,dely,delz,fpair;
  int mi;
  int mj;
  int singular;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double ****nodal_positions= atom->nodal_positions;
  double ****nodal_forces= atom->nodal_forces;
  int *element_type = atom->element_type;
  int *poly_count = atom->poly_count;
  int **node_types = atom->node_types;
  int **element_scale = atom->element_scale;
  int nlocal = atom->nlocal;
  int nodes_per_element;
  int nsurface;
  quadrature_counts = atom->quadrature_counts;
  e2quad_index = atom->e2quad_index;
  inner_quad_lists_index = atom->inner_quad_lists_index;
  inner_quad_lists_ucell = atom->inner_quad_lists_ucell;
  inner_quad_lists_counts = atom->inner_quad_lists_counts;
  outer_quad_lists_index = atom->outer_quad_lists_index;
  outer_quad_lists_ucell = atom->outer_quad_lists_ucell;
  outer_quad_lists_counts = atom->outer_quad_lists_counts;
  quadrature_point_data = atom->quadrature_point_data;
  
  int *nodes_count_list = atom->nodes_per_element_list;	
  quad_eflag = eflag;
  cutoff_skin = neighbor->skin;


  pqi = qi = 0;
  //preforce calculations
  if(pre_force_flag)
  pre_force_densities();
 
  pqi = qi = 0;
  //compute forces
  for (i = 0; i < nlocal; i++) {
    //compute this element's mass matrix
    compute_mass_matrix();
    for ( mi = 0; mi < max_nodes_per_element; mi++) {
      for (mj = 0; mj < max_nodes_per_element; mj++){
        mass_copy[mi][mj]=mass_matrix[mi][mj];
      }      
    }

    singular = LUPDecompose(mass_copy, max_nodes_per_element,0.00000000000001,  pivot);
    if(singular==0){
      error->one(FLERR,"LU matrix is degenerate");
    }
    atomic_flag = 0;
    current_element_type = element_type[i];
    current_element_scale = element_scale[i];
    current_poly_count = poly_count[i];
    type_array = node_types[i];
    current_x = x[i];
    if (eflag) {
      element_energy = 0;
    }
    //determine element type
    nodes_per_element = nodes_count_list[current_element_type];
    if (current_element_type == 0) {
      atomic_flag = 1;
    }
    //NOTE:might have to change matrices so they dont have zeros due to maximum node count; ill condition.
    if(atomic_flag){
      poly_counter = 0;
      current_nodal_positions = nodal_positions[i][poly_counter]; 
      compute_forcev(i);
      for (int dim = 0; dim < 3; dim++) {
      nodal_forces[i][0][0][dim] += force_column[0][dim];
      }
    }
    else{
      for (poly_counter = 0; poly_counter < current_poly_count; poly_counter++) {
        current_nodal_positions = nodal_positions[i][poly_counter]; 
        compute_forcev(i);
        for (int dim = 0; dim < 3; dim++) {
          for (mi = 0; mi < nodes_per_element; mi++) {
            current_force_column[mi] = force_column[mi][dim];
          }
          LUPSolve(mass_copy, pivot, current_force_column, nodes_per_element, current_nodal_forces);
          for (mi = 0; mi < nodes_per_element; mi++) {
            nodal_forces[i][poly_counter][mi][dim] += current_nodal_forces[mi];
          }
        }
      }
    }
    if (evflag) ev_tally_full(i,
          2 * element_energy, 0.0, fpair, delx, dely, delz);
  }
}

/* ----------------------------------------------------------------------
 allocate all arrays
 ------------------------------------------------------------------------- */

void PairCAC::allocate()
{
  allocated = 1;
  max_nodes_per_element = atom->nodes_per_element;

  memory->create(mass_matrix, max_nodes_per_element, max_nodes_per_element,"pairCAC:mass_matrix");
  memory->create(mass_copy, max_nodes_per_element, max_nodes_per_element,"pairCAC:copy_mass_matrix");
  memory->create(force_column, max_nodes_per_element,3,"pairCAC:force_residue");
  memory->create(current_force_column, max_nodes_per_element,"pairCAC:current_force_residue");
  memory->create(current_nodal_forces, max_nodes_per_element,"pairCAC:current_nodal_force");
  memory->create(pivot, max_nodes_per_element+1,"pairCAC:pivots");
  quadrature_init(2);
}

/* ----------------------------------------------------------------------
 global settings
 ------------------------------------------------------------------------- */

void PairCAC::settings(int narg, char **arg) {
  if (narg>1) error->all(FLERR,"Illegal CAC pair_style command");
 
  if(narg==1){
    if (strcmp(arg[0], "one") == 0) atom->one_layer_flag=one_layer_flag = 1;
    else error->all(FLERR, "Unexpected argument in CAC pair style invocation");
  }
  
  force->newton_pair=0;
}

/* ---------------------------------------------------------------------- */

void PairCAC::init_style()
{
  check_existence_flags();
  atom->max_neigh_inner_init = maxneigh_quad_inner = MAXNEIGHIN;
  atom->max_neigh_outer_init = maxneigh_quad_outer = MAXNEIGHOUT;
  int irequest = neighbor->request(this, instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->cac = 1;
}

//-----------------------------------------------------------------------

void PairCAC::check_existence_flags(){
//check existence of CAC atomvec here for convenience since children implement init_style
//but all implementations call this routine
if(!atom->CAC_flag)
error->all(FLERR,"CAC Pair style requires a CAC atom style");
//check if hybrid pair style is invoked
if(force->pair_match("hybrid",0,0)!=NULL)
error->all(FLERR,"CAC Pair styles cannot be invoked with hybrid; also don't use the word hybrid in your style name to avoid this error");
//check that the CAC comm style is defined
if (strcmp(comm->comm_style, "cac") != 0)
error->all(FLERR,"CAC Pair style requires a CAC comm style");
}

//-----------------------------------------------------------------------

void PairCAC::quadrature_init(int quadrature_rank){

  if(quadrature_rank==1){
  atom->quadrature_node_count=quadrature_node_count=1;
  memory->create(quadrature_weights,quadrature_node_count,"pairCAC:quadrature_weights");
  memory->create(quadrature_abcissae,quadrature_node_count,"pairCAC:quadrature_abcissae");
  quadrature_weights[0]=2;
  quadrature_abcissae[0]=0;
  }
  if(quadrature_rank==2){
  atom->quadrature_node_count=quadrature_node_count=2;
  memory->create(quadrature_weights,quadrature_node_count,"pairCAC:quadrature_weights");
  memory->create(quadrature_abcissae,quadrature_node_count,"pairCAC:quadrature_abcissae");
  quadrature_weights[0]=1;
  quadrature_weights[1]=1;
  quadrature_abcissae[0]=-0.5773502691896258;
  quadrature_abcissae[1]=0.5773502691896258;
  }
  if(quadrature_rank==3)
  {

  }
  if(quadrature_rank==4)
  {

  }
  if(quadrature_rank==5)
  {

  }

  memory->create(shape_quad_result, max_nodes_per_element, quadrature_node_count*quadrature_node_count*quadrature_node_count, "pairCAC:shape_quad_result");
}

//---------------------------------------------------------------------

//--------------------------------------------------------------------

void PairCAC::compute_mass_matrix()
{
  int j,k;
  double result=0;
  //precompute shape function quadrature abcissae
  for(int ii=0;ii<max_nodes_per_element;ii++){
    for (int i=0; i< quadrature_node_count;i++){
      for (int j=0; j< quadrature_node_count;j++){
        for ( k=0; k< quadrature_node_count;k++){
          shape_quad_result[ii][quadrature_node_count*quadrature_node_count*(i)+quadrature_node_count*(j)+k]=
          shape_function(quadrature_abcissae[i],quadrature_abcissae[j],quadrature_abcissae[k],2,ii+1);
        }
      }
    }
  }

//assemble matrix with quadrature array of values
  for ( j=0; j<max_nodes_per_element;j++){
    for ( k=j; k<max_nodes_per_element;k++){
      result=shape_product(j,k);
      mass_matrix[j][k]= result;
      if(k>j){
      mass_matrix[k][j]= result;
      }
    }
  }
}

double PairCAC::shape_product(int ii, int jj) {
  double result = 0;
  for (int i = 0; i< quadrature_node_count; i++) {
    for (int j = 0; j< quadrature_node_count; j++) {
      for (int k = 0; k< quadrature_node_count; k++) {
        result += quadrature_weights[i] * quadrature_weights[j] * quadrature_weights[k] *
          shape_quad_result[ii][quadrature_node_count*quadrature_node_count*(i)+quadrature_node_count*(j)+k] *
          shape_quad_result[jj][quadrature_node_count*quadrature_node_count*(i)+quadrature_node_count*(j)+k];
      }
    }
  }
  return result;
}


//----------------------------------------------------------------

void PairCAC::compute_forcev(int iii){
  double unit_cell_mapped[3];
  int *nodes_count_list = atom->nodes_per_element_list;	
  double coefficients;
  int nodes_per_element;
  double ****nodal_virial = atom->nodal_virial;
  double shape_func;
  //debug array 
  //stores neighbor information for computational weight assignment
  
  unit_cell_mapped[0] = 2 / double(current_element_scale[0]);
  unit_cell_mapped[1] = 2 / double(current_element_scale[1]);
  unit_cell_mapped[2] = 2 / double(current_element_scale[2]);
  
  double s, t, w;
  double sq, tq, wq;
  //compute interior quadrature point virtual neighbor lists
  double iso_volume=unit_cell_mapped[0]*unit_cell_mapped[1]*unit_cell_mapped[2];
  if(atomic_flag) iso_volume=1;
  nodes_per_element = nodes_count_list[current_element_type];
  for (int js = 0; js<nodes_per_element; js++) {
    for (int jj = 0; jj<3; jj++) {
      force_column[js][jj] = 0;
    }
  }

  double force_density[3];
  
    //sum over quadrature points to compute force density
   
  for (int quad_loop=0; quad_loop < quadrature_counts[iii] ; quad_loop++){
    pqi = e2quad_index[iii]+current_poly_count*quad_loop+poly_counter;
    quadrature_energy = 0;
    force_density[0] = 0;
    force_density[1] = 0;
    force_density[2] = 0;
    if(atom->CAC_virial){
    virial_density[0] = 0;
    virial_density[1] = 0;
    virial_density[2] = 0;
    virial_density[3] = 0;
    virial_density[4] = 0;
    virial_density[5] = 0;
    }
    if(!atomic_flag){
    s = quadrature_point_data[qi + quad_loop][0];
    t = quadrature_point_data[qi + quad_loop][1];
    w = quadrature_point_data[qi + quad_loop][2];
    sq = quadrature_point_data[qi + quad_loop][3];
    tq = quadrature_point_data[qi + quad_loop][4];
    wq = quadrature_point_data[qi + quad_loop][5];
  }
  coefficients = quadrature_point_data[qi + quad_loop][6];
  if(!atomic_flag)
  force_densities(iii, s, t, w, coefficients,
    force_density[0], force_density[1], force_density[2]);
  else
  force_densities(iii, current_x[0], current_x[1], current_x[2], coefficients,
    force_density[0], force_density[1], force_density[2]);
  if(!atomic_flag){
  for (int js = 0; js < nodes_per_element; js++) {
    shape_func = shape_function(sq, tq, wq, 2, js + 1);
    for (int jj = 0; jj < 3; jj++) {
    force_column[js][jj] += coefficients*force_density[jj] * shape_func;
    }
    if(atom->CAC_virial){
    for (int jj = 0; jj < 6; jj++)
    nodal_virial[iii][poly_counter][js][jj] += coefficients*virial_density[jj] * shape_func;
    }
  }
  }
  else{
    for (int jj = 0; jj < 3; jj++) {
    force_column[0][jj] = force_density[jj];
    }	
    if(atom->CAC_virial){
    for (int jj = 0; jj < 6; jj++) {
    nodal_virial[iii][0][0][jj] += virial_density[jj];
    }
    }
  }
  if (quad_eflag) {
  element_energy += coefficients*quadrature_energy/iso_volume;
  }
  }
  if(poly_counter == current_poly_count - 1) qi += quadrature_counts[iii];
}


//-------------------------------------------------------------------------

//3by3 solver
int PairCAC::mldivide3(const double m[3][3], const double *v, double *ans)
{
  // create augmented matrix for pivoting

  double aug[3][4];
  for (unsigned i = 0; i < 3; i++) {
    aug[i][3] = v[i];
    for (unsigned j = 0; j < 3; j++) aug[i][j] = m[i][j];
  }

  for (unsigned i = 0; i < 2; i++) {
    unsigned p = i;
    for (unsigned j = i+1; j < 3; j++) {
      if (fabs(aug[j][i]) > fabs(aug[i][i])) {
        double tempv[4];
        memcpy(tempv,aug[i],4*sizeof(double));
        memmove(aug[i],aug[j],4*sizeof(double));
        memcpy(aug[j],tempv,4*sizeof(double));
      }
    }

    while (p < 3 && aug[p][i] == 0.0) p++;

    if (p == 3) return 1;
    else
      if (p != i) {
        double tempv[4];
        memcpy(tempv,aug[i],4*sizeof(double));
        memmove(aug[i],aug[p],4*sizeof(double));
        memcpy(aug[p],tempv,4*sizeof(double));
      }

    for (unsigned j = i+1; j < 3; j++) {
      double m = aug[j][i]/aug[i][i];
      for (unsigned k=i+1; k<4; k++) aug[j][k]-=m*aug[i][k];
    }
  }

  if (aug[2][2] == 0.0) return 1;

  // back substitution

  ans[2] = aug[2][3]/aug[2][2];
  for (int i = 1; i >= 0; i--) {
    double sumax = 0.0;
    for (unsigned j = i+1; j < 3; j++) sumax += aug[i][j]*ans[j];
    ans[i] = (aug[i][3]-sumax) / aug[i][i];
  }

  return 0;
}

//--------------------------------------------------------------------------
//populate array of shape functions
void PairCAC::set_shape_functions(){
shape_functions[0]=&PairCAC::quad_shape_one;
shape_functions[1]=&PairCAC::quad_shape_two;
shape_functions[2]=&PairCAC::quad_shape_three;
shape_functions[3]=&PairCAC::quad_shape_four;
shape_functions[4]=&PairCAC::quad_shape_five;
shape_functions[5]=&PairCAC::quad_shape_six;
shape_functions[6]=&PairCAC::quad_shape_seven;
shape_functions[7]=&PairCAC::quad_shape_eight;
}

////-------------------------------------------------------------------------

void PairCAC::interpolation(int iii){

  double s, t, w, sp, sm, tp, tm, wp, wm;
  int *element_type = atom->element_type;
  double ****nodal_positions = atom->nodal_positions;
  double **x = atom->x;
  int neigh_max = inner_quad_lists_counts[pqi];
  int outer_neigh_max;
  double **outer_ucells;
  int outer_etypes[MAXNEIGHBUFF];
  int **outer_indices;
  if(neigh_max > MAXNEIGHBUFF)
  error->one(FLERR,"increase MAXNEIGHBUFF in pair_cac.cpp; this overestimates the maximum possible neighbors for a quadrature point");
  
  double **inner_ucells;
  int inner_etypes[MAXNEIGHBUFF];
  int **inner_indices;
  if(outer_neighflag){
  outer_neigh_max = outer_quad_lists_counts[pqi];
  if(outer_neigh_max > MAXNEIGHBUFF)
  error->one(FLERR,"increase MAXNEIGHBUFF in pair_cac.cpp; this overestimates the maximum possible neighbors for a quadrature point");	  
  }

  int shape_limit[MAXNEIGHBUFF];
  double shaperesult[MAXNEIGHBUFF][MAXESHAPE];
  double **nodes;
  inner_ucells = inner_quad_lists_ucell[pqi];
  inner_indices = inner_quad_lists_index[pqi];
  //init coords to zero etc.
  for (int l = 0; l < neigh_max; l++) {
  inner_neighbor_coords[l][0] = 0;
  inner_neighbor_coords[l][1] = 0;
  inner_neighbor_coords[l][2] = 0;
  }
  //consider elements
  //accumulate interpolation to coordinates
  for (int l = 0; l < neigh_max; l++) {
  inner_etypes[l] = element_type[inner_indices[l][0]];
  //Q8 interpolation scheme
  if(inner_etypes[l]==1){
  s = inner_ucells[l][0];
  t = inner_ucells[l][1];
  w = inner_ucells[l][2];
  sp = (1+s)/8;
  sm = (1-s)/8;
  tp = 1+t;
  tm = 1-t;
  wp = 1+w;
  wm = 1-w;
    nodes = nodal_positions[inner_indices[l][0]][inner_indices[l][1]];
    shaperesult[l][0]= sm*tm*wm;
    shaperesult[l][1]= sp*tm*wm;
    shaperesult[l][2]= sp*tp*wm;
    shaperesult[l][3]= sm*tp*wm;
    shaperesult[l][4]= sm*tm*wp;
    shaperesult[l][5]= sp*tm*wp;
    shaperesult[l][6]= sp*tp*wp;
    shaperesult[l][7]= sm*tp*wp;
  //#pragma vector always
  for(int r = 0; r < 8; r++){
    inner_neighbor_coords[l][0]+=nodes[r][0]*shaperesult[l][r];
    inner_neighbor_coords[l][1]+=nodes[r][1]*shaperesult[l][r];
    inner_neighbor_coords[l][2]+=nodes[r][2]*shaperesult[l][r];
  }
  }
  //add new shape function block here
  if (inner_etypes[l]==2){
  }
  //consider atoms
  if(inner_etypes[l]==0){	
  inner_neighbor_coords[l][0] = x[inner_indices[l][0]][0];
  inner_neighbor_coords[l][1] = x[inner_indices[l][0]][1];
  inner_neighbor_coords[l][2] = x[inner_indices[l][0]][2];
  }

  }

  if(outer_neighflag){
  outer_ucells = outer_quad_lists_ucell[pqi];
  outer_indices = outer_quad_lists_index[pqi];
  //init coords to zero etc.
  for (int l = 0; l < outer_neigh_max; l++) {
  outer_neighbor_coords[l][0] = 0;
  outer_neighbor_coords[l][1] = 0;
  outer_neighbor_coords[l][2] = 0;
  }
  //consider elements
  //accumulate interpolation to coordinates
  for (int l = 0; l < outer_neigh_max; l++) {
  outer_etypes[l] = element_type[outer_indices[l][0]];
  if(outer_etypes[l]==1){
  s = outer_ucells[l][0];
  t = outer_ucells[l][1];
  w = outer_ucells[l][2];
  sp = (1+s)/8;
  sm = (1-s)/8;
  tp = 1+t;
  tm = 1-t;
  wp = 1+w;
  wm = 1-w;
    nodes = nodal_positions[outer_indices[l][0]][outer_indices[l][1]];
    shaperesult[l][0]= sm*tm*wm;
    shaperesult[l][1]= sp*tm*wm;
    shaperesult[l][2]= sp*tp*wm;
    shaperesult[l][3]= sm*tp*wm;
    shaperesult[l][4]= sm*tm*wp;
    shaperesult[l][5]= sp*tm*wp;
    shaperesult[l][6]= sp*tp*wp;
    shaperesult[l][7]= sm*tp*wp;
  //#pragma vector always
  for(int r = 0; r < 8; r++){
  outer_neighbor_coords[l][0]+=nodes[r][0]*shaperesult[l][r];
  outer_neighbor_coords[l][1]+=nodes[r][1]*shaperesult[l][r];
  outer_neighbor_coords[l][2]+=nodes[r][2]*shaperesult[l][r];
  }
  }
  //add new shape function block here
  if (outer_etypes[l]==2){
  }
  //consider atoms
  if(outer_etypes[l]==0){
  outer_neighbor_coords[l][0] = x[outer_indices[l][0]][0];
  outer_neighbor_coords[l][1] = x[outer_indices[l][0]][1];
  outer_neighbor_coords[l][2] = x[outer_indices[l][0]][2];
  }

  }	  
  }
  
}

//-------------------------------------------------------------------------

double PairCAC::shape_function(double s, double t, double w, int flag, int index){
double shape_function;
if(flag==2){
    if(index==1){
    shape_function=(1-s)*(1-t)*(1-w)/8;
    }
    else if(index==2){
    shape_function=(1+s)*(1-t)*(1-w)/8;
    }
    else if(index==3){
    shape_function=(1+s)*(1+t)*(1-w)/8;
    }
    else if(index==4){
    shape_function=(1-s)*(1+t)*(1-w)/8;
    }
    else if(index==5){
    shape_function=(1-s)*(1-t)*(1+w)/8;
    }
    else if(index==6){
    shape_function=(1+s)*(1-t)*(1+w)/8;
    }
    else if(index==7){
    shape_function=(1+s)*(1+t)*(1+w)/8;
    }
    else if(index==8){
    shape_function=(1-s)*(1+t)*(1+w)/8;
    }


}
return shape_function;


}

//-------------------------------------------------------------------------

double PairCAC::shape_function_derivative(double s, double t, double w, int flag, int index, int derivative){
double shape_function;
//flag determines the element type and corresponding basis/shape functions
if(flag==2){

 if(derivative==1){
    if(index==1){
    shape_function=-(1-t)*(1-w)/8;
    }
    else if(index==2){
    shape_function=(1-t)*(1-w)/8;
    }
    else if(index==3){
    shape_function=(1+t)*(1-w)/8;
    }
    else if(index==4){
    shape_function=-(1+t)*(1-w)/8;
    }
    else if(index==5){
    shape_function=-(1-t)*(1+w)/8;
    }
    else if(index==6){
    shape_function=(1-t)*(1+w)/8;
    }
    else if(index==7){
    shape_function=(1+t)*(1+w)/8;
    }
    else if(index==8){
    shape_function=-(1+t)*(1+w)/8;
    }
 }
 else if(derivative==2){
            if(index==1){
    shape_function=-(1-s)*(1-w)/8;
    }
    else if(index==2){
    shape_function=-(1+s)*(1-w)/8;
    }
    else if(index==3){
    shape_function=(1+s)*(1-w)/8;
    }
    else if(index==4){
    shape_function=(1-s)*(1-w)/8;
    }
    else if(index==5){
    shape_function=-(1-s)*(1+w)/8;
    }
    else if(index==6){
    shape_function=-(1+s)*(1+w)/8;
    }
    else if(index==7){
    shape_function=(1+s)*(1+w)/8;
    }
    else if(index==8){
    shape_function=(1-s)*(1+w)/8;
    }

 }
 else if(derivative==3){
        if(index==1){
    shape_function=-(1-s)*(1-t)/8;
    }
    else if(index==2){
    shape_function=-(1+s)*(1-t)/8;
    }
    else if(index==3){
    shape_function=-(1+s)*(1+t)/8;
    }
    else if(index==4){
    shape_function=-(1-s)*(1+t)/8;
    }
    else if(index==5){
    shape_function=(1-s)*(1-t)/8;
    }
    else if(index==6){
    shape_function=(1+s)*(1-t)/8;
    }
    else if(index==7){
    shape_function=(1+s)*(1+t)/8;
    }
    else if(index==8){
    shape_function=(1-s)*(1+t)/8;
    }
 }
}

return shape_function;

}

//////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// File: doolittle_pivot.c                                                    //
// Routines:                                                                  //
//    Doolittle_LU_Decomposition_with_Pivoting                                //
//    Doolittle_LU_with_Pivoting_Solve                                        //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  int Doolittle_LU_Decomposition_with_Pivoting(double *A, int pivot[],      //
//                                                                    int n)  //
//                                                                            //
//  Description:                                                              //
//     This routine uses Doolittle's method to decompose the n x n matrix A   //
//     into a unit lower triangular matrix L and an upper triangular matrix U //
//     such that A = LU.                                                      //
//     The matrices L and U replace the matrix A so that the original matrix  //
//     A is destroyed.  Note!  In Doolittle's method the diagonal elements of //
//     L are 1 and are not stored.                                            //
//     The LU decomposition is convenient when one needs to solve the linear  //
//     equation Ax = B for the vector x while the matrix A is fixed and the   //
//     vector B is varied.  The routine for solving the linear system Ax = B  //
//     after performing the LU decomposition for A is                         //
//                      Doolittle_LU_with_Pivoting_Solve.                     //
//                                                                            //
//     The Doolittle method with partial pivoting is:  Determine the pivot    //
//     row and interchange the current row with the pivot row, then assuming  //
//     that row k is the current row, k = 0, ..., n - 1 evaluate in order the //
//     following pair of expressions                                          //
//       U[k][j] = A[k][j] - (L[k][0]*U[0][j] + ... + L[k][k-1]*U[k-1][j])    //
//                                 for j = k, k+1, ... , n-1                  //
//       L[i][k] = (A[i][k] - (L[i][0]*U[0][k] + . + L[i][k-1]*U[k-1][k]))    //
//                          / U[k][k]                                         //
//                                 for i = k+1, ... , n-1.                    //
//       The matrix U forms the upper triangular matrix, and the matrix L     //
//       forms the lower triangular matrix.                                   //
//                                                                            //
//  Arguments:                                                                //
//     double *A       Pointer to the first element of the matrix A[n][n].    //
//     int    pivot[]  The i-th element is the pivot row interchanged with    //
//                     row i.                                                 //
//     int     n       The number of rows or columns of the matrix A.         //
//                                                                            //
//  Return Values:                                                            //
//     0  Success                                                             //
//    -1  Failure - The matrix A is singular.                                 //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double A[N][N];                                                        //
//     int    pivot[N];                                                       //
//                                                                            //
//     (your code to create matrix A)                                         //
//     err = Doolittle_LU_Decomposition_with_Pivoting(&A[0][0], pivot, N);    //
//     if (err < 0) printf(" Matrix A is singular\n");                        //
//     else { printf(" The LU decomposition of A is \n");                     //
//           ...                                                              //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //


int PairCAC::LUPDecompose(double **A, int N, double Tol, int *P) {

    int i, j, k, imax;
    double maxA, *ptr, absA;

    for (i = 0; i <= N; i++)
        P[i] = i; //Unit permutation matrix, P[N] initialized with N

    for (i = 0; i < N; i++) {
        maxA = 0.0;
        imax = i;

        for (k = i; k < N; k++)
            if ((absA = fabs(A[k][i])) > maxA) {
                maxA = absA;
                imax = k;
            }

        if (maxA < Tol) {return 0; }//failure, matrix is degenerate

        if (imax != i) {
            //pivoting P
            j = P[i];
            P[i] = P[imax];
            P[imax] = j;

            //pivoting rows of A
            ptr = A[i];
            A[i] = A[imax];
            A[imax] = ptr;

            //counting pivots starting from N (for determinant)
            P[N]++;
        }

        for (j = i + 1; j < N; j++) {
            A[j][i] /= A[i][i];

            for (k = i + 1; k < N; k++)
                A[j][k] -= A[j][i] * A[i][k];
        }
    }

    return 1;  //decomposition done
}


/* INPUT: A,P filled in LUPDecompose; b - rhs vector; N - dimension
 * OUTPUT: x - solution vector of A*x=b
 */
void PairCAC::LUPSolve(double **A, int *P, double *b, int N, double *x) {

    for (int i = 0; i < N; i++) {
        x[i] = b[P[i]];

        for (int k = 0; k < i; k++)
            x[i] -= A[i][k] * x[k];
    }

    for (int i = N - 1; i >= 0; i--) {
        for (int k = i + 1; k < N; k++)
            x[i] -= A[i][k] * x[k];

        x[i] = x[i] / A[i][i];
    }
}


//memory usage due to quadrature point list memory structure
double PairCAC::memory_usage()
{
  double bytes_used = 0;
  return bytes_used;
}