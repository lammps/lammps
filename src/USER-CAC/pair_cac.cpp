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
#include "fix_cac_alloc_vector.h"
#include "fix_cac_flux_check.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "domain.h"
#include "asa_user.h"
#include "asa_data.h"

#define MAXESHAPE  10 //maximum number of shape functions per element
#define PLANE_EPSILON  1e-6 //error tolerance for surface flux calculation
#define MAXNEIGHOUT  500
#define MAXNEIGHIN  10
#define MAXNEIGHBUFF 6000
#define EXPAND 10
#define MAXLINE 1024
#define DELTA 4
using namespace LAMMPS_NS;


/* ---------------------------------------------------------------------- */

PairCAC::PairCAC(LAMMPS *lmp) : Pair(lmp)
{
  restartinfo = 0;
  manybody_flag = 1;
  pre_force_flag = 0;
  flux_enable = type_map = 0;
  nmax_surf = 0;
  cut_global_s = 0;
  quadrature_point_max = 0;
  quadrature_point_data = NULL;
  quadrature_counts = NULL;
  e2quad_index = NULL;
  old_atom_etype = NULL;
  old_all_atom_etype = NULL;
  flux_neigh_indices = NULL;
  flux_neigh_intercepts = NULL;
  flux_neigh_nintercepts = NULL;
  nplane_intersects = NULL;
  inner_neighbor_coords = NULL;
  inner_neighbor_types = NULL;
  inner_neighbor_charges = NULL;
  inner_neighbor_velocities = NULL;
  outer_neighbor_coords = NULL;
  outer_neighbor_types = NULL;
  outer_neighbor_charges = NULL;
  outer_neighbor_velocities = NULL;
  add_neighbor_coords = NULL;
  add_neighbor_types = NULL;
  add_neighbor_charges = NULL;
  add_neighbor_velocities = NULL;
  type_maps = NULL;
  quad_allocated = 0;
  local_inner_max = local_all_max = local_outer_max = local_add_max = 0;
  vel_inner_max = vel_outer_max = 0;
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
  double delx,dely,delz,fpair, v[6];
  int mi;
  int mj;
  int singular;
  ev_init(eflag,vflag);
  v[0]=0;v[1]=0;v[2]=0;v[3]=0;v[4]=0;v[5]=0;

  double **x = atom->x;
  double ****nodal_positions= atom->nodal_positions;
  double ****nodal_forces= atom->nodal_forces;
  double ****nodal_virial= atom->nodal_virial;
  int *element_type = atom->element_type;
  int *poly_count = atom->poly_count;
  int **node_types = atom->node_types;
  int **element_scale = atom->element_scale;
  int nlocal = atom->nlocal;
  int nodes_per_element;
  int nsurface;
  flux_compute = atom->flux_compute;
  quadrature_counts = atom->quadrature_counts;
  e2quad_index = atom->e2quad_index;
  inner_quad_lists_index = atom->inner_quad_lists_index;
  inner_quad_lists_ucell = atom->inner_quad_lists_ucell;
  inner_quad_lists_counts = atom->inner_quad_lists_counts;
  outer_quad_lists_index = atom->outer_quad_lists_index;
  outer_quad_lists_ucell = atom->outer_quad_lists_ucell;
  outer_quad_lists_counts = atom->outer_quad_lists_counts;
  quadrature_point_data = atom->quadrature_point_data;
  if(flux_compute){
    add_quad_lists_index = atom->add_quad_lists_index;
    add_quad_lists_ucell = atom->add_quad_lists_ucell;
    add_quad_lists_counts = atom->add_quad_lists_counts;
  }
  
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
    if(element_type[i]){
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
    //Tally energy
    if (evflag) ev_tally_full(i,
          2 * element_energy, 0.0, fpair, delx, dely, delz);

    //Tally virial, this requires looping over all nodes and atoms
    if(vflag){
      for (poly_counter = 0; poly_counter < current_poly_count; poly_counter++) {
            for (mi = 0; mi < nodes_per_element; mi++){
                for( int dim = 0; dim < 6; dim ++) v[dim] = nodal_virial[i][poly_counter][mi][dim];
                // For elements we have to multiply by the number of atoms assigned to node. Just assume that all the elements are rhombohedral for right now
                if(element_type[i] != 0){
                    for(int dim = 0; dim < 6; dim ++) v[dim] *= (element_scale[i][0]*element_scale[i][1]*element_scale[i][2])/nodes_count_list[element_type[i]];
                }
                v_tally_tensor(0,0,1,0,v[0],v[1],v[2],v[3],v[4],v[5]);
            }
      }
    }
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
  int irequest = neighbor->request(this, instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->cac = 1;

  // create fix needed for storing nodal flux vector if needed
  // deleted at the end of run by modify
  if(atom->cac_flux_flag){
    char **fixarg = new char*[4];
    fixarg[0] = (char *) "CACNODALFLUX";
    fixarg[1] = (char *) "all";
    fixarg[2] = (char *) "CAC/ALLOCVECTOR";
    fixarg[3] = (char *) "3";
    modify->add_fix(4,fixarg);
    
    
    fix_cac_alloc_vector = (FixCACAllocVector *) modify->fix[modify->nfix-1];

    fix_cac_alloc_vector->add_vector(24);
    atom->nodal_fluxes = fix_cac_alloc_vector->request_vector(0);

    fixarg[0] = (char *) "CACFLUXCHECK";
    fixarg[1] = (char *) "all";
    fixarg[2] = (char *) "CAC/FLUXCHECK";
    modify->add_fix(3,fixarg);
    delete [] fixarg;
  }

  
}

/* ---------------------------------------------------------------------- */

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
  if(atom->cac_flux_flag&&!flux_enable)
    error->all(FLERR,"This CAC Pair style does not support a surface averaged flux calculation");
}

/* ---------------------------------------------------------------------- */

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

/* ---------------------------------------------------------------------- */

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

/* ---------------------------------------------------------------------- */

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


/* ---------------------------------------------------------------------- */

void PairCAC::compute_forcev(int iii){
  double unit_cell_mapped[3];
  int *nodes_count_list = atom->nodes_per_element_list;	
  double coefficients;
  int nodes_per_element, ns;
  int **surface_counts = atom->surface_counts;
  double ****nodal_virial = atom->nodal_virial;
  double ****nodal_fluxes;
  double shape_func;
  if(flux_compute){
    nodal_fluxes = fix_cac_alloc_vector->request_vector(0);
    atom->nodal_fluxes = nodal_fluxes;
  }
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
    if(atom->cac_flux_flag==2&&flux_compute) quad_flux_flag = 1;
    else if(atom->cac_flux_flag==1&&flux_compute){
      //check to make sure this quadrature point is the closest one to a node
      if(quadrature_point_data[qi + quad_loop][7]) quad_flux_flag = 1;
      else quad_flux_flag = 0;
    }
    else quad_flux_flag = 0;
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
    if(quad_flux_flag)
      for(int init = 0; init < 24; init++)
        flux_density[init] = 0;
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
      for (int jj = 0; jj < 6; jj++){
          nodal_virial[iii][poly_counter][js][jj] += coefficients*virial_density[jj] * shape_func;
        }
    }
    if(quad_flux_flag&&atom->cac_flux_flag==2){
    for (int jj = 0; jj < 24; jj++)
      nodal_fluxes[iii][poly_counter][js][jj] += coefficients*flux_density[jj] * shape_func;
    }
  }
    if(quad_flux_flag&&atom->cac_flux_flag==1){
      ns = quadrature_point_data[qi + quad_loop][7];
      for (int jj = 0; jj < 24; jj++)
        nodal_fluxes[iii][poly_counter][ns][jj] = flux_density[jj];
    }
  }
  else{
    for (int jj = 0; jj < 3; jj++) {
      force_column[0][jj] = force_density[jj];
    }	
    if(atom->CAC_virial){
    for (int jj = 0; jj < 6; jj++) {
      nodal_virial[iii][0][0][jj] = virial_density[jj];
    }
    }
    if(quad_flux_flag){
    for (int jj = 0; jj < 24; jj++) {
      nodal_fluxes[iii][0][0][jj] = flux_density[jj];
    }
    }
  }
  if (quad_eflag) {
  element_energy += coefficients*quadrature_energy/iso_volume;
  }
  }
  if(poly_counter == current_poly_count - 1) qi += quadrature_counts[iii];
}


/* ---------------------------------------------------------------------- 
 3by3 solver
---------------------------------------------------------------------- */

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

/* ---------------------------------------------------------------------- */
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

/* ---------------------------------------------------------------------- */

void PairCAC::allocate_quad_memory(){

  int neigh_max_inner = inner_quad_lists_counts[pqi];
  int neigh_max_outer,neigh_max_add, neigh_max_all;
  if(outer_neighflag) neigh_max_outer = outer_quad_lists_counts[pqi];
  if(quad_flux_flag) {
    neigh_max_add = add_quad_lists_counts[pqi];
    if(outer_neighflag) neigh_max_all = neigh_max_inner + neigh_max_outer + neigh_max_add;
    else neigh_max_all = neigh_max_inner + neigh_max_add;
  }

  if(neigh_max_inner>local_inner_max){
    memory->grow(inner_neighbor_types, neigh_max_inner+EXPAND, "Pair_CAC:inner_neighbor_types");
    memory->grow(inner_neighbor_coords, neigh_max_inner+EXPAND, 3, "Pair_CAC:inner_neighbor_coords");
    if(atom->q_flag)
    memory->grow(inner_neighbor_charges, neigh_max_inner+EXPAND, "Pair_CAC:neighbor_charges");
    local_inner_max=neigh_max_inner+EXPAND;
  }
  if(neigh_max_inner>vel_inner_max&&quad_flux_flag){
      memory->grow(inner_neighbor_velocities, neigh_max_inner+EXPAND, 3, "Pair_CAC:inner_neighbor_velocities");
      vel_inner_max=neigh_max_inner+EXPAND;
    }
  if(outer_neighflag){
    if(neigh_max_outer>local_outer_max){
      memory->grow(outer_neighbor_coords, neigh_max_outer+EXPAND, 3, "Pair_CAC:outer_neighbor_coords");
      memory->grow(outer_neighbor_types, neigh_max_outer+EXPAND, "Pair_CAC:outer_neighbor_types");
      if(atom->q_flag)
      memory->grow(outer_neighbor_charges, neigh_max_outer+EXPAND, "Pair_CAC:outer_neighbor_charges");
      local_outer_max=neigh_max_outer+EXPAND;
    }
    if(neigh_max_outer>vel_outer_max&&quad_flux_flag){
      memory->grow(outer_neighbor_velocities, neigh_max_outer+EXPAND, 3, "Pair_CAC:outer_neighbor_velocities");
      vel_outer_max=neigh_max_outer+EXPAND;
    }
  }
  if(quad_flux_flag){
    //allocate material data arrays
    if(neigh_max_add>local_add_max){
      memory->grow(add_neighbor_coords, neigh_max_add+EXPAND, 3, "Pair_CAC:add_neighbor_coords");
      memory->grow(add_neighbor_types, neigh_max_add+EXPAND, "Pair_CAC:add_neighbor_types");
      if(atom->q_flag)
      memory->grow(outer_neighbor_charges, neigh_max_add+EXPAND, "Pair_CAC:neighbor_charges");
      memory->grow(add_neighbor_velocities, neigh_max_add+EXPAND, 3, "Pair_CAC:add_neighbor_velocities");
      local_add_max=neigh_max_add+EXPAND;
    }

    //allocate line-plane intersection arrays
    if(neigh_max_all>local_all_max){
      flux_neigh_nintercepts = 
        memory->grow(flux_neigh_nintercepts, neigh_max_all+EXPAND, "Pair CAC:flux_neigh_nintercepts");
      flux_neigh_indices = (int **) 
        memory->srealloc(flux_neigh_indices, (neigh_max_all+EXPAND)*sizeof(int *), "Pair CAC:flux_neigh_indices");
      nplane_intersects = (int **) 
        memory->srealloc(nplane_intersects, (neigh_max_all+EXPAND)*sizeof(int *), "Pair CAC:nplane_intersects");
      flux_neigh_intercepts = (int ***) 
        memory->srealloc(flux_neigh_intercepts, (neigh_max_all+EXPAND)*sizeof(int **), "Pair CAC:flux_neigh_intercepts");
      //init new array elements to NULL
      for (int alloc=local_all_max; alloc < neigh_max_all + EXPAND; alloc++){
        flux_neigh_indices[alloc] = NULL;
        nplane_intersects[alloc] = NULL;
        flux_neigh_intercepts[alloc] = NULL;
      }
      for (int alloc=0; alloc < neigh_max_all + EXPAND; alloc++){
        memory->grow(flux_neigh_indices[alloc], neigh_max_all+EXPAND, "Pair CAC:flux_neigh_indices[alloc]");
        memory->grow(nplane_intersects[alloc], neigh_max_all+EXPAND, "Pair CAC:nplane_intersects[alloc]");
        memory->grow(flux_neigh_intercepts[alloc], neigh_max_all+EXPAND, 2, "Pair CAC:flux_neigh_intercepts[alloc]");
      }
      local_all_max=neigh_max_all+EXPAND;
    }
  }

}

/* ---------------------------------------------------------------------- */

void PairCAC::init_quad_arrays(){

  int neigh_max_inner = inner_quad_lists_counts[pqi];
  int neigh_max_outer,neigh_max_add;
  int element_index, poly_index;
  int **node_types = atom->node_types;
  double **node_charges = atom->node_charges;
  int **inner_quad_indices = inner_quad_lists_index[pqi];
  int **outer_quad_indices, **add_quad_indices;

  if(outer_neighflag) {
    neigh_max_outer = outer_quad_lists_counts[pqi];
    outer_quad_indices = outer_quad_lists_index[pqi];
  }
  if(quad_flux_flag) {
    neigh_max_add = add_quad_lists_counts[pqi];
    add_quad_indices = add_quad_lists_index[pqi];
  }
  
  for (int l = 0; l < neigh_max_inner; l++) {
    element_index = inner_quad_indices[l][0];
    poly_index = inner_quad_indices[l][1];
    if(type_map)
    inner_neighbor_types[l] = type_maps[node_types[element_index][poly_index]];
    else
    inner_neighbor_types[l] = node_types[element_index][poly_index];
    if(atom->q_flag)
    inner_neighbor_charges[l] = node_charges[element_index][poly_index];
  }

  if(outer_neighflag)
    for (int l = 0; l < neigh_max_outer; l++) {
      element_index = outer_quad_indices[l][0];
      poly_index = outer_quad_indices[l][1];
      if(type_map)
        outer_neighbor_types[l] = type_maps[node_types[element_index][poly_index]];
      else
        outer_neighbor_types[l] = node_types[element_index][poly_index];
      if(atom->q_flag)
        outer_neighbor_charges[l] = node_charges[element_index][poly_index];
    }
  if(quad_flux_flag)
    for (int l = 0; l < neigh_max_add; l++) {
      element_index = add_quad_indices[l][0];
      poly_index = add_quad_indices[l][1];
      if(type_map)
        add_neighbor_types[l] = type_maps[node_types[element_index][poly_index]];
      else
        add_neighbor_types[l] = node_types[element_index][poly_index];
      if(atom->q_flag)
        add_neighbor_charges[l] = node_charges[element_index][poly_index];
    }
}

/* ---------------------------------------------------------------------- 
 Compute the current positions, velocities, etc. of all virtual/real atom
 neighbors for a quadrature point; uses the element shape functions and
 the parameters associated with each virtual atom neighbor position to do so.
 These parameters are stored in the ucell array.
---------------------------------------------------------------------- */

void PairCAC::interpolation(int iii, double sq, double tq, double wq){

  double s, t, w, sp, sm, tp, tm, wp, wm;
  int *element_type = atom->element_type;
  double ****nodal_positions = atom->nodal_positions;
  double ****nodal_velocities = atom->nodal_velocities;
  double **x = atom->x;
  double **v = atom->v;
  int neigh_max = inner_quad_lists_counts[pqi];
  int outer_neigh_max, add_neigh_max;
  double **inner_ucells, **outer_ucells, **add_ucells, shape_func;
  int inner_etypes[MAXNEIGHBUFF], outer_etypes[MAXNEIGHBUFF], add_etypes[MAXNEIGHBUFF];
  int **inner_indices, **outer_indices, **add_indices;

  if(neigh_max > MAXNEIGHBUFF)
    error->one(FLERR,"increase MAXNEIGHBUFF in pair_cac.cpp; this overestimates the maximum possible neighbors for a quadrature point");

  if(outer_neighflag){
    outer_neigh_max = outer_quad_lists_counts[pqi];
    if(outer_neigh_max > MAXNEIGHBUFF)
      error->one(FLERR,"increase MAXNEIGHBUFF in pair_cac.cpp; this overestimates the maximum possible neighbors for a quadrature point");	  
  }

  if(quad_flux_flag){
    add_neigh_max = add_quad_lists_counts[pqi];
    if(add_neigh_max > MAXNEIGHBUFF)
      error->one(FLERR,"increase MAXNEIGHBUFF in pair_cac.cpp; this overestimates the maximum possible neighbors for a quadrature point");	
  }

  double shaperesult[MAXNEIGHBUFF][MAXESHAPE];
  double **nodes, **node_velocities;
  inner_ucells = inner_quad_lists_ucell[pqi];
  inner_indices = inner_quad_lists_index[pqi];
  //init coords to zero etc.
  for (int l = 0; l < neigh_max; l++) {
    inner_neighbor_coords[l][0] = 0;
    inner_neighbor_coords[l][1] = 0;
    inner_neighbor_coords[l][2] = 0;
    if(quad_flux_flag){
      inner_neighbor_velocities[l][0] = 0;
      inner_neighbor_velocities[l][1] = 0;
      inner_neighbor_velocities[l][2] = 0;
    }
  }

  //compute interpolation for current quadrature point
  current_position[0]=0;
  current_position[1]=0;
  current_position[2]=0;
  if(quad_flux_flag){
    current_velocity[0]=0;
    current_velocity[1]=0;
    current_velocity[2]=0;
  }
  if(current_element_type == 1){
    s = sq;
    t = tq;
    w = wq;
    sp = (1+s)/8;
    sm = (1-s)/8;
    tp = 1+t;
    tm = 1-t;
    wp = 1+w;
    wm = 1-w;
    nodes = current_nodal_positions;
    if(quad_flux_flag)
    node_velocities = nodal_velocities[iii][poly_counter];
    shaperesult[0][0]= sm*tm*wm;
    shaperesult[0][1]= sp*tm*wm;
    shaperesult[0][2]= sp*tp*wm;
    shaperesult[0][3]= sm*tp*wm;
    shaperesult[0][4]= sm*tm*wp;
    shaperesult[0][5]= sp*tm*wp;
    shaperesult[0][6]= sp*tp*wp;
    shaperesult[0][7]= sm*tp*wp;
    for(int r = 0; r < 8; r++){
      current_position[0]+=nodes[r][0]*shaperesult[0][r];
      current_position[1]+=nodes[r][1]*shaperesult[0][r];
      current_position[2]+=nodes[r][2]*shaperesult[0][r];
      if(quad_flux_flag){
        current_velocity[0]+=node_velocities[r][0]*shaperesult[0][r];
        current_velocity[1]+=node_velocities[r][1]*shaperesult[0][r];
        current_velocity[2]+=node_velocities[r][2]*shaperesult[0][r];
      }
    }
  }
  else if (current_element_type == 2){
   //fill with other element type shape function
  }
  else 
  {
    current_position[0] = sq;
    current_position[1] = tq;
    current_position[2] = wq;
    if(quad_flux_flag){
      current_velocity[0] = v[iii][0];
      current_velocity[1] = v[iii][1];
      current_velocity[2] = v[iii][2];
    }
  }

  //compute interpolation for neighbor coordinates
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
    if(quad_flux_flag)
    node_velocities = nodal_velocities[inner_indices[l][0]][inner_indices[l][1]];
    shaperesult[l][0]= sm*tm*wm;
    shaperesult[l][1]= sp*tm*wm;
    shaperesult[l][2]= sp*tp*wm;
    shaperesult[l][3]= sm*tp*wm;
    shaperesult[l][4]= sm*tm*wp;
    shaperesult[l][5]= sp*tm*wp;
    shaperesult[l][6]= sp*tp*wp;
    shaperesult[l][7]= sm*tp*wp;
    for(int r = 0; r < 8; r++){
      inner_neighbor_coords[l][0]+=nodes[r][0]*shaperesult[l][r];
      inner_neighbor_coords[l][1]+=nodes[r][1]*shaperesult[l][r];
      inner_neighbor_coords[l][2]+=nodes[r][2]*shaperesult[l][r];
      if(quad_flux_flag){
        inner_neighbor_velocities[l][0]+=node_velocities[r][0]*shaperesult[l][r];
        inner_neighbor_velocities[l][1]+=node_velocities[r][1]*shaperesult[l][r];
        inner_neighbor_velocities[l][2]+=node_velocities[r][2]*shaperesult[l][r];
      }
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
    if(quad_flux_flag){
      inner_neighbor_velocities[l][0] = v[inner_indices[l][0]][0];
      inner_neighbor_velocities[l][1] = v[inner_indices[l][0]][1];
      inner_neighbor_velocities[l][2] = v[inner_indices[l][0]][2];
    }
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
    if(quad_flux_flag){
      outer_neighbor_velocities[l][0] = 0;
      outer_neighbor_velocities[l][1] = 0;
      outer_neighbor_velocities[l][2] = 0;
    }
  }

  //compute interpolation for outer neighbor coordinates
  for (int l = 0; l < outer_neigh_max; l++) {
    outer_etypes[l] = element_type[outer_indices[l][0]];
    //Q8 interpolation scheme
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
    if(quad_flux_flag)
    node_velocities = nodal_velocities[outer_indices[l][0]][outer_indices[l][1]];
    shaperesult[l][0]= sm*tm*wm;
    shaperesult[l][1]= sp*tm*wm;
    shaperesult[l][2]= sp*tp*wm;
    shaperesult[l][3]= sm*tp*wm;
    shaperesult[l][4]= sm*tm*wp;
    shaperesult[l][5]= sp*tm*wp;
    shaperesult[l][6]= sp*tp*wp;
    shaperesult[l][7]= sm*tp*wp;
    for(int r = 0; r < 8; r++){
      outer_neighbor_coords[l][0]+=nodes[r][0]*shaperesult[l][r];
      outer_neighbor_coords[l][1]+=nodes[r][1]*shaperesult[l][r];
      outer_neighbor_coords[l][2]+=nodes[r][2]*shaperesult[l][r];
      if(quad_flux_flag){
        outer_neighbor_velocities[l][0]+=node_velocities[r][0]*shaperesult[l][r];
        outer_neighbor_velocities[l][1]+=node_velocities[r][1]*shaperesult[l][r];
        outer_neighbor_velocities[l][2]+=node_velocities[r][2]*shaperesult[l][r];
      }
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
    if(quad_flux_flag){
      outer_neighbor_velocities[l][0] = v[outer_indices[l][0]][0];
      outer_neighbor_velocities[l][1] = v[outer_indices[l][0]][1];
      outer_neighbor_velocities[l][2] = v[outer_indices[l][0]][2];
    }
    }
  }	  
  }

  if(quad_flux_flag){
  add_ucells = add_quad_lists_ucell[pqi];
  add_indices = add_quad_lists_index[pqi];
  //init coords to zero etc.
  for (int l = 0; l < add_neigh_max; l++) {
  add_neighbor_coords[l][0] = 0;
  add_neighbor_coords[l][1] = 0;
  add_neighbor_coords[l][2] = 0;
  if(quad_flux_flag){
    add_neighbor_velocities[l][0] = 0;
    add_neighbor_velocities[l][1] = 0;
    add_neighbor_velocities[l][2] = 0;
  }
  }
  
  //compute interpolation for additional neighbor coordinates
  for (int l = 0; l < add_neigh_max; l++) {
    add_etypes[l] = element_type[add_indices[l][0]];
    //Q8 interpolation scheme
    if(add_etypes[l]==1){
    s = add_ucells[l][0];
    t = add_ucells[l][1];
    w = add_ucells[l][2];
    sp = (1+s)/8;
    sm = (1-s)/8;
    tp = 1+t;
    tm = 1-t;
    wp = 1+w;
    wm = 1-w;
    nodes = nodal_positions[add_indices[l][0]][add_indices[l][1]];
    if(quad_flux_flag)
    node_velocities = nodal_velocities[add_indices[l][0]][add_indices[l][1]];
    shaperesult[l][0]= sm*tm*wm;
    shaperesult[l][1]= sp*tm*wm;
    shaperesult[l][2]= sp*tp*wm;
    shaperesult[l][3]= sm*tp*wm;
    shaperesult[l][4]= sm*tm*wp;
    shaperesult[l][5]= sp*tm*wp;
    shaperesult[l][6]= sp*tp*wp;
    shaperesult[l][7]= sm*tp*wp;
    for(int r = 0; r < 8; r++){
      add_neighbor_coords[l][0]+=nodes[r][0]*shaperesult[l][r];
      add_neighbor_coords[l][1]+=nodes[r][1]*shaperesult[l][r];
      add_neighbor_coords[l][2]+=nodes[r][2]*shaperesult[l][r];
      if(quad_flux_flag){
        add_neighbor_velocities[l][0]+=node_velocities[r][0]*shaperesult[l][r];
        add_neighbor_velocities[l][1]+=node_velocities[r][1]*shaperesult[l][r];
        add_neighbor_velocities[l][2]+=node_velocities[r][2]*shaperesult[l][r];
      }
    }
  }
    //add new shape function block here
    if (add_etypes[l]==2){

    }
    //consider atoms
    if(add_etypes[l]==0){
    add_neighbor_coords[l][0] = x[add_indices[l][0]][0];
    add_neighbor_coords[l][1] = x[add_indices[l][0]][1];
    add_neighbor_coords[l][2] = x[add_indices[l][0]][2];
    if(quad_flux_flag){
      add_neighbor_velocities[l][0] = v[add_indices[l][0]][0];
      add_neighbor_velocities[l][1] = v[add_indices[l][0]][1];
      add_neighbor_velocities[l][2] = v[add_indices[l][0]][2];
    }
    }
  }
  }
  
}

/* ---------------------------------------------------------------------- 
 Find which surfaces are intersected by any given pair of virtual atoms 
 in the virtual atom neighbor list.
---------------------------------------------------------------------- */

void PairCAC::compute_intersections(){
  int all_neigh, is, isl, normal_flag, sign, scan_type1, scan_type2, index, jindex;
  int dim1, dim2, intersection_flag, neigh_intersect;
  double m, fpair, interaction_forceij[3], interaction_forceji[3], delxa[3];
  double scan_position1[3], scan_position2[3], q1, q2, delx, dely, delz, distancesq;  
  double intersection_point[3], proj, lparam, planecoord, plane_limits[2][2];
  double interactionx, interactiony, interactionz;
  double *box_center = atom->box_center;
  double *box_size = atom->box_size;
  int neigh_max = inner_quad_lists_counts[pqi];
  int neigh_max_outer;
  int neigh_add = add_quad_lists_counts[pqi];
  neigh_max_outer = inner_quad_lists_counts[pqi];
  all_neigh = neigh_max + neigh_max_outer + neigh_add;

  //consider plane intersections for the line between the quadrature point and neighbor l
  //reinitialize arrays
  for(int ineigh=0; ineigh < all_neigh; ineigh++){
    flux_neigh_nintercepts[ineigh] = 0;
    for(int jneigh=0; jneigh < all_neigh; jneigh++)
      nplane_intersects[ineigh][jneigh] = 0;
  }
  //determine which of the 6 planes of the atom box are intersected by this i-j pair
  for(int ineigh=0; ineigh < all_neigh; ineigh++){
    if(ineigh<neigh_max){
      scan_type1 = inner_neighbor_types[ineigh];
      if(atom->q_flag) q1 = inner_neighbor_charges[ineigh];
      scan_position1[0] = inner_neighbor_coords[ineigh][0];
      scan_position1[1] = inner_neighbor_coords[ineigh][1];
      scan_position1[2] = inner_neighbor_coords[ineigh][2];
    }
    else if(ineigh>=neigh_max && ineigh < neigh_max + neigh_max_outer){
      scan_type1 = outer_neighbor_types[ineigh-neigh_max];
      if(atom->q_flag) q1 = outer_neighbor_charges[ineigh-neigh_max];
      scan_position1[0] = outer_neighbor_coords[ineigh-neigh_max][0];
      scan_position1[1] = outer_neighbor_coords[ineigh-neigh_max][1];
      scan_position1[2] = outer_neighbor_coords[ineigh-neigh_max][2];
    }
    else{
      index = ineigh-neigh_max-neigh_max_outer;
      scan_type1 = add_neighbor_types[index];
      if(atom->q_flag) q1 = add_neighbor_charges[index];
      scan_position1[0] = add_neighbor_coords[index][0];
      scan_position1[1] = add_neighbor_coords[index][1];
      scan_position1[2] = add_neighbor_coords[index][2];
    }
    for(int jneigh=ineigh+1; jneigh < all_neigh; jneigh++){
      neigh_intersect = 0;
      if(jneigh<neigh_max){
        scan_type2 = inner_neighbor_types[jneigh];
        if(atom->q_flag) q2 = inner_neighbor_charges[jneigh];
        scan_position2[0] = inner_neighbor_coords[jneigh][0];
        scan_position2[1] = inner_neighbor_coords[jneigh][1];
        scan_position2[2] = inner_neighbor_coords[jneigh][2];
      }
      else if(jneigh>=neigh_max && jneigh < neigh_max + neigh_max_outer){
        scan_type2 = outer_neighbor_types[jneigh-neigh_max];
        if(atom->q_flag) q2 = outer_neighbor_charges[jneigh-neigh_max];
        scan_position2[0] = outer_neighbor_coords[jneigh-neigh_max][0];
        scan_position2[1] = outer_neighbor_coords[jneigh-neigh_max][1];
        scan_position2[2] = outer_neighbor_coords[jneigh-neigh_max][2];
      }
      else{
        jindex = jneigh-neigh_max-neigh_max_outer;
        scan_type2 = add_neighbor_types[jindex];
        if(atom->q_flag) q2 = add_neighbor_charges[jindex];
        scan_position2[0] = add_neighbor_coords[jindex][0];
        scan_position2[1] = add_neighbor_coords[jindex][1];
        scan_position2[2] = add_neighbor_coords[jindex][2];
      }

      delxa[0] = delx = scan_position1[0] - scan_position2[0];
      delxa[1] = dely = scan_position1[1] - scan_position2[1];
      delxa[2] = delz = scan_position1[2] - scan_position2[2];
      distancesq = delx*delx + dely*dely + delz*delz;
      if(distancesq > cutsq[scan_type1][scan_type2]) continue;

      for(int isl=0; isl < 2*domain->dimension; isl++){
        is = isl/2;
        if(isl%2==0) normal_flag = 1;
        else normal_flag = -1;
        if(is==0){
          dim1 = 1;
          dim2 = 2;
        }
        if(is==1){
          dim1 = 0;
          dim2 = 2;
        }
        if(is==2){
          dim1 = 0;
          dim2 = 1;
        }

        //test negative and positive sides of the box dimension
        if(isl%2==0) planecoord = current_position[is]-box_size[is]/2 + box_center[is];
        else planecoord = current_position[is]+box_size[is]/2 + box_center[is];
        plane_limits[0][0] = current_position[dim1]-box_size[dim1]/2 + box_center[dim1]-PLANE_EPSILON;
        plane_limits[0][1] = current_position[dim1]+box_size[dim1]/2 + box_center[dim1]+PLANE_EPSILON;
        plane_limits[1][0] = current_position[dim2]-box_size[dim2]/2 + box_center[dim2]-PLANE_EPSILON;
        plane_limits[1][1] = current_position[dim2]+box_size[dim2]/2 + box_center[dim2]+PLANE_EPSILON;

        intersection_flag = 1;
        //compute perpendicular projection of the line connecting i and j
        proj = scan_position1[is]-planecoord;

        //test if i-j normal coordinates are on opposing sides of the plane
        if((proj<0&&scan_position2[is]<planecoord)||((proj>0)&&scan_position2[is]>planecoord)) intersection_flag = 0;

        //use the ratio between this projection and the i-j displacement normal to the plane
        //to define the line parameter (0-1) at the point of intersection
        lparam = proj/(delxa[is]);
        if(delxa[is]==0) intersection_flag = 0;

        //use line parameter to extrapolate the possible intersection point between i-j
        intersection_point[dim1] = scan_position2[dim1]+delxa[dim1]*(1-lparam);
        intersection_point[dim2] = scan_position2[dim2]+delxa[dim2]*(1-lparam);

        //test the tangential coordinates to determine if the line through i-j crosses the finite sized plane
        if(intersection_point[dim1]<=plane_limits[0][0]||intersection_point[dim1]>plane_limits[0][1]) intersection_flag = 0;
        if(intersection_point[dim2]<=plane_limits[1][0]||intersection_point[dim2]>plane_limits[1][1]) intersection_flag = 0;
        if(intersection_flag){
          neigh_intersect = 1;
          int plane_index = nplane_intersects[ineigh][flux_neigh_nintercepts[ineigh]];
          flux_neigh_intercepts[ineigh][flux_neigh_nintercepts[ineigh]][plane_index] = isl;
          nplane_intersects[ineigh][flux_neigh_nintercepts[ineigh]]++;
        }
        if(nplane_intersects[ineigh][flux_neigh_nintercepts[ineigh]]<2){
          flux_neigh_indices[ineigh][flux_neigh_nintercepts[ineigh]] = jneigh;
        }
        //can intersect with box at most twice
        if(nplane_intersects[ineigh][flux_neigh_nintercepts[ineigh]]==2)
          break;
  
      }
      if(neigh_intersect)
      flux_neigh_nintercepts[ineigh]++;
    }
  }
}

/* ----------------------------------------------------------------------
 Compute the cac flux density due to a quadrature points immediate 
 virtual/real neighbors
---------------------------------------------------------------------- */

void PairCAC::current_quad_flux(int jneigh, double interactionx, double interactiony, double interactionz){
  int is, isl, normal_flag, sign, scan_type1, scan_type2,index, dim1, dim2, intersection_flag;
  double vix, viy, viz, vjx, vjy, vjz, m, fpair, interaction_forceij[3], interaction_forceji[3];
  double scan_position[3], intersection_point[3], proj, lparam, planecoord, plane_limits[2][2];
  double delx, dely, delz, delxa[3], *xj;
  double *box_center = atom->box_center;
  double *box_size = atom->box_size;
  xj = inner_neighbor_coords[jneigh];
  delxa[0] = current_position[0] - xj[0];
  delxa[1] = current_position[1] - xj[1];
  delxa[2] = current_position[2] - xj[2];
  vix = current_velocity[0];
  viy = current_velocity[1];
  viz = current_velocity[2];
  vjx = inner_neighbor_velocities[jneigh][0];
  vjy = inner_neighbor_velocities[jneigh][1];
  vjz = inner_neighbor_velocities[jneigh][2];
  if(flux_enable==1){
  interaction_forceij[0] = -interactionx;
  interaction_forceij[1] = -interactiony;
  interaction_forceij[2] = -interactionz;
  }
  else{
  interaction_forceij[0] = interactionx;
  interaction_forceij[1] = interactiony;
  interaction_forceij[2] = interactionz;
  }
  for(int isl=0; isl < 2*domain->dimension; isl++){
    is = isl/2;
    if(isl%2==0) normal_flag = 1;
    else normal_flag = -1;
    if(is==0){
      dim1 = 1;
      dim2 = 2;
    }
    if(is==1){
      dim1 = 0;
      dim2 = 2;
    }
    if(is==2){
      dim1 = 0;
      dim2 = 1;
    }

    //test negative and positive sides of the box dimension
    if(isl%2==0) planecoord = current_position[is]-box_size[is]/2 + box_center[is];
    else planecoord = current_position[is]+box_size[is]/2 + box_center[is];
    plane_limits[0][0] = current_position[dim1]-box_size[dim1]/2 + box_center[dim1]-PLANE_EPSILON;
    plane_limits[0][1] = current_position[dim1]+box_size[dim1]/2 + box_center[dim1]+PLANE_EPSILON;
    plane_limits[1][0] = current_position[dim2]-box_size[dim2]/2 + box_center[dim2]-PLANE_EPSILON;
    plane_limits[1][1] = current_position[dim2]+box_size[dim2]/2 + box_center[dim2]+PLANE_EPSILON;

    intersection_flag = 1;
    //compute perpendicular projection of the line connecting i and j
    proj = current_position[is]-planecoord;

    //test if i-j normal coordinates are on opposing sides of the plane
    if((proj<0&&xj[is]<planecoord)||((proj>0)&&xj[is]>planecoord)) intersection_flag = 0;

    //use the ratio between this projection and the i-j displacement normal to the plane
    //to define the line parameter (0-1) at the point of intersection
    lparam = proj/(delxa[is]);
    if(delxa[is]==0) intersection_flag = 0;

    //use line parameter to extrapolate the possible intersection point between i-j
    intersection_point[dim1] = xj[dim1]+delxa[dim1]*(1-lparam);
    intersection_point[dim2] = xj[dim2]+delxa[dim2]*(1-lparam);

    //test the tangential coordinates to determine if the line through i-j crosses the finite sized plane
    if(intersection_point[dim1]<=plane_limits[0][0]||intersection_point[dim1]>plane_limits[0][1]) intersection_flag = 0;
    if(intersection_point[dim2]<=plane_limits[1][0]||intersection_point[dim2]>plane_limits[1][1]) intersection_flag = 0;
    
    if(intersection_flag){
      if(current_position[is]<planecoord)
      sign=-normal_flag;
      else
      sign=normal_flag;
      //flux_enable is 1 in the case of pair forces, 2 in the case of many-body
      if(flux_enable==1){
        if(isl==0){
            //flux_contrib[jneigh] = -interaction_forceij[0]*sign;
          }
        flux_density[4*isl] += (interaction_forceij[0]*(vix+vjx) + 
        interaction_forceij[1]*(viy+vjy)+interaction_forceij[2]*(viz+vjz))*sign;
        flux_density[4*isl+1] -= interaction_forceij[0]*sign;
        flux_density[4*isl+2] -= interaction_forceij[1]*sign;
        flux_density[4*isl+3] -= interaction_forceij[2]*sign;
      }
      else{
        flux_density[4*isl] += (interaction_forceij[0]*vix+interaction_forceij[1]*viy+interaction_forceij[2]*viz
          -interaction_forceji[0]*vjx-interaction_forceji[1]*vjy-interaction_forceji[2]*vjz)*sign;
        flux_density[4*isl+1] -= (interaction_forceij[0]-interaction_forceji[0])*sign;
        flux_density[4*isl+2] -= (interaction_forceij[1]-interaction_forceji[1])*sign;
        flux_density[4*isl+3] -= (interaction_forceij[2]-interaction_forceji[2])*sign;
      }
      break;
    }
  }
}

/* ---------------------------------------------------------------------- 
 Compute the cac flux density due to virtual neighbors around a quadrature point
---------------------------------------------------------------------- */

void PairCAC::quad_neigh_flux(){
  int all_neigh, is, isl, normal_flag, sign, scan_type1, scan_type2, index, jindex;
  int dim1, dim2, intersection_flag, intersection_count, icontrib;
  double m, fpair, interaction_forceij[3], interaction_forceji[3], delxa[3];
  double scan_position1[3], scan_position2[3], q1, q2, delx, dely, delz, distancesq, fluxdistsq;  
  double intersection_point[3], proj, lparam, planecoord, plane_limits[2][2];
  double vix, viy, viz, vjx, vjy, vjz, interactionx, interactiony, interactionz;
  double *box_center = atom->box_center;
  double *box_size = atom->box_size;
  int neigh_max = inner_quad_lists_counts[pqi];
  int neigh_max_outer;
  int neigh_add = add_quad_lists_counts[pqi];
  double cut_add = atom->cut_add;
  double fluxcutsq = (cut_global_s+cut_add)*(cut_global_s+cut_add);
  if(outer_neighflag)
    neigh_max_outer = inner_quad_lists_counts[pqi];
  if(outer_neighflag)
    all_neigh = neigh_max + neigh_max_outer + neigh_add;
  else
    all_neigh = neigh_max + neigh_add;
  //icontrib = 0;

  //determine which of the 6 planes of the atom box are intersected by a given i-j pair
  for(int ineigh=0; ineigh < all_neigh; ineigh++){
    if(ineigh<neigh_max){
      scan_type1 = inner_neighbor_types[ineigh];
      if(atom->q_flag) q1 = inner_neighbor_charges[ineigh];
      scan_position1[0] = inner_neighbor_coords[ineigh][0];
      scan_position1[1] = inner_neighbor_coords[ineigh][1];
      scan_position1[2] = inner_neighbor_coords[ineigh][2];
      vix = inner_neighbor_velocities[ineigh][0];
      viy = inner_neighbor_velocities[ineigh][1];
      viz = inner_neighbor_velocities[ineigh][2];
    }
    else if(ineigh>=neigh_max && ineigh < neigh_max + neigh_max_outer){
      scan_type1 = outer_neighbor_types[ineigh-neigh_max];
      if(atom->q_flag) q1 = outer_neighbor_charges[ineigh-neigh_max];
      scan_position1[0] = outer_neighbor_coords[ineigh-neigh_max][0];
      scan_position1[1] = outer_neighbor_coords[ineigh-neigh_max][1];
      scan_position1[2] = outer_neighbor_coords[ineigh-neigh_max][2];
      vix = outer_neighbor_velocities[ineigh-neigh_max][0];
      viy = outer_neighbor_velocities[ineigh-neigh_max][1];
      viz = outer_neighbor_velocities[ineigh-neigh_max][2];
    }
    else{
      if(outer_neighflag)
        index = ineigh-neigh_max-neigh_max_outer;
      else
        index = ineigh-neigh_max;
      scan_type1 = add_neighbor_types[index];
      if(atom->q_flag) q1 = add_neighbor_charges[index];
      scan_position1[0] = add_neighbor_coords[index][0];
      scan_position1[1] = add_neighbor_coords[index][1];
      scan_position1[2] = add_neighbor_coords[index][2];
      vix = add_neighbor_velocities[index][0];
      viy = add_neighbor_velocities[index][1];
      viz = add_neighbor_velocities[index][2];
    }

    delxa[0] = delx = current_position[0] - scan_position1[0];
    delxa[1] = dely = current_position[1] - scan_position1[1];
    delxa[2] = delz = current_position[2] - scan_position1[2];
    fluxdistsq = delx*delx + dely*dely + delz*delz;
    if(fluxdistsq > fluxcutsq) continue;

    for(int jneigh=ineigh+1; jneigh < all_neigh; jneigh++){
      intersection_count = 0;
      if(jneigh<neigh_max){
        scan_type2 = inner_neighbor_types[jneigh];
        if(atom->q_flag) q2 = inner_neighbor_charges[jneigh];
        scan_position2[0] = inner_neighbor_coords[jneigh][0];
        scan_position2[1] = inner_neighbor_coords[jneigh][1];
        scan_position2[2] = inner_neighbor_coords[jneigh][2];
        vjx = inner_neighbor_velocities[jneigh][0];
        vjy = inner_neighbor_velocities[jneigh][1];
        vjz = inner_neighbor_velocities[jneigh][2];
      }
      else if(jneigh>=neigh_max && jneigh < neigh_max + neigh_max_outer){
        jindex = jneigh - neigh_max;
        scan_type2 = outer_neighbor_types[jindex];
        if(atom->q_flag) q2 = outer_neighbor_charges[jindex];
        scan_position2[0] = outer_neighbor_coords[jindex][0];
        scan_position2[1] = outer_neighbor_coords[jindex][1];
        scan_position2[2] = outer_neighbor_coords[jindex][2];
        vjx = outer_neighbor_velocities[jindex][0];
        vjy = outer_neighbor_velocities[jindex][1];
        vjz = outer_neighbor_velocities[jindex][2];
      }
      else{
        if(outer_neighflag)
          jindex = jneigh-neigh_max-neigh_max_outer;
        else
          jindex = jneigh-neigh_max;
        scan_type2 = add_neighbor_types[jindex];
        if(atom->q_flag) q2 = add_neighbor_charges[jindex];
        scan_position2[0] = add_neighbor_coords[jindex][0];
        scan_position2[1] = add_neighbor_coords[jindex][1];
        scan_position2[2] = add_neighbor_coords[jindex][2];
        vjx = add_neighbor_velocities[jindex][0];
        vjy = add_neighbor_velocities[jindex][1];
        vjz = add_neighbor_velocities[jindex][2];
      }

      delxa[0] = delx = scan_position1[0] - scan_position2[0];
      delxa[1] = dely = scan_position1[1] - scan_position2[1];
      delxa[2] = delz = scan_position1[2] - scan_position2[2];
      distancesq = delx*delx + dely*dely + delz*delz;
      if(distancesq > cutsq[scan_type1][scan_type2]) continue;
      
      if(flux_enable==1){
          if(atom->q_flag) fpair = pair_interaction_q(distancesq, scan_type1, scan_type2, q1, q2);
          else fpair = pair_interaction(distancesq, scan_type1, scan_type2);
          interaction_forceij[0] = -delx*fpair;
          interaction_forceij[1] = -dely*fpair;
          interaction_forceij[2] = -delz*fpair;
      }
      for(int isl=0; isl < 2*domain->dimension; isl++){
        is = isl/2;
        if(is==0){
          dim1 = 1;
          dim2 = 2;
        }
        if(is==1){
          dim1 = 0;
          dim2 = 2;
        }
        if(is==2){
          dim1 = 0;
          dim2 = 1;
        }

        //test negative and positive sides of the box dimension
        if(isl%2==0) planecoord = current_position[is]-box_size[is]/2 + box_center[is];
        else planecoord = current_position[is]+box_size[is]/2 + box_center[is];
        plane_limits[0][0] = current_position[dim1]-box_size[dim1]/2 + box_center[dim1]-PLANE_EPSILON;
        plane_limits[0][1] = current_position[dim1]+box_size[dim1]/2 + box_center[dim1]+PLANE_EPSILON;
        plane_limits[1][0] = current_position[dim2]-box_size[dim2]/2 + box_center[dim2]-PLANE_EPSILON;
        plane_limits[1][1] = current_position[dim2]+box_size[dim2]/2 + box_center[dim2]+PLANE_EPSILON;

        intersection_flag = 1;
        //compute perpendicular projection of the line connecting i and j
        proj = scan_position1[is]-planecoord;

        //test if i-j normal coordinates are on opposing sides of the plane
        if((proj<0&&scan_position2[is]<planecoord)||((proj>0)&&scan_position2[is]>planecoord)) intersection_flag = 0;

        //use the ratio between this projection and the i-j displacement normal to the plane
        //to define the line parameter (0-1) at the point of intersection
        lparam = proj/(delxa[is]);
        if(delxa[is]==0) intersection_flag = 0;

        //use line parameter to extrapolate the possible intersection point between i-j
        intersection_point[dim1] = scan_position2[dim1]+delxa[dim1]*(1-lparam);
        intersection_point[dim2] = scan_position2[dim2]+delxa[dim2]*(1-lparam);

        //test the tangential coordinates to determine if the line through i-j crosses the finite sized plane
        if(intersection_point[dim1]<=plane_limits[0][0]||intersection_point[dim1]>plane_limits[0][1]) intersection_flag = 0;
        if(intersection_point[dim2]<=plane_limits[1][0]||intersection_point[dim2]>plane_limits[1][1]) intersection_flag = 0;
        if(intersection_flag){
          intersection_count++;
          if(isl%2==0) normal_flag = 1;
          else normal_flag = -1;

          if(scan_position1[is]<planecoord) sign=-normal_flag;
          else sign=normal_flag;
          //flux_enable is 1 in the case of pair forces, 2 in the case of many-body
        if(flux_enable==1){
          if(isl==0){
            //flux_contrib[icontrib][0] = -interaction_forceij[0]*sign;
            //flux_contrib[icontrib][1] = ineigh;
            //flux_contrib[icontrib][2] = jneigh;
            //flux_contrib[icontrib][3] = isl;
            //icontrib++;
          }
          flux_density[4*isl] += (interaction_forceij[0]*(vix+vjx) + 
          interaction_forceij[1]*(viy+vjy)+interaction_forceij[2]*(viz+vjz))*sign;
          flux_density[4*isl+1] -= interaction_forceij[0]*sign;
          flux_density[4*isl+2] -= interaction_forceij[1]*sign;
          flux_density[4*isl+3] -= interaction_forceij[2]*sign;
        }
        }

        //can intersect with box at most twice
        if(intersection_count==2)
          break;
  
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

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

/* ---------------------------------------------------------------------- */

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
