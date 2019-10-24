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
#define MAXNEIGH1  500
#define MAXNEIGH2  10
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
  cutoff_skin = 2;
  max_expansion_count_inner = 0;
  max_expansion_count_outer = 0;
  quadrature_point_max = 0;
  quadrature_point_data = NULL;
  quadrature_counts = NULL;
  interior_scales = NULL;
  surface_counts = NULL;
  old_atom_etype = NULL;
  old_all_atom_etype = NULL;
  quad_allocated = 0;
  surface_counts_max[0] = 1;
  surface_counts_max[1] = 1;
  surface_counts_max[2] = 1;
  surface_counts_max_old[0] = 1;
  surface_counts_max_old[1] = 1;
  surface_counts_max_old[2] = 1;
  local_inner_max = 0;
  local_outer_max = 0;
  outer_neighflag = 0;
  sector_flag = 0;
	ghost_quad = 0;
  one_layer_flag = 0;
  old_quad_minima= NULL;
  old_minima_neighbors= NULL;
  shape_quad_result = NULL;
  asa_pointer = NULL;
  densemax=0;
  neighbor->pgsize=10;
  neighbor->oneatom=1;
  old_atom_count=0;
  old_quad_count=0;
  atom->CAC_pair_flag=1;
  //allocate shape function pointers
  shape_functions= (Shape_Functions *) memory->smalloc(sizeof(Shape_Functions)*MAXESHAPE, "Pair CAC:shape_functions");
  set_shape_functions();

  //instance asa interface object
   asa_pointer = new Asa_Data(lmp, this);
}

/* ---------------------------------------------------------------------- */

PairCAC::~PairCAC() {


  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(mass_matrix);
    memory->destroy(surface_counts);
    memory->destroy(interior_scales);
    memory->destroy(quadrature_abcissae);
    memory->destroy(quadrature_weights);
    memory->destroy(force_column);
    memory->destroy(current_force_column);
    memory->destroy(current_nodal_forces);
    memory->destroy(pivot);
    memory->destroy(surf_set);
    memory->destroy(dof_set);
    memory->destroy(sort_surf_set);
    memory->destroy(sort_dof_set);
    memory->destroy(shape_quad_result);
	}

  if (quad_allocated) {
    for (int init = 0; init < old_atom_count; init++) {

      if (old_atom_etype[init] == 0) {

         memory->destroy(inner_quad_lists_ucell[init][0]);
         memory->destroy(inner_quad_lists_index[init][0]);
         memory->sfree(inner_quad_lists_ucell[init]);
         memory->sfree(inner_quad_lists_index[init]);
         memory->destroy(inner_quad_lists_counts[init]);
         if (outer_neighflag) {
           memory->destroy(outer_quad_lists_ucell[init][0]);
           memory->destroy(outer_quad_lists_index[init][0]);
           memory->sfree(outer_quad_lists_ucell[init]);
           memory->sfree(outer_quad_lists_index[init]);
           memory->destroy(outer_quad_lists_counts[init]);
         }
      }
      else {
        for (int neigh_loop = 0; neigh_loop < old_quad_count; neigh_loop++) {
           memory->destroy(inner_quad_lists_ucell[init][neigh_loop]);
           memory->destroy(inner_quad_lists_index[init][neigh_loop]);
         }
         memory->sfree(inner_quad_lists_ucell[init]);
         memory->sfree(inner_quad_lists_index[init]);
         memory->destroy(inner_quad_lists_counts[init]);
         if (outer_neighflag) {
           for (int neigh_loop = 0; neigh_loop < old_quad_count; neigh_loop++) {
               memory->destroy(outer_quad_lists_ucell[init][neigh_loop]);
               memory->destroy(outer_quad_lists_index[init][neigh_loop]);
             }
             memory->sfree(outer_quad_lists_ucell[init]);
             memory->sfree(outer_quad_lists_index[init]);
             memory->destroy(outer_quad_lists_counts[init]);
         }
      }
	
    }

    memory->sfree(inner_quad_lists_ucell);
    memory->sfree(inner_quad_lists_index);
    memory->sfree(inner_quad_lists_counts);
    memory->sfree(outer_quad_lists_ucell);
    memory->sfree(outer_quad_lists_index);
    memory->sfree(outer_quad_lists_counts);
    memory->destroy(old_quad_minima);
    memory->destroy(old_minima_neighbors);
  }
  delete asa_pointer;
  if(quadrature_point_max)
  memory->destroy(quadrature_point_data);
  if(nmax)
  memory->destroy(quadrature_counts);




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
	
  int *nodes_count_list = atom->nodes_per_element_list;	
  quad_eflag = eflag;
  cutoff_skin = neighbor->skin;
  
  reneighbor_time = neighbor->lastcall;
  if (update->ntimestep == reneighbor_time||update->whichflag==2){
    atom->neigh_weight_flag=1;
    atom->weight_count=nlocal;
  }
 
  // loop over neighbors of my atoms
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
  if (update->ntimestep == reneighbor_time||update->whichflag==2) {
  //count number of pure atoms in the local domain
  natomic = 0;
  for (i = 0; i < nlocal; i++) {
    current_element_type = element_type[i];
  if (current_element_type == 0) natomic += 1;
  }

  // initialize or grow surface counts array for quadrature scheme
  // along with interior scaling for the quadrature domain
  if(nlocal+atom->nghost > nmax_surf)
    allocate_surface_counts();
  if (nlocal  > nmax) {
    memory->grow(atom->neighbor_weights, nlocal,3, "Pair CAC:neighbor_weights");
	  memory->grow(quadrature_counts, nlocal, "Pair CAC:quadrature_counts");
    nmax = nlocal;
  }
			
  surface_counts_max_old[0] = surface_counts_max[0];
  surface_counts_max_old[1] = surface_counts_max[1];
  surface_counts_max_old[2] = surface_counts_max[2];
  quad_list_counter=0;
	if(ghost_quad) nsurface = atom->nlocal + atom->nghost;
	else nsurface = atom->nlocal;
  for (i = 0; i < nsurface; i++) {
				
  current_element_scale = element_scale[i];
  current_element_type = element_type[i];
  current_poly_count = poly_count[i];
  if (current_element_type == 0) atomic_counter += 1;
  if (current_element_type != 0) {

    for (poly_counter = 0; poly_counter < poly_count[i]; poly_counter++) {
    current_nodal_positions = nodal_positions[i][poly_counter];
    int poly_surface_count[3];
    compute_surface_depths(interior_scale[0], interior_scale[1], interior_scale[2],
    poly_surface_count[0], poly_surface_count[1], poly_surface_count[2], 1);

    if(poly_counter==0){
      surface_counts[i][0]=poly_surface_count[0];
      surface_counts[i][1]=poly_surface_count[1];
      surface_counts[i][2]=poly_surface_count[2];
      interior_scales[i][0]=interior_scale[0];
      interior_scales[i][1]=interior_scale[1];
      interior_scales[i][2]=interior_scale[2];
    }
    else{
      if (poly_surface_count[0] > surface_counts[i][0]){ 
        surface_counts[i][0] = poly_surface_count[0];
        interior_scales[i][0]=interior_scale[0]; }
      if (poly_surface_count[1] > surface_counts[i][1]){ 
        surface_counts[i][1] = poly_surface_count[1];
        interior_scales[i][1]=interior_scale[1]; }
      if (poly_surface_count[2] > surface_counts[i][2]){ 
        surface_counts[i][2] = poly_surface_count[2];
        interior_scales[i][2]=interior_scale[2]; }
      }
      }

      if (surface_counts[i][0] > surface_counts_max[0]) surface_counts_max[0] = surface_counts[i][0];
      if (surface_counts[i][1] > surface_counts_max[1]) surface_counts_max[1] = surface_counts[i][1];
      if (surface_counts[i][2] > surface_counts_max[2]) surface_counts_max[2] = surface_counts[i][2];
    }
  }

  // initialize or grow memory for the neighbor list of virtual atoms at quadrature points

  if (nlocal)
  allocate_quad_neigh_list(surface_counts_max[0], surface_counts_max[1], surface_counts_max[2]);
  }
  atomic_counter = 0;
  int **neighbor_weights = atom-> neighbor_weights;

  //quadrature point neighbor list construction
  if (update->ntimestep == reneighbor_time||update->whichflag==2){
  for (i = 0; i < nlocal; i++) {
    atomic_flag = 0;
    current_list_index = i;
    current_element_type = element_type[i];
    current_element_scale = element_scale[i];
    current_poly_count = poly_count[i];
    type_array = node_types[i];
    
    current_x = x[i];
    
    neighbor_weights[i][0]=0;
    neighbor_weights[i][1]=0;
    neighbor_weights[i][2]=0;
  
    //determine element type
    if (current_element_type == 0) atomic_counter += 1;
    if (current_element_type == 0) {
      atomic_flag = 1;
    }
    neigh_quad_counter = 0;
    //NOTE:might have to change matrices so they dont have zeros due to maximum node count; ill condition.
    if(atomic_flag){
      poly_counter = 0;
	    current_nodal_positions = nodal_positions[i][poly_counter];
      compute_quad_neighbors(i);
    }
    else{
      for (poly_counter = 0; poly_counter < current_poly_count; poly_counter++){
	      current_nodal_positions = nodal_positions[i][poly_counter]; 
        compute_quad_neighbors(i);
	  }
    }
	if(!atomic_flag)
	quadrature_counts[i]=neigh_quad_counter/current_poly_count;
	else 
	quadrature_counts[i]=1;
  }
  }
  
  quad_list_counter=0;
  //preforce calculations
  if(pre_force_flag)
  pre_force_densities();
 
  quad_list_counter=0;
  //compute forces
  for (i = 0; i < nlocal; i++) {
    atomic_flag = 0;
    current_list_index = i;
    current_element_type = element_type[i];
    current_element_scale = element_scale[i];
    current_poly_count = poly_count[i];
    type_array = node_types[i];
    current_x = x[i];
    if (eflag) {
      element_energy = 0;
    }
    //determine element type
    if (current_element_type == 0) atomic_counter += 1;
      nodes_per_element = nodes_count_list[current_element_type];
    if (current_element_type == 0) {
      atomic_flag = 1;
    }
    neigh_quad_counter = 0;
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

  if(update->whichflag==2)
  copy_vectors(1);
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
  memory->create(surf_set, 6, 2, "pairCAC:surf_set");
  memory->create(dof_set, 6, 4, "pairCAC:surf_set");
  memory->create(sort_surf_set, 6, 2, "pairCAC:surf_set");
  memory->create(sort_dof_set, 6, 4, "pairCAC:surf_set");
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
  maxneigh_quad_inner = MAXNEIGH2;
  maxneigh_quad_outer = MAXNEIGH1;
  int irequest = neighbor->request(this, instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->cac = 1;
  
  //surface selection array 
  surf_set[0][0] = 1;
  surf_set[0][1] = -1;
  surf_set[1][0] = 1;
  surf_set[1][1] = 1;
  surf_set[2][0] = 2;
  surf_set[2][1] = -1;
  surf_set[3][0] = 2;
  surf_set[3][1] = 1;
  surf_set[4][0] = 3;
  surf_set[4][1] = -1;
  surf_set[5][0] = 3;
  surf_set[5][1] = 1;

  //surface DOF array

  dof_set[0][0] = 0;
  dof_set[0][1] = 3;
  dof_set[0][2] = 4;
  dof_set[0][3] = 7;

  dof_set[1][0] = 1;
  dof_set[1][1] = 2;
  dof_set[1][2] = 5;
  dof_set[1][3] = 6;

  dof_set[2][0] = 0;
  dof_set[2][1] = 1;
  dof_set[2][2] = 4;
  dof_set[2][3] = 5;

  dof_set[3][0] = 2;
  dof_set[3][1] = 3;
  dof_set[3][2] = 6;
  dof_set[3][3] = 7;

  dof_set[4][0] = 0;
  dof_set[4][1] = 1;
  dof_set[4][2] = 2;
  dof_set[4][3] = 3;

  dof_set[5][0] = 4;
  dof_set[5][1] = 5;
  dof_set[5][2] = 6;
  dof_set[5][3] = 7;
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

void PairCAC::compute_quad_neighbors(int iii){
	double unit_cell_mapped[3];
	//int nodes_per_element;
	double coefficients;
	int signs, signt, signw;
	unit_cell_mapped[0] = 2 / double(current_element_scale[0]);
	unit_cell_mapped[1] = 2 / double(current_element_scale[1]);
	unit_cell_mapped[2] = 2 / double(current_element_scale[2]);

	double interior_scale[3];
	
	int isurface_counts[3];
	interior_scale[0] = interior_scales[iii][0];
	interior_scale[1] = interior_scales[iii][1];
	interior_scale[2] = interior_scales[iii][2];
	isurface_counts[0] = surface_counts[iii][0];
	isurface_counts[1] = surface_counts[iii][1];
	isurface_counts[2] = surface_counts[iii][2];

	double s, t, w;
	s = t = w = 0;
	double sq, tq, wq;
	//compute interior quadrature point virtual neighbor lists
	interior_flag = 1;
	if(quadrature_point_max==0) {
		quadrature_point_max = 1000*atom->maxpoly;
		memory->create(quadrature_point_data,quadrature_point_max,7,"pairCAC:quadrature_point_data");
	}
	if (!atomic_flag) {
		for (int i = 0; i < quadrature_node_count; i++) {
			for (int j = 0; j < quadrature_node_count; j++) {
				for (int k = 0; k < quadrature_node_count; k++) {
					sq=s = interior_scale[0] * quadrature_abcissae[i];
					tq=t = interior_scale[1] * quadrature_abcissae[j];
					wq=w = interior_scale[2] * quadrature_abcissae[k];
					signs=signt=signw=1;
					if(sq<0) signs=-1;
					if(tq<0) signt=-1;
					if(wq<0) signw=-1;
					s = unit_cell_mapped[0] * (int((s+signs) / unit_cell_mapped[0]))-signs;
					t = unit_cell_mapped[1] * (int((t+signt) / unit_cell_mapped[1]))-signt;
					w = unit_cell_mapped[2] * (int((w+signw) / unit_cell_mapped[2]))-signw;

					if (quadrature_abcissae[i] < 0)
						s = s - 0.5*unit_cell_mapped[0];
					else
						s = s + 0.5*unit_cell_mapped[0];

					if (quadrature_abcissae[j] < 0)
						t = t - 0.5*unit_cell_mapped[1];
					else
						t = t + 0.5*unit_cell_mapped[1];

					if (quadrature_abcissae[k] < 0)
						w = w - 0.5*unit_cell_mapped[2];
					else
						w = w + 0.5*unit_cell_mapped[2];

                    quadrature_point_data[quad_list_counter][0]=s;
					quadrature_point_data[quad_list_counter][1]=t;
					quadrature_point_data[quad_list_counter][2]=w;
					quadrature_point_data[quad_list_counter][3]=sq;
					quadrature_point_data[quad_list_counter][4]=tq;
					quadrature_point_data[quad_list_counter][5]=wq;
					quadrature_point_data[quad_list_counter][6]=
					  interior_scale[0] * interior_scale[1] * interior_scale[2] *
					  quadrature_weights[i] * quadrature_weights[j] * quadrature_weights[k];
					quad_list_build(iii, s, t, w);
					neigh_quad_counter = neigh_quad_counter + 1;
          			quad_list_counter+=1;
					if(quad_list_counter==quadrature_point_max)
					grow_quad_data();
				}
			}
		}



		//compute surface contributions to element

		int sign[2];
		sign[0] = -1;
		sign[1] = 1;
		interior_flag = 0;

		// s axis surface contributions
		for (int sc = 0; sc < 2; sc++) {
			for (int i = 0; i < isurface_counts[0]; i++) {
				for (int j = 0; j < quadrature_node_count; j++) {
					for (int k = 0; k < quadrature_node_count; k++) {
						s = sign[sc] - i*unit_cell_mapped[0] * sign[sc];
						sq = s = s - 0.5*unit_cell_mapped[0] * sign[sc];
						tq = t = interior_scale[1] * quadrature_abcissae[j];
						wq = w = interior_scale[2] * quadrature_abcissae[k];
						signs=signt=signw=1;
					    if(tq<0) signt=-1;
					    if(wq<0) signw=-1;
						t = unit_cell_mapped[1] * (int((t+signt) / unit_cell_mapped[1]))-signt;
						w = unit_cell_mapped[2] * (int((w+signw) / unit_cell_mapped[2]))-signw;

						if (quadrature_abcissae[k] < 0)
							w = w - 0.5*unit_cell_mapped[2];
						else
							w = w + 0.5*unit_cell_mapped[2];

						if (quadrature_abcissae[j] < 0)
							t = t - 0.5*unit_cell_mapped[1];
						else
							t = t + 0.5*unit_cell_mapped[1];

						quadrature_point_data[quad_list_counter][0]=s;
					    quadrature_point_data[quad_list_counter][1]=t;
					    quadrature_point_data[quad_list_counter][2]=w;
						quadrature_point_data[quad_list_counter][3]=sq;
					    quadrature_point_data[quad_list_counter][4]=tq;
					    quadrature_point_data[quad_list_counter][5]=wq;
						quadrature_point_data[quad_list_counter][6]=
					      unit_cell_mapped[0] * interior_scale[1] *
						  interior_scale[2] * quadrature_weights[j] * quadrature_weights[k];	
						quad_list_build(iii, s, t, w);
						neigh_quad_counter = neigh_quad_counter + 1;
           				quad_list_counter+=1;  
						if(quad_list_counter==quadrature_point_max)
					    grow_quad_data();      
					}
				}
			}
		}

		// t axis contributions

		for (int sc = 0; sc < 2; sc++) {
			for (int i = 0; i < isurface_counts[1]; i++) {
				for (int j = 0; j < quadrature_node_count; j++) {
					for (int k = 0; k < quadrature_node_count; k++) {
						sq = s = interior_scale[0] * quadrature_abcissae[j];
						t = sign[sc] - i*unit_cell_mapped[1] * sign[sc];
						tq=t = t - 0.5*unit_cell_mapped[1] * sign[sc];
						wq = w = interior_scale[2] * quadrature_abcissae[k];
            signs=signt=signw=1;
					  if(sq<0) signs=-1;
					  if(wq<0) signw=-1;
						s = unit_cell_mapped[0] * (int((s+signs) / unit_cell_mapped[0]))-signs;
						w = unit_cell_mapped[2] * (int((w+signw) / unit_cell_mapped[2]))-signw;

						if (quadrature_abcissae[j] < 0)
							s = s - 0.5*unit_cell_mapped[0];
						else
							s = s + 0.5*unit_cell_mapped[0];

						if (quadrature_abcissae[k] < 0)
							w = w - 0.5*unit_cell_mapped[2];
						else
							w = w + 0.5*unit_cell_mapped[2];
						
						quadrature_point_data[quad_list_counter][0]=s;
					    quadrature_point_data[quad_list_counter][1]=t;
					    quadrature_point_data[quad_list_counter][2]=w;
						quadrature_point_data[quad_list_counter][3]=sq;
					    quadrature_point_data[quad_list_counter][4]=tq;
					    quadrature_point_data[quad_list_counter][5]=wq;
						quadrature_point_data[quad_list_counter][6]=
					      unit_cell_mapped[1] * interior_scale[0] *
						  interior_scale[2] * quadrature_weights[j] * quadrature_weights[k];	
						quad_list_build(iii, s, t, w);
						neigh_quad_counter = neigh_quad_counter + 1;
            			quad_list_counter+=1;
						if(quad_list_counter==quadrature_point_max)
					    grow_quad_data(); 
					}
				}
			}
		}
        //w axis contributions
		for (int sc = 0; sc < 2; sc++) {
			for (int i = 0; i < isurface_counts[2]; i++) {
				for (int j = 0; j < quadrature_node_count; j++) {
					for (int k = 0; k < quadrature_node_count; k++) {

						sq = s = interior_scale[0] * quadrature_abcissae[j];
						tq = t = interior_scale[1] * quadrature_abcissae[k];
						signs=signt=signw=1;
					  if(sq<0) signs=-1;
					  if(tq<0) signt=-1;
						s = unit_cell_mapped[0] * (int((s+signs) / unit_cell_mapped[0]))-signs;
						t = unit_cell_mapped[1] * (int((t+signt) / unit_cell_mapped[1]))-signt;
						w = sign[sc] - i*unit_cell_mapped[2] * sign[sc];
						wq = w = w - 0.5*unit_cell_mapped[2] * sign[sc];

						if (quadrature_abcissae[j] < 0)
							s = s - 0.5*unit_cell_mapped[0];
						else
							s = s + 0.5*unit_cell_mapped[0];

						if (quadrature_abcissae[k] < 0)
							t = t - 0.5*unit_cell_mapped[1];
						else
							t = t + 0.5*unit_cell_mapped[1];

                        quadrature_point_data[quad_list_counter][0]=s;
					    quadrature_point_data[quad_list_counter][1]=t;
					    quadrature_point_data[quad_list_counter][2]=w;
						quadrature_point_data[quad_list_counter][3]=sq;
					    quadrature_point_data[quad_list_counter][4]=tq;
					    quadrature_point_data[quad_list_counter][5]=wq;
						quadrature_point_data[quad_list_counter][6]=
					      unit_cell_mapped[2] * interior_scale[1] *
						  interior_scale[0] * quadrature_weights[j] * quadrature_weights[k];	
						quad_list_build(iii, s, t, w);
						neigh_quad_counter = neigh_quad_counter + 1;
            			quad_list_counter+=1;
						if(quad_list_counter==quadrature_point_max)
					    grow_quad_data(); 
					}
				}
			}
		}

		int surface_countsx;
		int surface_countsy;
		double unit_mappedx;
		double unit_mappedy;
    double interior_scalez;
		//compute edge contributions

		for (int sc = 0; sc < 12; sc++) {
			if (sc == 0 || sc == 1 || sc == 2 || sc == 3) {
        interior_scalez = interior_scale[2];
				surface_countsx = isurface_counts[0];
				surface_countsy = isurface_counts[1];
				unit_mappedx = unit_cell_mapped[0];
				unit_mappedy = unit_cell_mapped[1];
			}
			else if (sc == 4 || sc == 5 || sc == 6 || sc == 7) {
        interior_scalez = interior_scale[0];
				surface_countsx = isurface_counts[1];
				surface_countsy = isurface_counts[2];
				unit_mappedx = unit_cell_mapped[1];
				unit_mappedy = unit_cell_mapped[2];
			}
			else if (sc == 8 || sc == 9 || sc == 10 || sc == 11) {
        interior_scalez = interior_scale[1];
				surface_countsx = isurface_counts[0];
				surface_countsy = isurface_counts[2];
				unit_mappedx = unit_cell_mapped[0];
				unit_mappedy = unit_cell_mapped[2];
			}
           
			for (int i = 0; i < surface_countsx; i++) {
				for (int j = 0; j < surface_countsy; j++) {
					for (int k = 0; k < quadrature_node_count; k++) {
						if (sc == 0) {

							sq = s = -1 + (i + 0.5)*unit_cell_mapped[0];
							tq = t = -1 + (j + 0.5)*unit_cell_mapped[1];
							wq = w = interior_scale[2] * quadrature_abcissae[k];
							signs=signt=signw=1;
					    if(wq<0) signw=-1;
							w = unit_cell_mapped[2] * (int((w+signw) / unit_cell_mapped[2]))-signw;
							if (quadrature_abcissae[k] < 0)
								w = w - 0.5*unit_cell_mapped[2];
							else
								w = w + 0.5*unit_cell_mapped[2];
						}
						else if (sc == 1) {
							sq = s = 1 - (i + 0.5)*unit_cell_mapped[0];
							tq = t = -1 + (j + 0.5)*unit_cell_mapped[1];
							wq = w = interior_scale[2] * quadrature_abcissae[k];
							signs=signt=signw=1;
					    if(wq<0) signw=-1;
							w = unit_cell_mapped[2] * (int((w+signw) / unit_cell_mapped[2]))-signw;
							if (quadrature_abcissae[k] < 0)
								w = w - 0.5*unit_cell_mapped[2];
							else
								w = w + 0.5*unit_cell_mapped[2];
						}
						else if (sc == 2) {
							sq = s = -1 + (i + 0.5)*unit_cell_mapped[0];
							tq = t = 1 - (j + 0.5)*unit_cell_mapped[1];
							wq = w = interior_scale[2] * quadrature_abcissae[k];
							signs=signt=signw=1;
					    if(wq<0) signw=-1;
							w = unit_cell_mapped[2] * (int((w+signw) / unit_cell_mapped[2]))-signw;
							if (quadrature_abcissae[k] < 0)
								w = w - 0.5*unit_cell_mapped[2];
							else
								w = w + 0.5*unit_cell_mapped[2];
						}
						else if (sc == 3) {
							sq = s = 1 - (i + 0.5)*unit_cell_mapped[0];
							tq = t = 1 - (j + 0.5)*unit_cell_mapped[1];
							wq = w = interior_scale[2] * quadrature_abcissae[k];
							signs=signt=signw=1;
					    if(wq<0) signw=-1;
							w = unit_cell_mapped[2] * (int((w+signw) / unit_cell_mapped[2]))-signw;
							if (quadrature_abcissae[k] < 0)
								w = w - 0.5*unit_cell_mapped[2];
							else
								w = w + 0.5*unit_cell_mapped[2];
						}
						else if (sc == 4) {
							sq = s = interior_scale[0] * quadrature_abcissae[k];
							signs=signt=signw=1;
					    if(sq<0) signs=-1;
							s = unit_cell_mapped[0] * (int((s+signs) / unit_cell_mapped[0]))-signs;
							tq = t = -1 + (i + 0.5)*unit_cell_mapped[1];
							wq = w = -1 + (j + 0.5)*unit_cell_mapped[2];
							if (quadrature_abcissae[k] < 0)
								s = s - 0.5*unit_cell_mapped[0];
							else
								s = s + 0.5*unit_cell_mapped[0];

						}
						else if (sc == 5) {
							sq = s = interior_scale[0] * quadrature_abcissae[k];
							signs=signt=signw=1;
					    if(sq<0) signs=-1;
							s = unit_cell_mapped[0] * (int((s+signs) / unit_cell_mapped[0]))-signs;
							tq = t = 1 - (i + 0.5)*unit_cell_mapped[1];
							wq = w = -1 + (j + 0.5)*unit_cell_mapped[2];
							if (quadrature_abcissae[k] < 0)
								s = s - 0.5*unit_cell_mapped[0];
							else
								s = s + 0.5*unit_cell_mapped[0];
						}
						else if (sc == 6) {
							sq = s = interior_scale[0] * quadrature_abcissae[k];
							signs=signt=signw=1;
					    if(sq<0) signs=-1;
							s = unit_cell_mapped[0] * (int((s+signs) / unit_cell_mapped[0]))-signs;
							tq = t = -1 + (i + 0.5)*unit_cell_mapped[1];
							wq = w = 1 - (j + 0.5)*unit_cell_mapped[2];
							if (quadrature_abcissae[k] < 0)
								s = s - 0.5*unit_cell_mapped[0];
							else
								s = s + 0.5*unit_cell_mapped[0];
						}
						else if (sc == 7) {
							sq = s = interior_scale[0] * quadrature_abcissae[k];
							signs=signt=signw=1;
					    if(sq<0) signs=-1;
							s = unit_cell_mapped[0] * (int((s+signs) / unit_cell_mapped[0]))-signs;
							tq = t = 1 - (i + 0.5)*unit_cell_mapped[1];
							wq = w = 1 - (j + 0.5)*unit_cell_mapped[2];
							if (quadrature_abcissae[k] < 0)
								s = s - 0.5*unit_cell_mapped[0];
							else
								s = s + 0.5*unit_cell_mapped[0];
						}
						else if (sc == 8) {
							sq = s = -1 + (i + 0.5)*unit_cell_mapped[0];
							tq = t = interior_scale[1] * quadrature_abcissae[k];
							signs=signt=signw=1;
					    if(tq<0) signt=-1;
							t = unit_cell_mapped[1] * (int((t+signt) / unit_cell_mapped[1]))-signt;
							wq = w = -1 + (j + 0.5)*unit_cell_mapped[2];
							if (quadrature_abcissae[k] < 0)
								t = t - 0.5*unit_cell_mapped[1];
							else
								t = t + 0.5*unit_cell_mapped[1];

						}
						else if (sc == 9) {
							sq = s = 1 - (i + 0.5)*unit_cell_mapped[0];
							tq = t = interior_scale[1] * quadrature_abcissae[k];
							signs=signt=signw=1;
					    if(tq<0) signt=-1;
							t = unit_cell_mapped[1] * (int((t+signt) / unit_cell_mapped[1]))-signt;
							wq = w = -1 + (j + 0.5)*unit_cell_mapped[2];
							if (quadrature_abcissae[k] < 0)
								t = t - 0.5*unit_cell_mapped[1];
							else
								t = t + 0.5*unit_cell_mapped[1];
						}
						else if (sc == 10) {
							sq = s = -1 + (i + 0.5)*unit_cell_mapped[0];
							tq = t = interior_scale[1] * quadrature_abcissae[k];
							signs=signt=signw=1;
					    if(tq<0) signt=-1;
							t = unit_cell_mapped[1] * (int((t+signt) / unit_cell_mapped[1]))-signt;
							wq = w = 1 - (j + 0.5)*unit_cell_mapped[2];
							if (quadrature_abcissae[k] < 0)
								t = t - 0.5*unit_cell_mapped[1];
							else
								t = t + 0.5*unit_cell_mapped[1];
						}
						else if (sc == 11) {
							sq = s = 1 - (i + 0.5)*unit_cell_mapped[0];
							tq = t = interior_scale[1] * quadrature_abcissae[k];
							signs=signt=signw=1;
					    if(tq<0) signt=-1;
							t = unit_cell_mapped[1] * (int((t+signt) / unit_cell_mapped[1]))-signt;
							wq = w = 1 - (j + 0.5)*unit_cell_mapped[2];
							if (quadrature_abcissae[k] < 0)
								t = t - 0.5*unit_cell_mapped[1];
							else
								t = t + 0.5*unit_cell_mapped[1];
						}
                        
						quadrature_point_data[quad_list_counter][0]=s;
					  quadrature_point_data[quad_list_counter][1]=t;
					  quadrature_point_data[quad_list_counter][2]=w;
						quadrature_point_data[quad_list_counter][3]=sq;
					  quadrature_point_data[quad_list_counter][4]=tq;
					  quadrature_point_data[quad_list_counter][5]=wq;
						quadrature_point_data[quad_list_counter][6]=
					  coefficients = unit_mappedx * unit_mappedy * interior_scalez *
						quadrature_weights[k];	
						quad_list_build(iii, s, t, w);
						neigh_quad_counter = neigh_quad_counter + 1;
						quad_list_counter+=1;
						if(quad_list_counter==quadrature_point_max)
					    grow_quad_data(); 
					}

				}
			}
		}


		//compute corner contributions
		for (int sc = 0; sc < 8; sc++) {
			for (int i = 0; i < isurface_counts[0]; i++) {
				for (int j = 0; j < isurface_counts[1]; j++) {
					for (int k = 0; k < isurface_counts[2]; k++) {
					
						if (sc == 0) {
							s = -1 + (i + 0.5)*unit_cell_mapped[0];
							t = -1 + (j + 0.5)*unit_cell_mapped[1];
							w = -1 + (k + 0.5)*unit_cell_mapped[2];
						}
						else if (sc == 1) {
							s = 1 - (i + 0.5)*unit_cell_mapped[0];
							t = -1 + (j + 0.5)*unit_cell_mapped[1];
							w = -1 + (k + 0.5)*unit_cell_mapped[2];
						}
						else if (sc == 2) {
							s = 1 - (i + 0.5)*unit_cell_mapped[0];
							t = 1 - (j + 0.5)*unit_cell_mapped[1];
							w = -1 + (k + 0.5)*unit_cell_mapped[2];
						}
						else if (sc == 3) {
							s = -1 + (i + 0.5)*unit_cell_mapped[0];
							t = 1 - (j + 0.5)*unit_cell_mapped[1];
							w = -1 + (k + 0.5)*unit_cell_mapped[2];
						}
						else if (sc == 4) {
							s = -1 + (i + 0.5)*unit_cell_mapped[0];
							t = -1 + (j + 0.5)*unit_cell_mapped[1];
							w = 1 - (k + 0.5)*unit_cell_mapped[2];

						}
						else if (sc == 5) {
							s = 1 - (i + 0.5)*unit_cell_mapped[0];
							t = -1 + (j + 0.5)*unit_cell_mapped[1];
							w = 1 - (k + 0.5)*unit_cell_mapped[2];
						}
						else if (sc == 6) {
							s = 1 - (i + 0.5)*unit_cell_mapped[0];
							t = 1 - (j + 0.5)*unit_cell_mapped[1];
							w = 1 - (k + 0.5)*unit_cell_mapped[2];
						}
						else if (sc == 7) {
							s = -1 + (i + 0.5)*unit_cell_mapped[0];
							t = 1 - (j + 0.5)*unit_cell_mapped[1];
							w = 1 - (k + 0.5)*unit_cell_mapped[2];
						}
                        
            quadrature_point_data[quad_list_counter][0]=s;
					  quadrature_point_data[quad_list_counter][1]=t;
					  quadrature_point_data[quad_list_counter][2]=w;
						quadrature_point_data[quad_list_counter][3]=s;
					  quadrature_point_data[quad_list_counter][4]=t;
					  quadrature_point_data[quad_list_counter][5]=w;
						quadrature_point_data[quad_list_counter][6]=
					  unit_cell_mapped[0] * unit_cell_mapped[1] * unit_cell_mapped[2];	
						quad_list_build(iii, s, t, w);
						neigh_quad_counter = neigh_quad_counter + 1;
						quad_list_counter+=1;
						if(quad_list_counter==quadrature_point_max)
					    grow_quad_data(); 
					}

				}
			}
		}
	}
	else {
		quadrature_point_data[quad_list_counter][0]=current_x[0];
		quadrature_point_data[quad_list_counter][1]=current_x[1];
		quadrature_point_data[quad_list_counter][2]=current_x[2];
		quadrature_point_data[quad_list_counter][3]=current_x[0];
		quadrature_point_data[quad_list_counter][4]=current_x[1];
		quadrature_point_data[quad_list_counter][5]=current_x[2];
		quadrature_point_data[quad_list_counter][6]= 1;
		
		quad_list_build(iii, current_x[0], current_x[1], current_x[2]);
		quad_list_counter+=1;
		if(quad_list_counter==quadrature_point_max)
		grow_quad_data(); 
		neigh_quad_counter = 1;
	}
}

//----------------------------------------------------------------

void PairCAC::compute_forcev(int iii){
	double unit_cell_mapped[3];
	int *nodes_count_list = atom->nodes_per_element_list;	
	double coefficients;
	int nodes_per_element;
	int init_quad_list_counter = quad_list_counter;
	double ****nodal_virial = atom->nodal_virial;
	double shape_func;
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
	  s = quadrature_point_data[init_quad_list_counter+quad_loop][0];
	  t = quadrature_point_data[init_quad_list_counter+quad_loop][1];
	  w = quadrature_point_data[init_quad_list_counter+quad_loop][2];
	  sq = quadrature_point_data[init_quad_list_counter+quad_loop][3];
	  tq = quadrature_point_data[init_quad_list_counter+quad_loop][4];
	  wq = quadrature_point_data[init_quad_list_counter+quad_loop][5];
	}
	coefficients = quadrature_point_data[init_quad_list_counter+quad_loop][6];
	if(!atomic_flag)
	force_densities(iii, s, t, w, coefficients,
	  force_density[0], force_density[1], force_density[2]);
	else
	force_densities(iii, current_x[0], current_x[1], current_x[2], coefficients,
	  force_density[0], force_density[1], force_density[2]);
	  neigh_quad_counter = neigh_quad_counter + 1;
	if(!atomic_flag){
	for (int js = 0; js < nodes_per_element; js++) {
		shape_func = shape_function(sq, tq, wq, 2, js + 1);
	  for (int jj = 0; jj < 3; jj++) {
		force_column[js][jj] += coefficients*force_density[jj] * shape_func;
	  }
		if(atom->CAC_virial){
		for (int jj = 0; jj < 6; jj++) {
			
		nodal_virial[iii][poly_counter][js][jj] += coefficients*virial_density[jj] * shape_func;
		}
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
	quad_list_counter++;
  }
}

//--------------------------------------------------------------------

void PairCAC::quad_list_build(int iii, double s, double t, double w) {

  double delx, dely, delz;
	double shape_func;
	double shape_func2;
	double boxmap_matrix[3][3];
	int neighborflag = 0;
	int outofbounds = 0;
	double unit_cell_mapped[3];
  double ****nodal_positions = atom->nodal_positions;
	double scanning_unit_cell[3];
	double unit_cell[3];
	double distancesq;
	double current_position[3];
	double scan_position[3];
	double rcut;
  int **neighbor_weights = atom-> neighbor_weights;

	int inner_neigh_index = 0;
	int outer_neigh_index = 0;
	int nodes_per_element;
	int *nodes_count_list = atom->nodes_per_element_list;	
	expansion_count_inner = 0;
	expansion_count_outer = 0;
	neighbor_element_positions = nodal_positions[iii];
	if (!atomic_flag) {
		//equivalent isoparametric cutoff range for a cube of rcut

		unit_cell_mapped[0] = 2 / double(current_element_scale[0]);
		unit_cell_mapped[1] = 2 / double(current_element_scale[1]);
		unit_cell_mapped[2] = 2 / double(current_element_scale[2]);

		unit_cell[0] = s;
		unit_cell[1] = t;
		unit_cell[2] = w;

		//scan the surrounding unit cell locations in a cartesian grid
		//of isoparametric space until the cutoff is exceeded
		//for each grid scan

		scanning_unit_cell[0] = unit_cell[0];
		scanning_unit_cell[1] = unit_cell[1];
		scanning_unit_cell[2] = unit_cell[2];

		current_position[0] = 0;
		current_position[1] = 0;
		current_position[2] = 0;
	
    nodes_per_element=nodes_count_list[current_element_type];
		
		for (int kkk = 0; kkk < nodes_per_element; kkk++) {
			shape_func = shape_function(unit_cell[0], unit_cell[1], unit_cell[2], 2, kkk + 1);
			current_position[0] += current_nodal_positions[kkk][0] * shape_func;
			current_position[1] += current_nodal_positions[kkk][1] * shape_func;
			current_position[2] += current_nodal_positions[kkk][2] * shape_func;
		}

		if (outer_neighflag) { rcut = 2 * cut_global_s + cutoff_skin; }
		else { rcut = cut_global_s + cutoff_skin; }

		//try making a boxmap matrix for every type later
		for (int id = 0; id < 3; id++) {
			for (int jd = 0; jd < 3; jd++) {
				boxmap_matrix[id][jd] = 0;
				for (int n = 0; n < nodes_per_element; n++) {
					boxmap_matrix[id][jd] += current_nodal_positions[n][id] * shape_function_derivative
					(s, t, w, 2, n + 1, jd + 1);
				}
			}
		}

		// initialize local lattice vector approximation
		double a[3];
		a[0] = unit_cell_mapped[0] * boxmap_matrix[0][0];
		a[1] = unit_cell_mapped[0] * boxmap_matrix[1][0];
		a[2] = unit_cell_mapped[0] * boxmap_matrix[2][0];

		double b[3];
		b[0] = unit_cell_mapped[1] * boxmap_matrix[0][1];
		b[1] = unit_cell_mapped[1] * boxmap_matrix[1][1];
		b[2] = unit_cell_mapped[1] * boxmap_matrix[2][1];

		double c[3];
		c[0] = unit_cell_mapped[2] * boxmap_matrix[0][2];
		c[1] = unit_cell_mapped[2] * boxmap_matrix[1][2];
		c[2] = unit_cell_mapped[2] * boxmap_matrix[2][2];

		//perform gram schmidt orthogonalization of the three vectors
		double norm_a = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
		double norm_b = sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);

		double proj_b2a = (b[0] * a[0] + b[1] * a[1] + b[2] * a[2]) / norm_a;

		//compute component of b normal to a
		double b_orth[3];
		b_orth[0] = b[0] - proj_b2a*a[0] / norm_a;
		b_orth[1] = b[1] - proj_b2a*a[1] / norm_a;
		b_orth[2] = b[2] - proj_b2a*a[2] / norm_a;

		double proj_c2a = (c[0] * a[0] + c[1] * a[1] + c[2] * a[2]) / norm_a;
		double norm_b_orth = sqrt(b_orth[0] * b_orth[0] + b_orth[1] * b_orth[1]
			+ b_orth[2] * b_orth[2]);
		double proj_c2b_orth = (b_orth[0] * c[0] + b_orth[1] * c[1] + b_orth[2] * c[2]) / norm_b_orth;

		//compute component of c normal to a and b
		double c_orth[3];
		c_orth[0] = c[0] - proj_c2a*a[0] / norm_a - proj_c2b_orth*b_orth[0] / norm_b_orth;
		c_orth[1] = c[1] - proj_c2a*a[1] / norm_a - proj_c2b_orth*b_orth[1] / norm_b_orth;
		c_orth[2] = c[2] - proj_c2a*a[2] / norm_a - proj_c2b_orth*b_orth[2] / norm_b_orth;
		double norm_c_orth = sqrt(c_orth[0] * c_orth[0] + c_orth[1] * c_orth[1]
			+ c_orth[2] * c_orth[2]);

		int w_span = int(rcut / norm_c_orth) + 1;

		int t_upper_limit;
		int t_lower_limit;
		int s_upper_limit;
		int s_lower_limit;
		for (int polyscan = 0; polyscan < current_poly_count; polyscan++) {
			for (int wcount = -w_span; wcount < w_span + 1; wcount++) {
				t_lower_limit = -int((rcut + proj_c2b_orth*wcount) / norm_b_orth) - 1;
				t_upper_limit = int((rcut - proj_c2b_orth*wcount) / norm_b_orth) + 1;

				for (int tcount = t_lower_limit; tcount < t_upper_limit + 1; tcount++) {
					s_lower_limit = -int((rcut + proj_c2a*wcount + proj_b2a*tcount) / norm_a) - 1;
					s_upper_limit = int((rcut - proj_c2a*wcount - proj_b2a*tcount) / norm_a) + 1;

					for (int scount = s_lower_limit; scount < s_upper_limit + 1; scount++) {
						//scanning around atom
						outofbounds = 0;
						if (((scount == 0 && tcount == 0) && wcount == 0) && polyscan == poly_counter) { continue; }

						scanning_unit_cell[0] = scount*unit_cell_mapped[0] + unit_cell[0];
						scanning_unit_cell[1] = tcount*unit_cell_mapped[1] + unit_cell[1];
						scanning_unit_cell[2] = wcount*unit_cell_mapped[2] + unit_cell[2];

						scan_position[0] = 0;
						scan_position[1] = 0;
						scan_position[2] = 0;

						if (scanning_unit_cell[0] < -1 || scanning_unit_cell[1] < -1
							|| scanning_unit_cell[2] < -1) {
							neighborflag = 1;
							outofbounds = 1;
						}
						if (scanning_unit_cell[0] > 1 || scanning_unit_cell[1] > 1
							|| scanning_unit_cell[2] > 1) {
							neighborflag = 1;
							outofbounds = 1;
						}
	
						if (outofbounds == 0) {
							for (int kk = 0; kk < nodes_per_element; kk++) {
								shape_func2 = shape_function(scanning_unit_cell[0], scanning_unit_cell[1], scanning_unit_cell[2], 2, kk + 1);
								scan_position[0] += neighbor_element_positions[polyscan][kk][0] * shape_func2;
								scan_position[1] += neighbor_element_positions[polyscan][kk][1] * shape_func2;
								scan_position[2] += neighbor_element_positions[polyscan][kk][2] * shape_func2;
							}
							delx = current_position[0] - scan_position[0];
							dely = current_position[1] - scan_position[1];
							delz = current_position[2] - scan_position[2];
							distancesq = delx*delx + dely*dely + delz*delz;

							if (distancesq < (cut_global_s + cutoff_skin) * (cut_global_s + cutoff_skin)) {
								if (inner_neigh_index == maxneigh_quad_inner + expansion_count_inner*EXPAND) {
									//expand neighborlist memory structure for additional virtual atoms
									expansion_count_inner += 1;
							
									memory->grow(inner_quad_lists_ucell[iii][neigh_quad_counter], maxneigh_quad_inner + expansion_count_inner*EXPAND, 3, "Pair CAC:cell coords expand");
									if(sector_flag)
									memory->grow(inner_quad_lists_index[iii][neigh_quad_counter], maxneigh_quad_inner + expansion_count_inner*EXPAND, 3, "Pair CAC:cell indexes expand");
									else
									memory->grow(inner_quad_lists_index[iii][neigh_quad_counter], maxneigh_quad_inner + expansion_count_inner*EXPAND, 2, "Pair CAC:cell indexes expand");

								}
								inner_quad_lists_ucell[iii][neigh_quad_counter][inner_neigh_index][0] = scanning_unit_cell[0];
								inner_quad_lists_ucell[iii][neigh_quad_counter][inner_neigh_index][1] = scanning_unit_cell[1];
								inner_quad_lists_ucell[iii][neigh_quad_counter][inner_neigh_index][2] = scanning_unit_cell[2];
								//quad_list_container[iii].inner_list2ucell[neigh_quad_counter].cell_indexes[inner_neigh_index][0] = 0;
								inner_quad_lists_index[iii][neigh_quad_counter][inner_neigh_index][0] = current_list_index;
								inner_quad_lists_index[iii][neigh_quad_counter][inner_neigh_index][1] = polyscan;
                if(sector_flag) inner_quad_lists_index[iii][neigh_quad_counter][inner_neigh_index][2] = 
                  quad_sector_select(scanning_unit_cell[0],scanning_unit_cell[1],scanning_unit_cell[2],current_list_index, polyscan);

								inner_neigh_index++;

								inner_quad_lists_counts[iii][neigh_quad_counter] = inner_neigh_index;
                neighbor_weights[iii][1]++;

							}

							else if (distancesq <  (2*cut_global_s + cutoff_skin)  * (2*cut_global_s + cutoff_skin)) {
								if (outer_neighflag) {
									if (outer_neigh_index == maxneigh_quad_outer + expansion_count_outer*EXPAND) {
										expansion_count_outer += 1;
									
										memory->grow(outer_quad_lists_ucell[iii][neigh_quad_counter], maxneigh_quad_outer + expansion_count_outer*EXPAND, 3, "Pair CAC:cell coords expand");
										if(sector_flag)
										memory->grow(outer_quad_lists_index[iii][neigh_quad_counter], maxneigh_quad_outer + expansion_count_outer*EXPAND, 3, "Pair CAC:cell indexes expand");
										else
										memory->grow(outer_quad_lists_index[iii][neigh_quad_counter], maxneigh_quad_outer + expansion_count_outer*EXPAND, 2, "Pair CAC:cell indexes expand");
										
									}
									outer_quad_lists_ucell[iii][neigh_quad_counter][outer_neigh_index][0] = scanning_unit_cell[0];
									outer_quad_lists_ucell[iii][neigh_quad_counter][outer_neigh_index][1] = scanning_unit_cell[1];
									outer_quad_lists_ucell[iii][neigh_quad_counter][outer_neigh_index][2] = scanning_unit_cell[2];
									//quad_list_container[iii].outer_list2ucell[neigh_quad_counter].cell_indexes[outer_neigh_index][0] = 0;
									outer_quad_lists_index[iii][neigh_quad_counter][outer_neigh_index][0] = current_list_index;
									outer_quad_lists_index[iii][neigh_quad_counter][outer_neigh_index][1] = polyscan;
                  if(sector_flag) outer_quad_lists_index[iii][neigh_quad_counter][outer_neigh_index][2] = 
                    quad_sector_select(scanning_unit_cell[0],scanning_unit_cell[1],scanning_unit_cell[2],current_list_index, polyscan);

									outer_neigh_index++;

									outer_quad_lists_counts[iii][neigh_quad_counter] = outer_neigh_index;
									neighbor_weights[iii][1]++;
								}
							}
						}

					}
				}

			}
		}
		if (neighborflag == 1) {
			neighbor_accumulate(current_position[0], current_position[1]
				, current_position[2], iii, inner_neigh_index, outer_neigh_index);
		}
	}
	else {

		neighbor_accumulate(current_x[0], current_x[1]
			, current_x[2], iii, inner_neigh_index, outer_neigh_index);
	}
	//accumulate weights for outer neigh flag cases
	if(outer_neighflag)
    neighbor_weights[iii][2]+=(outer_quad_lists_counts[iii][neigh_quad_counter]+inner_quad_lists_counts[iii][neigh_quad_counter])*
	inner_quad_lists_counts[iii][neigh_quad_counter]+inner_quad_lists_counts[iii][neigh_quad_counter];
}

//contribute force density from neighboring elements of surface quadrature point
//------------------------------------------------------------------------
//this method is designed for 8 node parallelpiped elements; IT IS NOT GENERAL!!.
void PairCAC::neighbor_accumulate(double x,double y,double z,int iii,int inner_neigh_initial, int outer_neigh_initial){
  int i,j,jj,jnum,sign1,sign2;
  double delx,dely,delz;
  int *jlist,*numneigh,**firstneigh;
  double ****nodal_positions= atom->nodal_positions;
	int *element_type = atom->element_type;
	int **element_scale = atom->element_scale;
	int *poly_count = atom->poly_count;
	int *nodes_count_list = atom->nodes_per_element_list;	
  double distancesq;
  //double surface_normal[3];
  double min_point[3];
  double boxmap_matrix[3][3];
  int outofbounds=0;
  int **neighbor_weights = atom-> neighbor_weights;
  double unit_cell_mapped[3];
  double scanning_unit_cell[3];
	int complete = 0;
	// asa_cg work arrays
	double Work[100];
	double **coords = atom->x;
	long iWork[2];
	double min_distance;
  double shape_func2;
  double scan_position[3];
  double rcut;

	int inner_neigh_index = inner_neigh_initial;
	int outer_neigh_index = outer_neigh_initial;

  //initialize quadrature position vector
  quad_r[0]=x;
  quad_r[1]=y;
  quad_r[2]=z;

  ASA_INT n ;
  n = 2 ; /* problem dimension */
  double xm[2], lo[2], hi[2] ;
	int swap_dof_set[6];
	double swap_distancesq;
	double swap_distancesq_set[6];

  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

    jlist = firstneigh[quad_list_counter];
    jnum = numneigh[quad_list_counter];

	for (jj = 0; jj < jnum; jj++) {
		j = jlist[jj];
		j &= NEIGHMASK;

		complete = 0;
		neighbor_element_positions = nodal_positions[j];
		neighbor_element_type = element_type[j];
		neighbor_element_scale = element_scale[j];
		neigh_poly_count=poly_count[j];
		if (neighbor_element_type != 0) {
			unit_cell_mapped[0] = 2 / double(neighbor_element_scale[0]);
			unit_cell_mapped[1] = 2 / double(neighbor_element_scale[1]);
			unit_cell_mapped[2] = 2 / double(neighbor_element_scale[2]);
      //neighbor_weights[iii][2]+=neigh_poly_count;
			neigh_nodes_per_element=nodes_count_list[neighbor_element_type];

			xm[0] = 0;
			xm[1] = 0;

		//preliminary sorting of surfaces according to initial distances of starting point 0,0
			for (int si = 0; si < 6; si++) {
				dof_surf_list[0] = sort_dof_set[si][0];
				dof_surf_list[1] = sort_dof_set[si][1];
				dof_surf_list[2] = sort_dof_set[si][2];
				dof_surf_list[3] = sort_dof_set[si][3];
				surf_select[0] = sort_surf_set[si][0];
				surf_select[1] = sort_surf_set[si][1];

				if (surf_select[0] == 1 && surf_select[1] == -1) {
					shape_args[0] = -1 + unit_cell_mapped[0] / 2;
					shape_args[1] = xm[0];
					shape_args[2] = xm[1];
				}
				else if (surf_select[0] == 1 && surf_select[1] == 1) {
					shape_args[0] = 1 - unit_cell_mapped[0] / 2;
					shape_args[1] = xm[0];
					shape_args[2] = xm[1];
				}
				else if (surf_select[0] == 2 && surf_select[1] == -1) {
					shape_args[0] = xm[0];
					shape_args[1] = -1 + unit_cell_mapped[1] / 2;
					shape_args[2] = xm[1];
				}
				else if (surf_select[0] == 2 && surf_select[1] == 1) {
					shape_args[0] = xm[0];
					shape_args[1] = 1 - unit_cell_mapped[1] / 2;
					shape_args[2] = xm[1];
				}
				else if (surf_select[0] == 3 && surf_select[1] == -1) {
					shape_args[0] = xm[0];
					shape_args[1] = xm[1];
					shape_args[2] = -1 + unit_cell_mapped[2] / 2;
				}
				else if (surf_select[0] == 3 && surf_select[1] == 1) {
					shape_args[0] = xm[0];
					shape_args[1] = xm[1];
					shape_args[2] = 1 - unit_cell_mapped[2] / 2;
				}

				min_point[0] = 0;
				min_point[1] = 0;
				min_point[2] = 0;
				for (int kk = 0; kk < neigh_nodes_per_element; kk++) {
					shape_func2 = shape_function(shape_args[0], shape_args[1], shape_args[2], 2, kk + 1);
					min_point[0] += neighbor_element_positions[0][kk][0] * shape_func2;
					min_point[1] += neighbor_element_positions[0][kk][1] * shape_func2;
					min_point[2] += neighbor_element_positions[0][kk][2] * shape_func2;//check for error of minimum point later
				}
				delx = x - min_point[0];
				dely = y - min_point[1];
				delz = z - min_point[2];
				distancesq = delx*delx + dely*dely + delz*delz;
				//reshuffle the sorted list if a new minimum surface is found
				swap_distancesq_set[si] = distancesq;
			}

			//double current_distancesq_min = swap_distancesq_set[0];
			
			int scan_found;
			for (int si = 1; si < 6; si++) {
				int swap_min = 0;
				for (int scan = 0; scan < si; scan++)
				{
					if (swap_distancesq_set[si] < swap_distancesq_set[scan]) {
						swap_min = 1;
						scan_found = scan;
						break;
					}
				}
				if (swap_min) {
					for (int swap_surface = si-scan_found; swap_surface > 0; swap_surface--) {
						    swap_distancesq = swap_distancesq_set[si - swap_surface];
							swap_dof_set[0] = sort_dof_set[si - swap_surface][0];
							swap_dof_set[1] = sort_dof_set[si - swap_surface][1];
							swap_dof_set[2] = sort_dof_set[si - swap_surface][2];
							swap_dof_set[3] = sort_dof_set[si - swap_surface][3];
							swap_dof_set[4] = sort_surf_set[si - swap_surface][0];
							swap_dof_set[5] = sort_surf_set[si - swap_surface][1];

							swap_distancesq_set[si - swap_surface] = swap_distancesq_set[si];
							sort_dof_set[si - swap_surface][0] = sort_dof_set[si][0];
							sort_dof_set[si - swap_surface][1] = sort_dof_set[si][1];
							sort_dof_set[si - swap_surface][2] = sort_dof_set[si][2];
							sort_dof_set[si - swap_surface][3] = sort_dof_set[si][3];
							sort_surf_set[si - swap_surface][0] = sort_surf_set[si][0];
							sort_surf_set[si - swap_surface][1] = sort_surf_set[si][1];

							swap_distancesq_set[si]= swap_distancesq;
							sort_dof_set[si][0] = swap_dof_set[0];
							sort_dof_set[si][1] = swap_dof_set[1];
							sort_dof_set[si][2] = swap_dof_set[2];
							sort_dof_set[si][3] = swap_dof_set[3];
							sort_surf_set[si][0] = swap_dof_set[4];
							sort_surf_set[si][1] = swap_dof_set[5];
						
					}
				}
			}

		//test the 6 surfaces for proximity satisfying the cutoff radius

			for (int si = 0; si < 6; si++) {
				if (complete == 0) {
					dof_surf_list[0] = sort_dof_set[si][0];
					dof_surf_list[1] = sort_dof_set[si][1];
					dof_surf_list[2] = sort_dof_set[si][2];
					dof_surf_list[3] = sort_dof_set[si][3];
					surf_select[0] = sort_surf_set[si][0];
					surf_select[1] = sort_surf_set[si][1];
					//find the minimum distance on the surface

					/* allocate arrays for problem solution and bounds */

					xm[0] = 0;
					xm[1] = 0;
					for (i = 0; i < n; i++) lo[i] = (double)-1;
					for (i = 0; i < n; i++) hi[i] = (double)1;

					//clock_t tforce_density_min_e = clock();
					iWork[0] = 0;
					iWork[1] = 0;
					for (int Workcounter = 0; Workcounter < 100; Workcounter++) {
						Work[Workcounter] = 0;
					}

					double unit_cell_min = unit_cell_mapped[0];
					if (unit_cell_min > unit_cell_mapped[1]) unit_cell_min = unit_cell_mapped[1];
					if (unit_cell_min > unit_cell_mapped[2]) unit_cell_min = unit_cell_mapped[2];
					//loop minimum for every poly DOF to ensure minimum
					// run the minimization code
					for (poly_min = 0; poly_min < neigh_poly_count; poly_min++) {
						asa_pointer->call_asa_cg(xm, lo, hi, n, 
							1.e-2*unit_cell_min, NULL, Work, iWork);

						double tol = 0.00001*unit_cell_min;
						if (xm[0] > 1 + tol || xm[1] > 1 + tol || xm[0] < -1 - tol || xm[1] < -1 - tol) {
							error->one(FLERR, "minimum points exceed element domain");
						}
						if (surf_select[0] == 1 && surf_select[1] == -1) {
							shape_args[0] = -1 + unit_cell_mapped[0] / 2;
							shape_args[1] = xm[0];
							shape_args[2] = xm[1];
						}
						else if (surf_select[0] == 1 && surf_select[1] == 1) {
							shape_args[0] = 1 - unit_cell_mapped[0] / 2;
							shape_args[1] = xm[0];
							shape_args[2] = xm[1];
						}
						else if (surf_select[0] == 2 && surf_select[1] == -1) {
							shape_args[0] = xm[0];
							shape_args[1] = -1 + unit_cell_mapped[1] / 2;
							shape_args[2] = xm[1];
						}
						else if (surf_select[0] == 2 && surf_select[1] == 1) {
							shape_args[0] = xm[0];
							shape_args[1] = 1 - unit_cell_mapped[1] / 2;
							shape_args[2] = xm[1];
						}
						else if (surf_select[0] == 3 && surf_select[1] == -1) {
							shape_args[0] = xm[0];
							shape_args[1] = xm[1];
							shape_args[2] = -1 + unit_cell_mapped[2] / 2;
						}
						else if (surf_select[0] == 3 && surf_select[1] == 1) {
							shape_args[0] = xm[0];
							shape_args[1] = xm[1];
							shape_args[2] = 1 - unit_cell_mapped[2] / 2;
						}
						min_point[0] = 0;
						min_point[1] = 0;
						min_point[2] = 0;
						for (int kk = 0; kk < neigh_nodes_per_element; kk++) {
							shape_func2 = shape_function(shape_args[0], shape_args[1], shape_args[2], 2, kk + 1);
							min_point[0] += neighbor_element_positions[poly_min][kk][0] * shape_func2;
							min_point[1] += neighbor_element_positions[poly_min][kk][1] * shape_func2;
							min_point[2] += neighbor_element_positions[poly_min][kk][2] * shape_func2;//check for error of minimum point later
						}

						delx = x - min_point[0];
						dely = y - min_point[1];
						delz = z - min_point[2];
						distancesq = delx*delx + dely*dely + delz*delz;
						if (outer_neighflag) {
							if (distancesq <  (2*cut_global_s + cutoff_skin) *(2*cut_global_s + cutoff_skin)) {
								complete = 1;
								break;
							}
						}
						else {
							if (distancesq < (cut_global_s + cutoff_skin) *(cut_global_s + cutoff_skin)) {
								complete = 1;
								break;
							}
						}
					}
				}
				else { break; }
			}

			if (complete == 1) {
			//compute position of minimum mapped point to obtain geometric values
      sign1=sign2=1;
			if(xm[0]<0) sign1=-1;
			if(xm[1]<0) sign2=-1;
			if (surf_select[0] == 1 && surf_select[1] == -1) {
                
				xm[0] = (int)((xm[0]+sign1) / unit_cell_mapped[1]);
				xm[0] = xm[0] * unit_cell_mapped[1]-sign1;
				xm[1] = (int)((xm[1]+sign2) / unit_cell_mapped[2]);
				xm[1] = xm[1] * unit_cell_mapped[2]-sign2;

				if (xm[0] < 0 && xm[0]>-1)
					xm[0] = xm[0] - 0.5*unit_cell_mapped[1];
				else if (xm[0] < 0)
					xm[0] = xm[0] + 0.5*unit_cell_mapped[1];
				else if (xm[0] > 0 && xm[0] < 1)
					xm[0] = xm[0] + 0.5*unit_cell_mapped[1];
				else
					xm[0] = xm[0] - 0.5*unit_cell_mapped[1];

				if (xm[1] < 0 && xm[1]>-1)
					xm[1] = xm[1] - 0.5*unit_cell_mapped[2];
				else if (xm[1] < 0)
					xm[1] = xm[1] + 0.5*unit_cell_mapped[2];
				else if (xm[1] > 0 && xm[1] < 1)
					xm[1] = xm[1] + 0.5*unit_cell_mapped[2];
				else
					xm[1] = xm[1] - 0.5*unit_cell_mapped[2];

				shape_args[0] = -1 + unit_cell_mapped[0] / 2;
				shape_args[1] = xm[0];
				shape_args[2] = xm[1];

			}
			else if (surf_select[0] == 1 && surf_select[1] == 1) {

				xm[0] = (int)((xm[0]+sign1) / unit_cell_mapped[1]);
				xm[0] = xm[0] * unit_cell_mapped[1]-sign1;
				xm[1] = (int)((xm[1]+sign2) / unit_cell_mapped[2]);
				xm[1] = xm[1] * unit_cell_mapped[2]-sign2;

				if (xm[0] < 0 && xm[0]>-1)
					xm[0] = xm[0] - 0.5*unit_cell_mapped[1];
				else if (xm[0] < 0)
					xm[0] = xm[0] + 0.5*unit_cell_mapped[1];
				else if (xm[0] > 0 && xm[0] < 1)
					xm[0] = xm[0] + 0.5*unit_cell_mapped[1];
				else
					xm[0] = xm[0] - 0.5*unit_cell_mapped[1];

				if (xm[1] < 0 && xm[1]>-1)
					xm[1] = xm[1] - 0.5*unit_cell_mapped[2];
				else if (xm[1] < 0)
					xm[1] = xm[1] + 0.5*unit_cell_mapped[2];
				else if (xm[1] > 0 && xm[1] < 1)
					xm[1] = xm[1] + 0.5*unit_cell_mapped[2];
				else
					xm[1] = xm[1] - 0.5*unit_cell_mapped[2];

				shape_args[0] = 1 - unit_cell_mapped[0] / 2;
				shape_args[1] = xm[0];
				shape_args[2] = xm[1];
			}
			else if (surf_select[0] == 2 && surf_select[1] == -1) {

				xm[0] = (int)((xm[0]+sign1) / unit_cell_mapped[0]);
				xm[0] = xm[0] * unit_cell_mapped[0]-sign1;
				xm[1] = (int)((xm[1]+sign2)/ unit_cell_mapped[2]);
				xm[1] = xm[1] * unit_cell_mapped[2]-sign2;

				if (xm[0] < 0 && xm[0]>-1)
					xm[0] = xm[0] - 0.5*unit_cell_mapped[0];
				else if (xm[0] < 0)
					xm[0] = xm[0] + 0.5*unit_cell_mapped[0];
				else if (xm[0] > 0 && xm[0] < 1)
					xm[0] = xm[0] + 0.5*unit_cell_mapped[0];
				else
					xm[0] = xm[0] - 0.5*unit_cell_mapped[0];

				if (xm[1] < 0 && xm[1]>-1)
					xm[1] = xm[1] - 0.5*unit_cell_mapped[2];
				else if (xm[1] < 0)
					xm[1] = xm[1] + 0.5*unit_cell_mapped[2];
				else if (xm[1] > 0 && xm[1] < 1)
					xm[1] = xm[1] + 0.5*unit_cell_mapped[2];
				else
					xm[1] = xm[1] - 0.5*unit_cell_mapped[2];

				shape_args[0] = xm[0];
				shape_args[1] = -1 + unit_cell_mapped[1] / 2;
				shape_args[2] = xm[1];
			}
			else if (surf_select[0] == 2 && surf_select[1] == 1) {

				xm[0] = (int)((xm[0]+sign1) / unit_cell_mapped[0]);
				xm[0] = xm[0] * unit_cell_mapped[0]-sign1;
				xm[1] = (int)((xm[1]+sign2) / unit_cell_mapped[2]);
				xm[1] = xm[1] * unit_cell_mapped[2]-sign2;

				if (xm[0] < 0 && xm[0]>-1)
					xm[0] = xm[0] - 0.5*unit_cell_mapped[0];
				else if (xm[0] < 0)
					xm[0] = xm[0] + 0.5*unit_cell_mapped[0];
				else if (xm[0] > 0 && xm[0] < 1)
					xm[0] = xm[0] + 0.5*unit_cell_mapped[0];
				else
					xm[0] = xm[0] - 0.5*unit_cell_mapped[0];

				if (xm[1] < 0 && xm[1]>-1)
					xm[1] = xm[1] - 0.5*unit_cell_mapped[2];
				else if (xm[1] < 0)
					xm[1] = xm[1] + 0.5*unit_cell_mapped[2];
				else if (xm[1] > 0 && xm[1] < 1)
					xm[1] = xm[1] + 0.5*unit_cell_mapped[2];
				else
					xm[1] = xm[1] - 0.5*unit_cell_mapped[2];

				shape_args[0] = xm[0];
				shape_args[1] = 1 - unit_cell_mapped[1] / 2;
				shape_args[2] = xm[1];
			}
			else if (surf_select[0] == 3 && surf_select[1] == -1) {

				xm[0] = (int)((xm[0]+sign1) / unit_cell_mapped[0]);
				xm[0] = xm[0] * unit_cell_mapped[0]-sign1;
				xm[1] = (int)((xm[1]+sign2) / unit_cell_mapped[1]);
				xm[1] = xm[1] * unit_cell_mapped[1]-sign2;

				if (xm[0] < 0 && xm[0]>-1)
					xm[0] = xm[0] - 0.5*unit_cell_mapped[0];
				else if (xm[0] < 0)
					xm[0] = xm[0] + 0.5*unit_cell_mapped[0];
				else if (xm[0] > 0 && xm[0] < 1)
					xm[0] = xm[0] + 0.5*unit_cell_mapped[0];
				else
					xm[0] = xm[0] - 0.5*unit_cell_mapped[0];

				if (xm[1] < 0 && xm[1]>-1)
					xm[1] = xm[1] - 0.5*unit_cell_mapped[1];
				else if (xm[1] < 0)
					xm[1] = xm[1] + 0.5*unit_cell_mapped[1];
				else if (xm[1] > 0 && xm[1] < 1)
					xm[1] = xm[1] + 0.5*unit_cell_mapped[1];
				else
					xm[1] = xm[1] - 0.5*unit_cell_mapped[1];

				shape_args[0] = xm[0];
				shape_args[1] = xm[1];
				shape_args[2] = -1 + unit_cell_mapped[2] / 2;
			}
			else if (surf_select[0] == 3 && surf_select[1] == 1) {
				xm[0] = (int)((xm[0]+sign1) / unit_cell_mapped[0]);
				xm[0] = xm[0] * unit_cell_mapped[0]-sign1;
				xm[1] = (int)((xm[1]+sign2) / unit_cell_mapped[1]);
				xm[1] = xm[1] * unit_cell_mapped[1]-sign2;

				if (xm[0] < 0 && xm[0]>-1)
					xm[0] = xm[0] - 0.5*unit_cell_mapped[0];
				else if (xm[0] < 0)
					xm[0] = xm[0] + 0.5*unit_cell_mapped[0];
				else if (xm[0] > 0 && xm[0] < 1)
					xm[0] = xm[0] + 0.5*unit_cell_mapped[0];
				else
					xm[0] = xm[0] - 0.5*unit_cell_mapped[0];

				if (xm[1] < 0 && xm[1]>-1)
					xm[1] = xm[1] - 0.5*unit_cell_mapped[1];
				else if (xm[1] < 0)
					xm[1] = xm[1] + 0.5*unit_cell_mapped[1];
				else if (xm[1] > 0 && xm[1] < 1)
					xm[1] = xm[1] + 0.5*unit_cell_mapped[1];
				else
					xm[1] = xm[1] - 0.5*unit_cell_mapped[1];

				shape_args[0] = xm[0];
				shape_args[1] = xm[1];
				shape_args[2] = 1 - unit_cell_mapped[2] / 2;
			}

			if (xm[0] > 1 || xm[1] > 1 || xm[0] < -1 || xm[1] < -1) {
				error->one(FLERR, "minimum points exceed element domain");
			}
			min_point[0] = 0;
			min_point[1] = 0;
			min_point[2] = 0;
			for (int kk = 0; kk < neigh_nodes_per_element; kk++) {
				shape_func2 = shape_function(shape_args[0], shape_args[1], shape_args[2], 2, kk + 1);
				min_point[0] += neighbor_element_positions[poly_min][kk][0] * shape_func2;
				min_point[1] += neighbor_element_positions[poly_min][kk][1] * shape_func2;
				min_point[2] += neighbor_element_positions[poly_min][kk][2] * shape_func2;//check for error of minimum point later
			}

			delx = x - min_point[0];
			dely = y - min_point[1];
			delz = z - min_point[2];
			distancesq = delx*delx + dely*dely + delz*delz;
			min_distance = sqrt(distancesq);
			if (outer_neighflag) { rcut = (2*cut_global_s + cutoff_skin); }
			else { rcut = (cut_global_s + cutoff_skin); }
			rcut = rcut + min_distance;

			for (int polyscan = 0; polyscan < neigh_poly_count; polyscan++) {

				//try making a boxmap matrix for every type later
				for (int id = 0; id < 3; id++) {
					for (int jd = 0; jd < 3; jd++) {
						boxmap_matrix[id][jd] = 0;
						for (int n = 0; n < neigh_nodes_per_element; n++) {

							boxmap_matrix[id][jd] += neighbor_element_positions[polyscan][n][id]
								* shape_function_derivative(shape_args[0], shape_args[1], shape_args[2],
									2, n + 1, jd + 1);
						}
					}
				}

				// initialize local lattice vector approximation
				double a[3];
				a[0] = unit_cell_mapped[0] * boxmap_matrix[0][0];
				a[1] = unit_cell_mapped[0] * boxmap_matrix[1][0];
				a[2] = unit_cell_mapped[0] * boxmap_matrix[2][0];

				double b[3];
				b[0] = unit_cell_mapped[1] * boxmap_matrix[0][1];
				b[1] = unit_cell_mapped[1] * boxmap_matrix[1][1];
				b[2] = unit_cell_mapped[1] * boxmap_matrix[2][1];

				double c[3];
				c[0] = unit_cell_mapped[2] * boxmap_matrix[0][2];
				c[1] = unit_cell_mapped[2] * boxmap_matrix[1][2];
				c[2] = unit_cell_mapped[2] * boxmap_matrix[2][2];

				//perform gram schmidt orthogonalization of the three vectors
				double norm_a = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
				double norm_b = sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);

				double proj_b2a = (b[0] * a[0] + b[1] * a[1] + b[2] * a[2]) / norm_a;

				//compute component of b normal to a
				double b_orth[3];
				b_orth[0] = b[0] - proj_b2a*a[0] / norm_a;
				b_orth[1] = b[1] - proj_b2a*a[1] / norm_a;
				b_orth[2] = b[2] - proj_b2a*a[2] / norm_a;

				double proj_c2a = (b[0] * a[0] + b[1] * a[1] + b[2] * a[2]) / norm_a;
				double norm_b_orth = sqrt(b_orth[0] * b_orth[0] + b_orth[1] * b_orth[1]
					+ b_orth[2] * b_orth[2]);
				double proj_c2b_orth = (b_orth[0] * c[0] + b_orth[1] * c[1] + b_orth[2] * c[2]) / norm_b_orth;

				//compute component of c normal to a and b
				double c_orth[3];
				c_orth[0] = c[0] - proj_c2a*a[0] / norm_a - proj_c2b_orth*b_orth[0] / norm_b_orth;
				c_orth[1] = c[1] - proj_c2a*a[1] / norm_a - proj_c2b_orth*b_orth[1] / norm_b_orth;
				c_orth[2] = c[2] - proj_c2a*a[2] / norm_a - proj_c2b_orth*b_orth[2] / norm_b_orth;
				double norm_c_orth = sqrt(c_orth[0] * c_orth[0] + c_orth[1] * c_orth[1]
					+ c_orth[2] * c_orth[2]);

				int w_span = int(rcut / norm_c_orth) + 1;
				int t_upper_limit;
				int t_lower_limit;
				int s_upper_limit;
				int s_lower_limit;

				for (int wcount = -w_span; wcount < w_span + 1; wcount++) {
					t_lower_limit = -int((rcut + proj_c2b_orth*wcount) / norm_b_orth) - 1;
					t_upper_limit = int((rcut - proj_c2b_orth*wcount) / norm_b_orth) + 1;

					for (int tcount = t_lower_limit; tcount < t_upper_limit + 1; tcount++) {
						s_lower_limit = -int((rcut + proj_c2a*wcount + proj_b2a*tcount) / norm_a) - 1;
						s_upper_limit = int((rcut - proj_c2a*wcount - proj_b2a*tcount) / norm_a) + 1;

						for (int scount = s_lower_limit; scount < s_upper_limit + 1; scount++) {
							//scanning around atom
							outofbounds = 0;

							scanning_unit_cell[0] = scount*unit_cell_mapped[0] + shape_args[0];
							scanning_unit_cell[1] = tcount*unit_cell_mapped[1] + shape_args[1];
							scanning_unit_cell[2] = wcount*unit_cell_mapped[2] + shape_args[2];

							scan_position[0] = 0;
							scan_position[1] = 0;
							scan_position[2] = 0;

							if (scanning_unit_cell[0] < -1 || scanning_unit_cell[1] < -1
								|| scanning_unit_cell[2] < -1) {
								outofbounds = 1;
							}
							if (scanning_unit_cell[0] > 1 || scanning_unit_cell[1] > 1
								|| scanning_unit_cell[2] > 1) {
								outofbounds = 1;
							}

							if (outofbounds == 0) {
								for (int kk = 0; kk < neigh_nodes_per_element; kk++) {
									shape_func2 = shape_function(scanning_unit_cell[0], scanning_unit_cell[1], scanning_unit_cell[2], 2, kk + 1);
									scan_position[0] += neighbor_element_positions[polyscan][kk][0] * shape_func2;
									scan_position[1] += neighbor_element_positions[polyscan][kk][1] * shape_func2;
									scan_position[2] += neighbor_element_positions[polyscan][kk][2] * shape_func2;
								}
								delx = x - scan_position[0];
								dely = y - scan_position[1];
								delz = z - scan_position[2];
								distancesq = delx*delx + dely*dely + delz*delz;
								if (distancesq < (cut_global_s + cutoff_skin)  * (cut_global_s + cutoff_skin)) {
									if (inner_neigh_index == maxneigh_quad_inner + expansion_count_inner*EXPAND) {
									//expand neighborlist memory structure for additional virtual atoms
									expansion_count_inner += 1;
						
									memory->grow(inner_quad_lists_ucell[iii][neigh_quad_counter], maxneigh_quad_inner + expansion_count_inner*EXPAND, 3, "Pair CAC:cell coords expand");
									if(sector_flag)
									memory->grow(inner_quad_lists_index[iii][neigh_quad_counter], maxneigh_quad_inner + expansion_count_inner*EXPAND, 3, "Pair CAC:cell indexes expand");
									else
									memory->grow(inner_quad_lists_index[iii][neigh_quad_counter], maxneigh_quad_inner + expansion_count_inner*EXPAND, 2, "Pair CAC:cell indexes expand");

								}
								inner_quad_lists_ucell[iii][neigh_quad_counter][inner_neigh_index][0] = scanning_unit_cell[0];
								inner_quad_lists_ucell[iii][neigh_quad_counter][inner_neigh_index][1] = scanning_unit_cell[1];
								inner_quad_lists_ucell[iii][neigh_quad_counter][inner_neigh_index][2] = scanning_unit_cell[2];
								//quad_list_container[iii].inner_list2ucell[neigh_quad_counter].cell_indexes[inner_neigh_index][0] = 0;
								inner_quad_lists_index[iii][neigh_quad_counter][inner_neigh_index][0] = j;
								inner_quad_lists_index[iii][neigh_quad_counter][inner_neigh_index][1] = polyscan;
                if(sector_flag) inner_quad_lists_index[iii][neigh_quad_counter][inner_neigh_index][2] = 
                  quad_sector_select(scanning_unit_cell[0],scanning_unit_cell[1],scanning_unit_cell[2],j, polyscan);

								inner_neigh_index++;

								inner_quad_lists_counts[iii][neigh_quad_counter] = inner_neigh_index;
                neighbor_weights[iii][1]++;
								}
								else if (distancesq < (2*cut_global_s + cutoff_skin)  * (2*cut_global_s + cutoff_skin)) {
									//complete = 1;
									if (outer_neighflag) {
								if (outer_neigh_index == maxneigh_quad_outer + expansion_count_outer*EXPAND) {
										expansion_count_outer += 1;
									
										memory->grow(outer_quad_lists_ucell[iii][neigh_quad_counter], maxneigh_quad_outer + expansion_count_outer*EXPAND, 3, "Pair CAC:cell coords expand");
										if(sector_flag)
										memory->grow(outer_quad_lists_index[iii][neigh_quad_counter], maxneigh_quad_outer + expansion_count_outer*EXPAND, 3, "Pair CAC:cell indexes expand");
										else
										memory->grow(outer_quad_lists_index[iii][neigh_quad_counter], maxneigh_quad_outer + expansion_count_outer*EXPAND, 2, "Pair CAC:cell indexes expand");
										
									}
									outer_quad_lists_ucell[iii][neigh_quad_counter][outer_neigh_index][0] =  scanning_unit_cell[0];
									outer_quad_lists_ucell[iii][neigh_quad_counter][outer_neigh_index][1] =  scanning_unit_cell[1];
									outer_quad_lists_ucell[iii][neigh_quad_counter][outer_neigh_index][2] =  scanning_unit_cell[2];
									//quad_list_container[iii].outer_list2ucell[neigh_quad_counter].cell_indexes[outer_neigh_index][0] = 0;
									outer_quad_lists_index[iii][neigh_quad_counter][outer_neigh_index][0] = j;
									outer_quad_lists_index[iii][neigh_quad_counter][outer_neigh_index][1] = polyscan;
                  if(sector_flag) outer_quad_lists_index[iii][neigh_quad_counter][outer_neigh_index][2] = 
                    quad_sector_select(scanning_unit_cell[0],scanning_unit_cell[1],scanning_unit_cell[2],j, polyscan);

									outer_neigh_index++;
                 

									outer_quad_lists_counts[iii][neigh_quad_counter] = outer_neigh_index;
									neighbor_weights[iii][1]++;
									}
								}
							}

						}
					}

				}
			}
		}

	}
		else {
			delx = x - coords[j][0];
			dely = y - coords[j][1];
			delz = z - coords[j][2];
			distancesq = delx*delx + dely*dely + delz*delz;
				if (distancesq < (cut_global_s + cutoff_skin)   * (cut_global_s + cutoff_skin)) {
								if (inner_neigh_index == maxneigh_quad_inner + expansion_count_inner*EXPAND) {
									//expand neighborlist memory structure for additional virtual atoms
									expansion_count_inner += 1;
							
									memory->grow(inner_quad_lists_ucell[iii][neigh_quad_counter], maxneigh_quad_inner + expansion_count_inner*EXPAND, 3, "Pair CAC:cell coords expand");
									if(sector_flag)
									memory->grow(inner_quad_lists_index[iii][neigh_quad_counter], maxneigh_quad_inner + expansion_count_inner*EXPAND, 3, "Pair CAC:cell indexes expand");
									else
									memory->grow(inner_quad_lists_index[iii][neigh_quad_counter], maxneigh_quad_inner + expansion_count_inner*EXPAND, 2, "Pair CAC:cell indexes expand");

								}
								inner_quad_lists_ucell[iii][neigh_quad_counter][inner_neigh_index][0] = coords[j][0];
								inner_quad_lists_ucell[iii][neigh_quad_counter][inner_neigh_index][1] = coords[j][1];
								inner_quad_lists_ucell[iii][neigh_quad_counter][inner_neigh_index][2] = coords[j][2];
								//quad_list_container[iii].inner_list2ucell[neigh_quad_counter].cell_indexes[inner_neigh_index][0] = 0;
								inner_quad_lists_index[iii][neigh_quad_counter][inner_neigh_index][0] = j;
								inner_quad_lists_index[iii][neigh_quad_counter][inner_neigh_index][1] = 0;
                if(sector_flag) inner_quad_lists_index[iii][neigh_quad_counter][inner_neigh_index][2] = 0;

								inner_neigh_index++;

								inner_quad_lists_counts[iii][neigh_quad_counter] = inner_neigh_index;
								neighbor_weights[iii][0]++;
			}

			else if (distancesq < (2*cut_global_s + cutoff_skin)  * (2*cut_global_s + cutoff_skin)) {
				if (outer_neighflag) {
						if (outer_neigh_index == maxneigh_quad_outer + expansion_count_outer*EXPAND) {
										expansion_count_outer += 1;
										memory->grow(outer_quad_lists_ucell[iii][neigh_quad_counter], maxneigh_quad_outer + expansion_count_outer*EXPAND, 3, "Pair CAC:cell coords expand");
										if(sector_flag)
										memory->grow(outer_quad_lists_index[iii][neigh_quad_counter], maxneigh_quad_outer + expansion_count_outer*EXPAND, 3, "Pair CAC:cell indexes expand");
										else
										memory->grow(outer_quad_lists_index[iii][neigh_quad_counter], maxneigh_quad_outer + expansion_count_outer*EXPAND, 2, "Pair CAC:cell indexes expand");
									}
									outer_quad_lists_ucell[iii][neigh_quad_counter][outer_neigh_index][0] = coords[j][0];
									outer_quad_lists_ucell[iii][neigh_quad_counter][outer_neigh_index][1] = coords[j][1];
									outer_quad_lists_ucell[iii][neigh_quad_counter][outer_neigh_index][2] = coords[j][2];
									//quad_list_container[iii].outer_list2ucell[neigh_quad_counter].cell_indexes[outer_neigh_index][0] = 0;
									outer_quad_lists_index[iii][neigh_quad_counter][outer_neigh_index][0] = j;
									outer_quad_lists_index[iii][neigh_quad_counter][outer_neigh_index][1] = 0;
                  if(sector_flag) outer_quad_lists_index[iii][neigh_quad_counter][outer_neigh_index][2] = 0;

									outer_neigh_index++;

									outer_quad_lists_counts[iii][neigh_quad_counter] = outer_neigh_index;
									neighbor_weights[iii][0]++;
				}
			}
		}
		
	}
}

//associate a quadrature point sector (region closest to quadrature point in mapped undeformed space) with a coordinate

int PairCAC::quad_sector_select(double s, double t, double w, int eindex,int pindex){
int *element_scale = atom->element_scale[eindex];
int element_type = atom->element_type[eindex];
int init_quad, quad_index;
//Q8 implementation
if(element_type == 1){
double unit_cell_mapped[3];
int n1,n2,n3,ns,nt,nw,s1,s2,s3;
int signs,signt,signw,surfs,surft,surfw;
int node_index[3];
int quad = quadrature_node_count;
n1 = surface_counts[eindex][0];
n2 = surface_counts[eindex][1];
n3 = surface_counts[eindex][2];
int quad_count = quad*quad*quad + 2 * n1*quad*quad + 2 * n2*quad*quad +
		+ 2 * n3*quad*quad + 4 * n1*n2*quad + 4 * n3*n2*quad + 4 * n1*n3*quad
		+ 8 * n1*n2*n3;
init_quad = quad_count * pindex;
unit_cell_mapped[0] = 2 / double(element_scale[0]);
unit_cell_mapped[1] = 2 / double(element_scale[1]);
unit_cell_mapped[2] = 2 / double(element_scale[2]);
signs=signt=signw=1;
if(s<0) signs=-1;
if(t<0) signt=-1;
if(w<0) signw=-1;
//unit cell count of this point
ns = (int((s+signs) / unit_cell_mapped[0]))+signs;
nt = (int((t+signt) / unit_cell_mapped[1]))+signt;
nw = (int((w+signw) / unit_cell_mapped[2]))+signw;
if(ns<0) ns = -ns;
if(nt<0) nt = -nt;
if(nw<0) nw = -nw;

//start by determining if this point lies in a surface quadrature region
//then determine if its within an edge or a corner; assign quad index accordingly
//depending on which quad point is closest
surfs=surft=surfw=0;

if(signs<0&&ns>(element_scale[0]-n1)){
  surfs = -1;
  s1 = element_scale[0]-ns;
  }
if(signs>=0&&ns>(element_scale[0]-n1)){ 
  surfs = 1;
  s1 = element_scale[0]-ns;
  }
if(signt<0&&nt>(element_scale[1]-n2)){
  surft = -1;
  s2 = element_scale[1]-nt;
  }
if(signt>=0&&nt>(element_scale[1]-n2)){ 
  surft = 1;
  s2 = element_scale[1]-nt;
  }
if(signw<0&&nw>(element_scale[2]-n3)){
  surfw = -1;
  s3 = element_scale[2]-nw;
  }
if(signw>=0&&nw>(element_scale[2]-n3)){
  surfw = 1;
  s3 = element_scale[2]-nw;
  }
//determine quadrature point proximity for solely two point quadrature at this point
//error in place in case restriction forgotten when code updated
if(quadrature_node_count!=2)
error->one(FLERR,"Remember to generalize the quadrature proximity calculation in quad_sector_select in pair_cac.cpp (only quad=2 now)");
//interior cell point
if(surfs==0&&surft==0&&surfw==0){
quad_index = init_quad+quad*quad*(signw+1)/2+quad*(signt+1)/2+(signs+1)/2;
return quad_index;
}
//if s surface point
if(surfs!=0&&surft==0&&surfw==0){
	if(surfs==-1) quad_index = init_quad + quad*quad*quad + quad*quad*s1+quad*(signt+1)/2+(signw+1)/2;
  if(surfs==1) quad_index = init_quad + quad*quad*quad + quad*quad*n1+quad*quad*s1+quad*(signt+1)/2+(signw+1)/2;
	return quad_index;
}
//if t surface point
if(surfs==0&&surft!=0&&surfw==0){
	if(surft==-1) quad_index = init_quad + quad*quad*quad + 2*quad*quad*n1+quad*quad*s2+quad*(signs+1)/2+(signw+1)/2;
	if(surft==1) quad_index = init_quad + quad*quad*quad + 2*quad*quad*n1+quad*quad*n2+quad*quad*s2+quad*(signs+1)/2+(signw+1)/2;
	return quad_index;
}
//if w surface point
if(surfs==0&&surft==0&&surfw!=0){
	if(surfw==-1) quad_index = init_quad + quad*quad*quad + 2*quad*quad*(n1+n2)+quad*quad*s3+quad*(signs+1)/2+(signt+1)/2;
	if(surfw==1) quad_index = init_quad + quad*quad*quad + 2*quad*quad*(n1+n2)+quad*quad*n3+quad*quad*s3+quad*(signs+1)/2+(signt+1)/2;
	return quad_index;
}
//if w oriented edge point
if(surfs!=0&&surft!=0&&surfw==0){
quad_index = init_quad + quad*quad*quad + 2*quad*quad*(n1+n2+n3);
if(signs==-1&&signt==-1)
  quad_index += n2*quad*s1+quad*s2+(signw+1)/2;
if(signs==1&&signt==-1)
  quad_index += quad*n1*n2+n2*quad*s1+quad*s2+(signw+1)/2;
if(signs==-1&&signt==1)
  quad_index += 2*quad*n1*n2+n2*quad*s1+quad*s2+(signw+1)/2;
if(signs==1&&signt==1)
  quad_index += 3*quad*n1*n2+n2*quad*s1+quad*s2+(signw+1)/2;
	return quad_index;
}
//if s oriented edge point
if(surfs==0&&surft!=0&&surfw!=0){
quad_index = init_quad + quad*quad*quad + 2*quad*quad*(n1+n2+n3)+4*quad*n1*n2;
if(signw==-1&&signt==-1)
  quad_index += n3*quad*s2+quad*s3+(signs+1)/2;
if(signw==-1&&signt==1)
  quad_index += quad*n2*n3+n3*quad*s2+quad*s3+(signs+1)/2;
if(signw==1&&signt==-1)
  quad_index += 2*quad*n2*n3+n3*quad*s2+quad*s3+(signs+1)/2;
if(signw==1&&signt==1)
  quad_index += 3*quad*n2*n3+n3*quad*s2+quad*s3+(signs+1)/2;
	return quad_index;
}
//if t oriented edge point
if(surfs!=0&&surft==0&&surfw!=0){
quad_index = init_quad + quad*quad*quad + 2*quad*quad*(n1+n2+n3)+4*quad*n1*n2+4*quad*n2*n3;
if(signw==-1&&signt==-1)
  quad_index += n3*quad*s1+quad*s1+(signt+1)/2;
if(signw==-1&&signs==1)
  quad_index += quad*n1*n3+n3*quad*s1+quad*s1+(signt+1)/2;
if(signw==1&&signs==-1)
  quad_index += 2*quad*n1*n3+n3*quad*s1+quad*s1+(signt+1)/2;
if(signw==1&&signs==1)
  quad_index += 3*quad*n1*n3+n3*quad*s1+quad*s1+(signt+1)/2;
	return quad_index;
}
quad_index = init_quad + quad*quad*quad + 2*quad*quad*(n1+n2+n3)+4*quad*(n1*n2+n2*n3+n1*n3);
//check for corner points
if(surfs==-1&&surft==-1&&surfw==-1) quad_index += n3*n2*s1+n3*s2+s3;
if(surfs==1&&surft==-1&&surfw==-1) quad_index += n1*n2*n3+n3*n2*s1+n3*s2+s3;
if(surfs==1&&surft==1&&surfw==-1) quad_index += 2*n1*n2*n3+n3*n2*s1+n3*s2+s3;
if(surfs==-1&&surft==1&&surfw==-1) quad_index += 3*n1*n2*n3+n3*n2*s1+n3*s2+s3;
if(surfs==-1&&surft==-1&&surfw==1) quad_index += 4*n1*n2*n3+n3*n2*s1+n3*s2+s3;
if(surfs==1&&surft==-1&&surfw==1) quad_index += 5*n1*n2*n3+n3*n2*s1+n3*s2+s3;
if(surfs==1&&surft==1&&surfw==1) quad_index += 6*n1*n2*n3+n3*n2*s1+n3*s2+s3;
if(surfs==-1&&surft==1&&surfw==1) quad_index += 7*n1*n2*n3+n3*n2*s1+n3*s2+s3;
return quad_index;
}
return 0;
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
  int neigh_max = inner_quad_lists_counts[iii][neigh_quad_counter];
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
  outer_neigh_max = outer_quad_lists_counts[iii][neigh_quad_counter];
  if(outer_neigh_max > MAXNEIGHBUFF)
  error->one(FLERR,"increase MAXNEIGHBUFF in pair_cac.cpp; this overestimates the maximum possible neighbors for a quadrature point");	  
  }

  int shape_limit[MAXNEIGHBUFF];
  double shaperesult[MAXNEIGHBUFF][MAXESHAPE];
  double **nodes;
  inner_ucells = inner_quad_lists_ucell[iii][neigh_quad_counter];
  inner_indices = inner_quad_lists_index[iii][neigh_quad_counter];
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
  outer_ucells = outer_quad_lists_ucell[iii][neigh_quad_counter];
  outer_indices = outer_quad_lists_index[iii][neigh_quad_counter];
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

/////////////////////////////////////////////////////
void PairCAC::compute_surface_depths(double &scalex, double &scaley, double &scalez,
	int &countx, int &county, int &countz, int flag) {
	int poly = poly_counter;
double unit_cell_mapped[3];
double rcut;
if (outer_neighflag) { rcut = 2 * cut_global_s; }
else{ rcut =  cut_global_s; }
//flag determines the current element type and corresponding procedure to calculate parameters for 
//surface penetration depth to be used when computing force density with influences from neighboring
//elements

unit_cell_mapped[0] = 2 / double(current_element_scale[0]);
unit_cell_mapped[1] = 2 / double(current_element_scale[1]);
unit_cell_mapped[2] = 2 / double(current_element_scale[2]);
	double ds_x = (current_nodal_positions[0][0] - current_nodal_positions[1][0])*
		(current_nodal_positions[0][0] - current_nodal_positions[1][0]);
	double ds_y = (current_nodal_positions[0][1] - current_nodal_positions[1][1])*
		(current_nodal_positions[0][1] - current_nodal_positions[1][1]);
	double ds_z = (current_nodal_positions[0][2] - current_nodal_positions[1][2])*
		(current_nodal_positions[0][2] - current_nodal_positions[1][2]);
	double ds_surf = 2 * rcut / sqrt(ds_x + ds_y + ds_z);
	ds_surf = unit_cell_mapped[0] * (int)(ds_surf / unit_cell_mapped[0]) + unit_cell_mapped[0];

	double dt_x = (current_nodal_positions[0][0] - current_nodal_positions[3][0])*
		(current_nodal_positions[0][0] - current_nodal_positions[3][0]);
	double dt_y = (current_nodal_positions[0][1] - current_nodal_positions[3][1])*
		(current_nodal_positions[0][1] - current_nodal_positions[3][1]);
	double dt_z = (current_nodal_positions[0][2] - current_nodal_positions[3][2])*
		(current_nodal_positions[0][2] - current_nodal_positions[3][2]);

	double dt_surf = 2 * rcut / sqrt(dt_x + dt_y + dt_z);
	dt_surf = unit_cell_mapped[1] * (int)(dt_surf / unit_cell_mapped[1]) + unit_cell_mapped[1];

	double dw_x = (current_nodal_positions[0][0] - current_nodal_positions[4][0])*
		(current_nodal_positions[0][0] - current_nodal_positions[4][0]);
	double dw_y = (current_nodal_positions[0][1] - current_nodal_positions[4][1])*
		(current_nodal_positions[0][1] - current_nodal_positions[3][1]);
	double dw_z = (current_nodal_positions[0][2] - current_nodal_positions[4][2])*
		(current_nodal_positions[0][2] - current_nodal_positions[4][2]);

	double dw_surf = 2 * rcut / sqrt(dw_x + dw_y + dw_z);
	dw_surf = unit_cell_mapped[2] * (int)(dw_surf / unit_cell_mapped[2]) + unit_cell_mapped[2];
	if (ds_surf > 1 ) {
		ds_surf = 1;


	}
	if (dt_surf > 1 ) {
		
		dt_surf = 1;
		

	}
	if ( dw_surf > 1) {
	
		dw_surf = 1;

	}
	if (one_layer_flag) {
		ds_surf = unit_cell_mapped[0];
		dt_surf = unit_cell_mapped[1];
		dw_surf = unit_cell_mapped[2];
	}
	scalex = 1 - ds_surf;
	scaley = 1 - dt_surf;
	scalez = 1 - dw_surf;
	if(ds_surf==1)
	scalex = 0.0;
	if(dt_surf==1)
	scaley = 0.0;
	if(dw_surf==1)
	scalez = 0.0;

	countx = (int)(ds_surf / unit_cell_mapped[0]);
	county = (int)(dt_surf / unit_cell_mapped[1]);
	countz = (int)(dw_surf / unit_cell_mapped[2]);
  if(countx==0) countx=1;
	if(county==0) county=1;
	if(countz==0) countz=1;



}

//////////////////////////////////////////////////////////////

void PairCAC::allocate_quad_neigh_list(int n1,int n2,int n3) {
	int *element_type = atom->element_type;
  int quad = quadrature_node_count;
	int quad_count;
	max_quad_per_element = quad_count = quad*quad*quad + 2 * n1*quad*quad + 
	2 * n2*quad*quad + 2 * n3*quad*quad + 4 * n1*n2*quad + 4 * n3*n2*quad + 
	4 * n1*n3*quad + 8 * n1*n2*n3;
	
	maxneigh_quad_inner += max_expansion_count_inner*EXPAND;
	maxneigh_quad_outer += max_expansion_count_outer*EXPAND;
	max_expansion_count_inner = 0;
	max_expansion_count_outer = 0;
	// initialize quadrature point neighbor list vectors
	if (quad_allocated) {
		for (int init = 0; init < old_atom_count; init++) {

			if (old_atom_etype[init] == 0) {

				memory->destroy(inner_quad_lists_ucell[init][0]);
				memory->destroy(inner_quad_lists_index[init][0]);
				memory->sfree(inner_quad_lists_ucell[init]);
				memory->sfree(inner_quad_lists_index[init]);
				memory->destroy(inner_quad_lists_counts[init]);
				if (outer_neighflag) {

				memory->destroy(outer_quad_lists_ucell[init][0]);
				memory->destroy(outer_quad_lists_index[init][0]);
				memory->sfree(outer_quad_lists_ucell[init]);
				memory->sfree(outer_quad_lists_index[init]);
				memory->destroy(outer_quad_lists_counts[init]);
				}
			}
			else {

				for (int neigh_loop = 0; neigh_loop < old_quad_count; neigh_loop++) {
					memory->destroy(inner_quad_lists_ucell[init][neigh_loop]);
				  memory->destroy(inner_quad_lists_index[init][neigh_loop]);
				}
				memory->sfree(inner_quad_lists_ucell[init]);
				memory->sfree(inner_quad_lists_index[init]);
				memory->destroy(inner_quad_lists_counts[init]);
				if (outer_neighflag) {

					for (int neigh_loop = 0; neigh_loop < old_quad_count; neigh_loop++) {
					memory->destroy(outer_quad_lists_ucell[init][neigh_loop]);
				  memory->destroy(outer_quad_lists_index[init][neigh_loop]);
					}
				memory->sfree(outer_quad_lists_ucell[init]);
				memory->sfree(outer_quad_lists_index[init]);
				memory->destroy(outer_quad_lists_counts[init]);
				}
			}
	
		}

		
    memory->sfree(inner_quad_lists_ucell);
		memory->sfree(inner_quad_lists_index);
		memory->sfree(inner_quad_lists_counts);
		memory->sfree(outer_quad_lists_ucell);
		memory->sfree(outer_quad_lists_index);
		memory->sfree(outer_quad_lists_counts);
	}
	
		inner_quad_lists_ucell= (double ****) memory->smalloc(sizeof(double ***)*atom->nlocal, "Pair CAC:inner_quad_lists_ucell");
		inner_quad_lists_index= (int ****) memory->smalloc(sizeof(int ***)*atom->nlocal, "Pair CAC:inner_quad_lists_index");
		inner_quad_lists_counts= (int **) memory->smalloc(sizeof(int *)*atom->nlocal, "Pair CAC:inner_quad_lists_counts");
		outer_quad_lists_ucell= (double ****) memory->smalloc(sizeof(double ***)*atom->nlocal, "Pair CAC:outer_quad_lists_ucell");
		outer_quad_lists_index= (int ****) memory->smalloc(sizeof(int ***)*atom->nlocal, "Pair CAC:outer_quad_lists_index");
		outer_quad_lists_counts= (int **) memory->smalloc(sizeof(int *)*atom->nlocal, "Pair CAC:outer_quad_lists_counts");
		for (int init = 0; init < atom->nlocal; init++) {
			
			if (element_type[init] == 0) {
	
				memory->create(inner_quad_lists_counts[init],1, "Pair CAC:inner_quad_lists_counts");
				inner_quad_lists_ucell[init]= (double ***) memory->smalloc(sizeof(double **), "Pair CAC:inner_quad_lists_ucell");
		    inner_quad_lists_index[init]= (int ***) memory->smalloc(sizeof(int **), "Pair CAC:inner_quad_lists_index");
        memory->create(inner_quad_lists_ucell[init][0], maxneigh_quad_inner, 3,  "Pair CAC:inner_quad_lists_ucell");
        if(sector_flag)
				memory->create(inner_quad_lists_index[init][0], maxneigh_quad_inner, 3,  "Pair CAC:inner_quad_lists_index");
        else
        memory->create(inner_quad_lists_index[init][0], maxneigh_quad_inner, 2,  "Pair CAC:inner_quad_lists_index");
				if (outer_neighflag) {
				memory->create(outer_quad_lists_counts[init],1, "Pair CAC:outer_quad_lists_counts");
				outer_quad_lists_ucell[init]= (double ***) memory->smalloc(sizeof(double **), "Pair CAC:outer_quad_lists_ucell");
		    outer_quad_lists_index[init]= (int ***) memory->smalloc(sizeof(int **), "Pair CAC:outer_quad_lists_index");
        memory->create(outer_quad_lists_ucell[init][0], maxneigh_quad_outer, 3, "Pair CAC:outer_quad_lists_ucell");
        if(sector_flag)
				memory->create(outer_quad_lists_index[init][0], maxneigh_quad_outer, 3, "Pair CAC:outer_quad_lists_index");
        else
        memory->create(outer_quad_lists_index[init][0], maxneigh_quad_outer, 2, "Pair CAC:outer_quad_lists_index");
				}
			}
			else {
				memory->create(inner_quad_lists_counts[init],quad_count*atom->maxpoly, "Pair CAC:inner_quad_lists_counts");
				inner_quad_lists_ucell[init]= (double ***) memory->smalloc(sizeof(double **)*quad_count*atom->maxpoly, "Pair CAC:inner_quad_lists_ucell");
		    inner_quad_lists_index[init]= (int ***) memory->smalloc(sizeof(int **)*quad_count*atom->maxpoly, "Pair CAC:inner_quad_lists_index");
				 for (int neigh_loop = 0; neigh_loop < quad_count*atom->maxpoly; neigh_loop++) {
				   memory->create(inner_quad_lists_ucell[init][neigh_loop], maxneigh_quad_inner, 3, "Pair CAC:inner_quad_lists_ucell");
           if(sector_flag)
				   memory->create(inner_quad_lists_index[init][neigh_loop], maxneigh_quad_inner, 3, "Pair CAC:inner_quad_lists_index");
           else
           memory->create(inner_quad_lists_index[init][neigh_loop], maxneigh_quad_inner, 2, "Pair CAC:inner_quad_lists_index");
				}
				if (outer_neighflag) {
				 memory->create(outer_quad_lists_counts[init],quad_count*atom->maxpoly, "Pair CAC:inner_quad_lists_counts");
				 outer_quad_lists_ucell[init]= (double ***) memory->smalloc(sizeof(double **)*quad_count*atom->maxpoly, "Pair CAC:inner_quad_lists_ucell");
		     outer_quad_lists_index[init]= (int ***) memory->smalloc(sizeof(int **)*quad_count*atom->maxpoly, "Pair CAC:inner_quad_lists_index");
				 for (int neigh_loop = 0; neigh_loop < quad_count*atom->maxpoly; neigh_loop++) {
				   memory->create(outer_quad_lists_ucell[init][neigh_loop], maxneigh_quad_outer, 3, "Pair CAC:outer_quad_lists_ucell");
           if(sector_flag)
				   memory->create(outer_quad_lists_index[init][neigh_loop], maxneigh_quad_outer, 3, "Pair CAC:outer_quad_lists_index");
           else
           memory->create(outer_quad_lists_index[init][neigh_loop], maxneigh_quad_outer, 2, "Pair CAC:outer_quad_lists_index");
				}
				}
			}
		}
	
	//initialize counts to zero
	for (int init = 0; init < atom->nlocal; init++) {
		if (element_type[init] == 0) {
			inner_quad_lists_counts[init][0] = 0;
			if (outer_neighflag) {
				outer_quad_lists_counts[init][0] = 0;
			}
		
		}
		else {
			for (int init2 = 0; init2 < quad_count*atom->maxpoly; init2++) {
				inner_quad_lists_counts[init][init2] = 0;
				if (outer_neighflag) {
					outer_quad_lists_counts[init][init2] = 0;
				}
			}
		}
	}
	allocate_quad_attribute(n1,n2,n3);
	quad_allocated = 1;
	if(atom->nlocal>old_atom_count)
	memory->grow(old_atom_etype, atom->nlocal, "Pair CAC:old_element_type_map");
	old_atom_count = atom->nlocal;
	old_quad_count = quad_count*atom->maxpoly;

	if(ghost_quad){
	  memory->grow(old_all_atom_etype, atom->nlocal + atom->nghost, "Pair CAC:old_element_type_map");
	  old_all_atom_count = atom->nlocal + atom->nghost;
	}
	
	for (int init = 0; init < atom->nlocal; init++) {
		old_atom_etype[init]= element_type[init];
	}
}

void PairCAC::allocate_surface_counts() {
	if(!ghost_quad){
	memory->grow(surface_counts, atom->nlocal , 3, "Pair CAC:surface_counts");
	memory->grow(interior_scales, atom->nlocal , 3, "Pair CAC:interior_scales");
	nmax_surf = atom->nlocal;
	}
	else{
	memory->grow(surface_counts, atom->nlocal + atom->nghost, 3, "Pair CAC:surface_counts");
	memory->grow(interior_scales, atom->nlocal + atom->nghost, 3, "Pair CAC:interior_scales");
	nmax_surf = atom->nlocal + atom->nghost;	
	}
}

/* ----------------------------------------------------------------------
   grow quadrature data arrays
------------------------------------------------------------------------- */

void PairCAC::grow_quad_data(){
  quadrature_point_max += 1000*atom->maxpoly;	
  memory->grow(quadrature_point_data,quadrature_point_max,7,"pairCAC:quadrature_point_data");
}

/* ----------------------------------------------------------------------
   copy dense arrays to atomvec arrays for energy_force evaluation
------------------------------------------------------------------------- */

void PairCAC::copy_vectors(int copymode){
int *npoly = atom->poly_count;
  int *nodes_per_element_list = atom->nodes_per_element_list;
  int *element_type = atom->element_type;
  double ****nodal_positions = atom->nodal_positions;
  double ****nodal_velocities = atom->nodal_velocities;
  double ****nodal_forces = atom->nodal_forces;
  double *min_x = atom->min_x;
  double *min_v = atom->min_v;
  double *min_f = atom->min_f;
  atom->dense_count=0;

  
  //copy contents of min vectors to the avec arrays and vice versa
  int dense_count_x=0;
  int dense_count_v=0;
  int dense_count_f=0;
  if(copymode==0){
  for(int element_counter=0; element_counter < atom->nlocal; element_counter++){
    for(int poly_counter=0; poly_counter < npoly[element_counter]; poly_counter++){
      for(int node_counter=0; node_counter < nodes_per_element_list[element_type[element_counter]]; node_counter++){
         nodal_positions[element_counter][poly_counter][node_counter][0] = min_x[dense_count_x++];
         nodal_positions[element_counter][poly_counter][node_counter][1] = min_x[dense_count_x++];
         nodal_positions[element_counter][poly_counter][node_counter][2] = min_x[dense_count_x++];
         nodal_forces[element_counter][poly_counter][node_counter][0] = min_f[dense_count_f++];
         nodal_forces[element_counter][poly_counter][node_counter][1] = min_f[dense_count_f++];
         nodal_forces[element_counter][poly_counter][node_counter][2] = min_f[dense_count_f++];
       }
     }
  }
  }
  if(copymode==1){
   for(int element_counter=0; element_counter < atom->nlocal; element_counter++){
     atom->dense_count+=3*npoly[element_counter]*nodes_per_element_list[element_type[element_counter]];
  }
  
  //grow the dense aligned vectors
  if(atom->dense_count>densemax){
  min_x = memory->grow(atom->min_x,atom->dense_count,"min_CAC_cg:min_x");
  min_v = memory->grow(atom->min_v,atom->dense_count,"min_CAC_cg:min_x");
  min_f = memory->grow(atom->min_f,atom->dense_count,"min_CAC_cg:min_f");
  densemax=atom->dense_count;
  }

  for(int element_counter=0; element_counter < atom->nlocal; element_counter++){
    for(int poly_counter=0; poly_counter < npoly[element_counter]; poly_counter++){
      for(int node_counter=0; node_counter < nodes_per_element_list[element_type[element_counter]]; node_counter++){
         min_x[dense_count_x++] = nodal_positions[element_counter][poly_counter][node_counter][0];
         min_x[dense_count_x++] = nodal_positions[element_counter][poly_counter][node_counter][1];
         min_x[dense_count_x++] = nodal_positions[element_counter][poly_counter][node_counter][2];
		     min_v[dense_count_v++] = nodal_velocities[element_counter][poly_counter][node_counter][0];
         min_v[dense_count_v++] = nodal_velocities[element_counter][poly_counter][node_counter][1];
         min_v[dense_count_v++] = nodal_velocities[element_counter][poly_counter][node_counter][2];
         min_f[dense_count_f++] = nodal_forces[element_counter][poly_counter][node_counter][0];
         min_f[dense_count_f++] = nodal_forces[element_counter][poly_counter][node_counter][1];
         min_f[dense_count_f++] = nodal_forces[element_counter][poly_counter][node_counter][2];
       }
     }
  }
  }

}

//memory usage due to quadrature point list memory structure
double PairCAC::memory_usage()
{
	
    double bytes_used = 0;
    if (quad_allocated) {
		for (int init = 0; init < old_atom_count; init++) {

			if (old_atom_etype[init] == 0) {
				bytes_used +=memory->usage(inner_quad_lists_ucell[init][0], 3, maxneigh_quad_inner);
        if(sector_flag)
				bytes_used +=memory->usage(inner_quad_lists_index[init][0], 3, maxneigh_quad_inner);
        else
        bytes_used +=memory->usage(inner_quad_lists_index[init][0], 2, maxneigh_quad_inner);
				if(outer_neighflag){
				  bytes_used +=memory->usage(outer_quad_lists_ucell[init][0], 3, maxneigh_quad_outer);
          if(sector_flag)
				  bytes_used +=memory->usage(outer_quad_lists_index[init][0], 3, maxneigh_quad_outer);
          else
          bytes_used +=memory->usage(outer_quad_lists_index[init][0], 2, maxneigh_quad_outer);
				}
				//bytes_used +=memory->usage(quad_list_container[init],1);
				
			}
			else {

				for (int neigh_loop = 0; neigh_loop < old_quad_count; neigh_loop++) {
					bytes_used +=memory->usage(inner_quad_lists_ucell[init][neigh_loop], maxneigh_quad_inner, 3);
          if(sector_flag)
					bytes_used +=memory->usage(inner_quad_lists_index[init][neigh_loop], maxneigh_quad_inner, 3);
          else
          bytes_used +=memory->usage(inner_quad_lists_index[init][neigh_loop], maxneigh_quad_inner, 2);
					if(outer_neighflag){
				    bytes_used +=memory->usage(outer_quad_lists_ucell[init][neigh_loop], maxneigh_quad_outer, 3);
            if(sector_flag)
				    bytes_used +=memory->usage(outer_quad_lists_index[init][neigh_loop], maxneigh_quad_outer, 3);
            else
            bytes_used +=memory->usage(outer_quad_lists_index[init][neigh_loop], maxneigh_quad_outer, 2);
				    }
				}
				//bytes_used +=memory->usage(quad_list_container[init].list2ucell,1);
				
			
			}
		}
	}
    
  return bytes_used;

}