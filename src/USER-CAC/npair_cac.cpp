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
#include "npair_cac.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "atom.h"
#include "comm.h"
#include "atom_vec.h"
#include "molecule.h"
#include "domain.h"
#include "my_page.h"
#include "error.h"
#include "memory.h"
#include "asa_user.h"
#include "asa_data.h"
#include "nbin.h"

#define MAXNEIGH  10
#define EXPAND 50
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NPairCAC::NPairCAC(LAMMPS *lmp) : NPair(lmp) {
  nmax=0;
  maxneigh = MAXNEIGH;
  surface_counts = NULL;
  list_container = NULL;
  interior_scales = NULL;
  scan_flags = NULL;
  bin_scan_flags = NULL;
  bins_searched = NULL;
  interface_flags = NULL;
  quadrature_point_data = NULL;
  quadrature_counts = NULL;
  e2quad_index = NULL;
  current_quad_list = NULL;
  neighbor_weights = NULL;
  quad_scan_list_max = 0;
  neigh_allocated = quad_allocated = 0;
  quadrature_point_max = quadrature_poly_count = quadrature_poly_max = 0;
  quad_count_max = max_atom_count = max_quad_per_element = 0;
  maxbins_searched = 0;
  nbins_searched = 0;
  nmax = 0;

  surface_counts_max[0] = 1;
  surface_counts_max[1] = 1;
  surface_counts_max[2] = 1;
  surface_counts_max_old[0] = 1;
  surface_counts_max_old[1] = 1;
  surface_counts_max_old[2] = 1;
  //instance asa interface object
  asa_pointer = new Asa_Data(lmp, this);
  atom->npair_cac = this;

  //initiate arrays used in the neighbor_accumulate function
  memory->create(surf_set, 6, 2, "pairCAC:surf_set");
  memory->create(dof_set, 6, 4, "pairCAC:surf_set");
  memory->create(sort_surf_set, 6, 2, "pairCAC:surf_set");
  memory->create(sort_dof_set, 6, 4, "pairCAC:surf_set");
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

  for (int si = 0; si < 6; si++) {
    sort_dof_set[si][0] = dof_set[si][0];
    sort_dof_set[si][1] = dof_set[si][1];
    sort_dof_set[si][2] = dof_set[si][2];
    sort_dof_set[si][3] = dof_set[si][3];
    sort_surf_set[si][0] = surf_set[si][0];
    sort_surf_set[si][1] = surf_set[si][1];
  }
}

/* ---------------------------------------------------------------------- */

NPairCAC::~NPairCAC() 
{
  memory->destroy(surf_set);
  memory->destroy(dof_set);
  memory->destroy(sort_surf_set);
  memory->destroy(sort_dof_set);
  memory->destroy(surface_counts);
  memory->destroy(interior_scales);
  memory->destroy(scan_flags);
  memory->destroy(bin_scan_flags);
  memory->destroy(bins_searched);
  memory->destroy(interface_flags);
  if (neigh_allocated) {
  for (int init = 0; init < max_atom_count; init++)
    memory->destroy(list_container[init]);	
    memory->sfree(list_container);
  }
  if(quadrature_point_max)
  memory->destroy(quadrature_point_data);
  if(nmax)
  memory->destroy(quadrature_counts);
  if(quad_allocated){
    for(int neigh_loop = 0; neigh_loop < quad_count_max; neigh_loop++){
      memory->destroy(inner_quad_lists_ucell[neigh_loop]);
      memory->destroy(inner_quad_lists_index[neigh_loop]);
    }
    memory->sfree(inner_quad_lists_ucell);
    memory->sfree(inner_quad_lists_index);
    memory->sfree(inner_quad_lists_counts);
    memory->sfree(inner_quad_neigh_maxes);
    if(outer_neigh_flag){
      for(int neigh_loop = 0; neigh_loop < quad_count_max; neigh_loop++){
        memory->destroy(outer_quad_lists_ucell[neigh_loop]);
        memory->destroy(outer_quad_lists_index[neigh_loop]);
      }
      memory->sfree(outer_quad_lists_ucell);
      memory->sfree(outer_quad_lists_index);
      memory->sfree(outer_quad_lists_counts);
      memory->sfree(outer_quad_neigh_maxes);
    }
  }
  delete asa_pointer;
}

/* ----------------------------------------------------------------------
   binned neighbor list construction for all neighbors
   Neighboring is based on the overlap between element bounding boxes + cutoff bound
   and points (quadrature points or pure atoms)
------------------------------------------------------------------------- */

void NPairCAC::build(NeighList *list)
{
  int i, j, k, itype, jtype, ibin, iatom;
  double delx, dely, delz, rsq, s, t, w;
  int *neighptr;

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  tagint *tag = atom->tag;
  tagint *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  double interior_scale[3];
  if (includegroup) nlocal = atom->nfirst;
  double ****nodal_positions = atom->nodal_positions;
  double quad_position[3], shape_func;
  cutneighmax = neighbor->cutneighmax;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int *element_type = atom->element_type;
  int *poly_count = atom->poly_count;
  int **element_scale = atom->element_scale;
  int current_element_type;
  int neighbor_element_type;
  int nsum = atom->nlocal;
  int *nodes_per_element_list = atom->nodes_per_element_list;
  int nodes_per_element;
  double cutoff_skin = neighbor->skin;
  outer_neigh_flag = atom->outer_neigh_flag;
  sector_flag = atom->sector_flag;
  ghost_quad_flag = atom->ghost_quad_flag;
  if(ghost_quad_flag) nsum = atom->nlocal + atom->nghost;

  cut_global = cutneighmax - cutoff_skin;
  if(outer_neigh_flag) cut_global = (cutneighmax-cutoff_skin)/2;

  //settings for compute cac quad count; provides weight to fix balance
  atom->neigh_weight_flag=1;
  atom->weight_count=nlocal;

  quadrature_init(2);
  
  memory->grow(scan_flags,atom->nlocal + atom->nghost,"scan_flags");
  memory->grow(bin_scan_flags,mbins,"bin_scan_flags");
  memory->grow(interface_flags,nsum,"interface_flags");

  //initialize scan flags to 0
  for (int init=0; init < atom->nlocal + atom->nghost; init++){
    scan_flags[init]=0;
    if(init<nsum) interface_flags[init] = 0;
  }

  //initialize bin scan flags to 0
  for (int init=0; init < mbins; init++){
    bin_scan_flags[init]=0;
  }
  
  // initialize or grow memory for the neighbor list of virtual atoms at quadrature points
  if (atom->nlocal)
  allocate_neigh_list();

  /*compute neighbor list of atoms and elements first; use this to decide
  the distribution of quadrature points*/

  //loop over elements
  for (i = 0; i < nlocal; i++) {
  itype = type[i];
  current_element_type = element_type[i];
    
  // loop over all atoms in surrounding bins in stencil including self
  // skip i = j
  //loop over quadrature points of elements, by convention an atom is one quadrature point
  neigh_count = 0;
  nbins_searched = 0;
  for(int ocount=0; ocount < nbin_element_overlap[i]; ocount++){
    ibin = bin_element_overlap[i][ocount];
  for (k = 0; k < nstencil; k++) {
    if(bin_scan_flags[ibin + stencil[k]]) continue;

    if(nbins_searched==maxbins_searched){
      maxbins_searched+=EXPAND;
      memory->grow(bins_searched, maxbins_searched, "NPairCAC:bins_searched");
    }
    bin_scan_flags[ibin + stencil[k]] = 1;
    bins_searched[nbins_searched++] = ibin+stencil[k];

    for (int jj = 0; jj < bin_ncontent[ibin + stencil[k]]; jj++) {
      if(ibin + stencil[k]<0) error->one(FLERR," negative bin index");
      if(ibin + stencil[k]>=mbins) error->one(FLERR," excessive bin index");
      j = bin_content[ibin + stencil[k]][jj];
      neighbor_element_type = element_type[j];
      if (i == j) continue;
      //checks for a possible repeat element in this neighbor list
      if(neighbor_element_type!=0)
        if(scan_flags[j]) continue;
        
      jtype = type[j];
      
      if (exclude && exclusion(i, j, itype, jtype, mask, molecule)) continue;
      
      //check if array is large enough; allocate more if needed
      if (neigh_count == list_maxes[i]) {
          list_maxes[i] += EXPAND;
          memory->grow(list_container[i], list_maxes[i], "NPair CAC:list_container");
      }
      
      if (neighbor_element_type != 0||current_element_type != 0) {
      if(neighbor_element_type != 0 && current_element_type != 0){
        if(CAC_decide_element2element(i,j)){
        scan_flags[j]=1;
        list_container[i][neigh_count++] = j;
        }
      }
      else if(neighbor_element_type != 0){
        if(CAC_decide_element2element(j,i)){
        scan_flags[j]=1;
        list_container[i][neigh_count++] = j;
        }
      }
      else{
        if(CAC_decide_element2element(i,j))
        list_container[i][neigh_count++] = j;
      }
      }
      else{
      delx = x[i][0] - x[j][0];
      dely = x[i][1] - x[j][1];
      delz = x[i][2] - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
        if (rsq <= cutneighmax*cutneighmax)
        list_container[i][neigh_count++] = j;
      }
  }   
  }
  //reset bin scan flags
  for(int reset = 0; reset < nbins_searched; reset++)
  bin_scan_flags[bins_searched[reset]] = 0;
  }
  numneigh[i] = neigh_count;
  
  //reset scan flags for next quadrture point searching for its neighbors
  for(int reset=0; reset < neigh_count; reset++){
    scan_flags[list_container[i][reset]]=0;
  }
  
  }

  //determine which elements need more quadrature points due to interfacing with smaller resolutions
  for (i = 0; i < nlocal; i++) {
    if(!element_type[i]) continue;
    for (int jj = 0; jj < numneigh[i]; jj++){
      j = list_container[i][jj];
      //determination used by product of all three element scales (volumetric comparison)
      //not a rigorous choice but takes care of most obvious cases that will arise; problematic
      //for models with elements that are slender in a given dimension
      if(element_scale[i][0]*element_scale[i][1]*element_scale[i][2]>2*element_scale[j][0]*element_scale[j][1]*element_scale[j][2]){
        interface_flags[i] = 1;
        break;
      }
    }
  }

  //communicate interface flags for ghosts if ghost quadrature points are defined
  if(ghost_quad_flag)
  comm->forward_comm_npair(this,0);
  
  if (nsum > nmax) 
    allocate_local_arrays();

  //initialize the neighbor weights used to balance procs to 0
  for (i = 0; i < nsum; i++){  
    neighbor_weights[i][0]=0;
    neighbor_weights[i][1]=0;
    neighbor_weights[i][2]=0;
  }

  //compute surface counts for quadrature neighbor list allocation
  // initialize or grow surface counts array for quadrature scheme
  // along with interior scaling for the quadrature domain
  surface_counts_max_old[0] = surface_counts_max[0];
  surface_counts_max_old[1] = surface_counts_max[1];
  surface_counts_max_old[2] = surface_counts_max[2];
  for (i = 0; i < nsum; i++) {
    if (element_type[i] != 0) {
      for (poly_counter = 0; poly_counter < poly_count[i]; poly_counter++) {
        int poly_surface_count[3];
        current_nodal_positions = nodal_positions[i][poly_counter];
        compute_surface_depths(interior_scale[0], interior_scale[1], interior_scale[2],
        poly_surface_count[0], poly_surface_count[1], poly_surface_count[2], i, element_type[i]);
        if(poly_counter==0){
          surface_counts[i][0]=poly_surface_count[0];
          surface_counts[i][1]=poly_surface_count[1];
          surface_counts[i][2]=poly_surface_count[2];
          interior_scales[i][0]=interior_scale[0];
          interior_scales[i][1]=interior_scale[1];
          interior_scales[i][2]=interior_scale[2];
        }
        else{
          if (poly_surface_count[0] > surface_counts[i][0]){ surface_counts[i][0] = poly_surface_count[0];
            interior_scales[i][0]=interior_scale[0]; }
          if (poly_surface_count[1] > surface_counts[i][1]){ surface_counts[i][1] = poly_surface_count[1];
            interior_scales[i][1]=interior_scale[1]; }
          if (poly_surface_count[2] > surface_counts[i][2]){ surface_counts[i][2] = poly_surface_count[2];
            interior_scales[i][2]=interior_scale[2]; }
          }
        }
      if (surface_counts[i][0] > surface_counts_max[0]) surface_counts_max[0] = surface_counts[i][0];
      if (surface_counts[i][1] > surface_counts_max[1]) surface_counts_max[1] = surface_counts[i][1];
      if (surface_counts[i][2] > surface_counts_max[2]) surface_counts_max[2] = surface_counts[i][2];
    }
  }
  
  //compute quadrature point data for each element
  qi = 0;
  pqi = 0;
  for (i = 0; i < nsum; i++) {
    quadrature_counts[i] = compute_quad_points(i);
    e2quad_index[i] = pqi;
    pqi += quadrature_counts[i]*poly_count[i];
    if(quadrature_counts[i] > max_quad_per_element) max_quad_per_element = quadrature_counts[i];
  }
  quadrature_poly_count = pqi;
  atom->max_quad_per_element = max_quad_per_element;
  allocate_quad_neigh_list();

  /*compute quadrature point virtual atom list in several steps; first bin the quadrature point, search the stencil
  around it for atoms and elements, find virtual atoms on each neighboring element using asa cg*/
  qi = 0;
  pqi = 0;

  //loop over elements
  for (i = 0; i < nlocal; i++) {
  itype = type[i];
  current_element_type = element_type[i];
  nodes_per_element = nodes_per_element_list[current_element_type];

  // loop over all atoms in surrounding bins in stencil including self
  // skip i = j
  //loop over quadrature points of elements, by convention an atom is one quadrature point
  for (int iquad = 0; iquad < quadrature_counts[i]; iquad++){
    for (poly_counter = 0; poly_counter < poly_count[i]; poly_counter++){
    current_nodal_positions = nodal_positions[i][poly_counter];
    neigh_count = 0;
    //compute current position from the quadrature point data (mapped reference coordinate for Q8 case)
    quad_position[0] = 0;
    quad_position[1] = 0;
    quad_position[2] = 0;

    //Q8 case
    if(current_element_type == 1){
    s = quadrature_point_data[qi][0];
    t = quadrature_point_data[qi][1];
    w = quadrature_point_data[qi][2];
      for (int kkk = 0; kkk < nodes_per_element; kkk++) {
        shape_func = shape_function(s, t, w, 2, kkk + 1);
        quad_position[0] += current_nodal_positions[kkk][0] * shape_func;
        quad_position[1] += current_nodal_positions[kkk][1] * shape_func;
        quad_position[2] += current_nodal_positions[kkk][2] * shape_func;
      }
    }
    //atom case
    if(current_element_type==0){
      quad_position[0] = quadrature_point_data[qi][0];
      quad_position[1] = quadrature_point_data[qi][1];
      quad_position[2] = quadrature_point_data[qi][2];
    }
  
  //now that the current quadrature position is computed search its stencil to create the virtual atom list
  ibin = nb->coord2bin(quad_position);
  for (k = 0; k < nstencil; k++) {
    for (int jj = 0; jj < bin_ncontent[ibin + stencil[k]]; jj++) {
      if(ibin + stencil[k]<0) error->one(FLERR," negative bin index");
      if(ibin + stencil[k]>=mbins) error->one(FLERR," excessive bin index");
      j = bin_content[ibin + stencil[k]][jj];
      neighbor_element_type = element_type[j];
      if (i == j) continue;
      if(neighbor_element_type!=0){
        //checks for repeat elements
        if(scan_flags[j]) continue;
      }
        
      jtype = type[j];
      
      if (exclude && exclusion(i, j, itype, jtype, mask, molecule)) continue;
      
      //check if array is large enough; allocate more if needed
      if (neigh_count == quad_scan_list_max) {
          quad_scan_list_max += EXPAND;
          memory->grow(current_quad_list, 
          quad_scan_list_max, "NPair CAC:current_quad_list");
      }
      
      if (neighbor_element_type != 0) {
      if(CAC_decide_quad2element(quad_position,j)){
        scan_flags[j]=1;
        current_quad_list[neigh_count++] = j;
      }
      }
      else if(neighbor_element_type == 0){
      delx = quad_position[0] - x[j][0];
      dely = quad_position[1] - x[j][1];
      delz = quad_position[2] - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      if (rsq <= cutneighmax*cutneighmax)
       current_quad_list[neigh_count++] = j;
      }
    }
  }

  //use the current quad list to construct the virtual atom list with asa cg and local unit cell scans
  quad_list_build(i,quadrature_point_data[qi][0],quadrature_point_data[qi][1],quadrature_point_data[qi][2]);

  //reset scan flags for next quadrture point searching for its neighbors
  for(int reset=0; reset < neigh_count; reset++){
    scan_flags[current_quad_list[reset]]=0;
  }
  pqi++;
  }
qi++;
}
}
memory->destroy(quadrature_abcissae);
}

///////////////////////////////////////////
double NPairCAC::shape_function(double s, double t, double w, int flag, int index) {
  double shape_function = 0;
  if (flag == 2) {

    if (index == 1) {
      shape_function = (1 - s)*(1 - t)*(1 - w) / 8;
    }
    else if (index == 2) {
      shape_function = (1 + s)*(1 - t)*(1 - w) / 8;
    }
    else if (index == 3) {
      shape_function = (1 + s)*(1 + t)*(1 - w) / 8;
    }
    else if (index == 4) {
      shape_function = (1 - s)*(1 + t)*(1 - w) / 8;
    }
    else if (index == 5) {
      shape_function = (1 - s)*(1 - t)*(1 + w) / 8;
    }
    else if (index == 6) {
      shape_function = (1 + s)*(1 - t)*(1 + w) / 8;
    }
    else if (index == 7) {
      shape_function = (1 + s)*(1 + t)*(1 + w) / 8;
    }
    else if (index == 8) {
      shape_function = (1 - s)*(1 + t)*(1 + w) / 8;
    }


  }
  return shape_function;

}

//-------------------------------------------------------------------------

double NPairCAC::shape_function_derivative(double s, double t, double w, int flag, int index, int derivative){
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

//--------------------------------------------------------------------

void NPairCAC::quad_list_build(int iii, double s, double t, double w) {

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
  double cutoff_skin = neighbor->skin;
  int *element_type = atom->element_type;
  int *poly_count = atom->poly_count;
  double **x = atom->x;
  int inner_neigh_index = 0;
  int outer_neigh_index = 0;
  int nodes_per_element;
  int *nodes_count_list = atom->nodes_per_element_list;
  int current_element_type = element_type[iii];
  int current_poly_count = poly_count[iii];
  int *current_element_scale = atom->element_scale[iii];
  double *current_x = x[iii];
  neighbor_element_positions = nodal_positions[iii];
  if (current_element_type==1) {
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

    if (outer_neigh_flag) { rcut = 2 * cut_global + cutoff_skin; }
    else { rcut = cut_global + cutoff_skin; }

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

              if (distancesq < (cut_global + cutoff_skin) * (cut_global + cutoff_skin)) {
                if (inner_neigh_index == inner_quad_neigh_maxes[pqi]) {
                  //expand neighborlist memory structure for additional virtual atoms
                  inner_quad_neigh_maxes[pqi] += EXPAND;
              
                  memory->grow(inner_quad_lists_ucell[pqi], inner_quad_neigh_maxes[pqi], 3, "Pair CAC:cell coords expand");
                  if(sector_flag)
                  memory->grow(inner_quad_lists_index[pqi], inner_quad_neigh_maxes[pqi], 3, "Pair CAC:cell indexes expand");
                  else
                  memory->grow(inner_quad_lists_index[pqi], inner_quad_neigh_maxes[pqi], 2, "Pair CAC:cell indexes expand");

                }
                inner_quad_lists_ucell[pqi][inner_neigh_index][0] = scanning_unit_cell[0];
                inner_quad_lists_ucell[pqi][inner_neigh_index][1] = scanning_unit_cell[1];
                inner_quad_lists_ucell[pqi][inner_neigh_index][2] = scanning_unit_cell[2];
                inner_quad_lists_index[pqi][inner_neigh_index][0] = iii;
                inner_quad_lists_index[pqi][inner_neigh_index][1] = polyscan;
                if(sector_flag) inner_quad_lists_index[pqi][inner_neigh_index][2] = e2quad_index[iii] + current_poly_count*
                  quad_sector_select(scanning_unit_cell[0],scanning_unit_cell[1],scanning_unit_cell[2],iii, polyscan) + polyscan;

                inner_neigh_index++;

                inner_quad_lists_counts[pqi] = inner_neigh_index;
                neighbor_weights[iii][1]++;

              }

              else if (distancesq <  (2*cut_global + cutoff_skin)  * (2*cut_global + cutoff_skin)) {
                if (outer_neigh_flag) {
                  if (outer_neigh_index == outer_quad_neigh_maxes[pqi]) {
                    outer_quad_neigh_maxes[pqi] += EXPAND;
                  
                    memory->grow(outer_quad_lists_ucell[pqi], outer_quad_neigh_maxes[pqi], 3, "Pair CAC:cell coords expand");
                    if(sector_flag)
                    memory->grow(outer_quad_lists_index[pqi], outer_quad_neigh_maxes[pqi], 3, "Pair CAC:cell indexes expand");
                    else
                    memory->grow(outer_quad_lists_index[pqi], outer_quad_neigh_maxes[pqi], 2, "Pair CAC:cell indexes expand");
                    
                  }
                  outer_quad_lists_ucell[pqi][outer_neigh_index][0] = scanning_unit_cell[0];
                  outer_quad_lists_ucell[pqi][outer_neigh_index][1] = scanning_unit_cell[1];
                  outer_quad_lists_ucell[pqi][outer_neigh_index][2] = scanning_unit_cell[2];
                  outer_quad_lists_index[pqi][outer_neigh_index][0] = iii;
                  outer_quad_lists_index[pqi][outer_neigh_index][1] = polyscan;
                  if(sector_flag) outer_quad_lists_index[pqi][outer_neigh_index][2] = e2quad_index[iii] + current_poly_count*
                  quad_sector_select(scanning_unit_cell[0],scanning_unit_cell[1],scanning_unit_cell[2],iii, polyscan) + polyscan;

                  outer_neigh_index++;

                  outer_quad_lists_counts[pqi] = outer_neigh_index;
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
  if(outer_neigh_flag)
    neighbor_weights[iii][2]+=(outer_quad_lists_counts[pqi]+inner_quad_lists_counts[pqi])*
  inner_quad_lists_counts[pqi]+inner_quad_lists_counts[pqi];
}

//contribute force density from neighboring elements of surface quadrature point
//------------------------------------------------------------------------
//this method is designed for 8 node parallelpiped elements; IT IS NOT GENERAL!!.
void NPairCAC::neighbor_accumulate(double x,double y,double z,int iii,int inner_neigh_initial, int outer_neigh_initial){
  int i,j,jj,jnum,sign1,sign2,flag, check_flag;
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
  double unit_cell_mapped[3];
  double scanning_unit_cell[3];
  int complete = 0;
  // asa_cg work arrays
  double Work[100];
  double **coords = atom->x;
  long iWork[2];
  double min_distance;
  double shape_func2, shape_args[3];
  double scan_position[3];
  double rcut;
  double cutoff_skin = neighbor->skin;
  int current_element_type = element_type[iii];
  int inner_neigh_index = inner_neigh_initial;
  int outer_neigh_index = outer_neigh_initial;
  int neigh_poly_count;

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

  for (jj = 0; jj < neigh_count; jj++) {
    j = current_quad_list[jj];
    j &= NEIGHMASK;
    neigh_poly_count = poly_count[j];
    complete = 0;
    neighbor_element_positions = nodal_positions[j];
    neighbor_element_type = element_type[j];
    neighbor_element_scale = element_scale[j];
    neigh_poly_count=poly_count[j];
    check_flag = 0;
    if (neighbor_element_type == 1) {
      unit_cell_mapped[0] = 2 / double(neighbor_element_scale[0]);
      unit_cell_mapped[1] = 2 / double(neighbor_element_scale[1]);
      unit_cell_mapped[2] = 2 / double(neighbor_element_scale[2]);
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
            flag = asa_pointer->call_asa_cg(xm, lo, hi, n, 
              1.e-2*unit_cell_min, NULL, Work, iWork);
              if(flag==7) check_flag = 1;

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
            if (outer_neigh_flag) {
              if (distancesq <  (2*cut_global + cutoff_skin) *(2*cut_global + cutoff_skin)) {
                complete = 1;
                break;
              }
            }
            else {
              if (distancesq < (cut_global + cutoff_skin) *(cut_global + cutoff_skin)) {
                complete = 1;
                break;
              }
            }
          }
        }
        else { break; }
      }
      //check if flag was ever 7 and no completion was reached
      if(complete==0){
        if(check_flag)
          error->one(FLERR,"asa_cg iterations failed in finding a solution for quadrature neighboring processes. flag = 7");
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
      if (outer_neigh_flag) { rcut = (2*cut_global + cutoff_skin); }
      else { rcut = (cut_global + cutoff_skin); }
      rcut = rcut + min_distance;

      for (int polyscan = 0; polyscan < neigh_poly_count; polyscan++) {

        //try making a boxmap matrix for every type later
        for (int id = 0; id < 3; id++) {
          for (int jd = 0; jd < 3; jd++) {
            boxmap_matrix[id][jd] = 0;
            for (int ni = 0; ni < neigh_nodes_per_element; ni++) {

              boxmap_matrix[id][jd] += neighbor_element_positions[polyscan][ni][id]
                * shape_function_derivative(shape_args[0], shape_args[1], shape_args[2],
                  2, ni + 1, jd + 1);
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
                if (distancesq < (cut_global + cutoff_skin)  * (cut_global + cutoff_skin)) {
                  if (inner_neigh_index == inner_quad_neigh_maxes[pqi]) {
                  //expand neighborlist memory structure for additional virtual atoms
                  inner_quad_neigh_maxes[pqi] += EXPAND;
            
                  memory->grow(inner_quad_lists_ucell[pqi], inner_quad_neigh_maxes[pqi], 3, "Pair CAC:cell coords expand");
                  if(sector_flag)
                  memory->grow(inner_quad_lists_index[pqi], inner_quad_neigh_maxes[pqi], 3, "Pair CAC:cell indexes expand");
                  else
                  memory->grow(inner_quad_lists_index[pqi], inner_quad_neigh_maxes[pqi], 2, "Pair CAC:cell indexes expand");

                }
                inner_quad_lists_ucell[pqi][inner_neigh_index][0] = scanning_unit_cell[0];
                inner_quad_lists_ucell[pqi][inner_neigh_index][1] = scanning_unit_cell[1];
                inner_quad_lists_ucell[pqi][inner_neigh_index][2] = scanning_unit_cell[2];
                inner_quad_lists_index[pqi][inner_neigh_index][0] = j;
                inner_quad_lists_index[pqi][inner_neigh_index][1] = polyscan;
                if(sector_flag) inner_quad_lists_index[pqi][inner_neigh_index][2] = e2quad_index[j] + poly_count[j]*
                  quad_sector_select(scanning_unit_cell[0],scanning_unit_cell[1],scanning_unit_cell[2],j, polyscan) + polyscan;

                inner_neigh_index++;

                inner_quad_lists_counts[pqi] = inner_neigh_index;
                neighbor_weights[iii][1]++;
                }
                else if (distancesq < (2*cut_global + cutoff_skin)  * (2*cut_global + cutoff_skin)) {
                  //complete = 1;
                  if (outer_neigh_flag) {
                if (outer_neigh_index == outer_quad_neigh_maxes[pqi]) {
                    outer_quad_neigh_maxes[pqi] += EXPAND;
                  
                    memory->grow(outer_quad_lists_ucell[pqi], outer_quad_neigh_maxes[pqi], 3, "Pair CAC:cell coords expand");
                    if(sector_flag)
                    memory->grow(outer_quad_lists_index[pqi], outer_quad_neigh_maxes[pqi], 3, "Pair CAC:cell indexes expand");
                    else
                    memory->grow(outer_quad_lists_index[pqi], outer_quad_neigh_maxes[pqi], 2, "Pair CAC:cell indexes expand");
                    
                  }
                  outer_quad_lists_ucell[pqi][outer_neigh_index][0] =  scanning_unit_cell[0];
                  outer_quad_lists_ucell[pqi][outer_neigh_index][1] =  scanning_unit_cell[1];
                  outer_quad_lists_ucell[pqi][outer_neigh_index][2] =  scanning_unit_cell[2];
                  outer_quad_lists_index[pqi][outer_neigh_index][0] = j;
                  outer_quad_lists_index[pqi][outer_neigh_index][1] = polyscan;
                  if(sector_flag) outer_quad_lists_index[pqi][outer_neigh_index][2] = e2quad_index[j] + poly_count[j]*
                  quad_sector_select(scanning_unit_cell[0],scanning_unit_cell[1],scanning_unit_cell[2],j, polyscan) + polyscan;

                  outer_neigh_index++;

                  outer_quad_lists_counts[pqi] = outer_neigh_index;
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
        if (distancesq < (cut_global + cutoff_skin)   * (cut_global + cutoff_skin)) {
                if (inner_neigh_index == inner_quad_neigh_maxes[pqi]) {
                  //expand neighborlist memory structure for additional virtual atoms
                  inner_quad_neigh_maxes[pqi] += EXPAND;
              
                  memory->grow(inner_quad_lists_ucell[pqi], inner_quad_neigh_maxes[pqi], 3, "Pair CAC:cell coords expand");
                  if(sector_flag)
                  memory->grow(inner_quad_lists_index[pqi], inner_quad_neigh_maxes[pqi], 3, "Pair CAC:cell indexes expand");
                  else
                  memory->grow(inner_quad_lists_index[pqi], inner_quad_neigh_maxes[pqi], 2, "Pair CAC:cell indexes expand");

                }
                inner_quad_lists_ucell[pqi][inner_neigh_index][0] = coords[j][0];
                inner_quad_lists_ucell[pqi][inner_neigh_index][1] = coords[j][1];
                inner_quad_lists_ucell[pqi][inner_neigh_index][2] = coords[j][2];
                inner_quad_lists_index[pqi][inner_neigh_index][0] = j;
                inner_quad_lists_index[pqi][inner_neigh_index][1] = 0;
                if(sector_flag) inner_quad_lists_index[pqi][inner_neigh_index][2] = e2quad_index[j];

                inner_neigh_index++;

                inner_quad_lists_counts[pqi] = inner_neigh_index;
                neighbor_weights[iii][0]++;
      }

      else if (distancesq < (2*cut_global + cutoff_skin)  * (2*cut_global + cutoff_skin)) {
        if (outer_neigh_flag) {
            if (outer_neigh_index == outer_quad_neigh_maxes[pqi]) {
                    outer_quad_neigh_maxes[pqi] += EXPAND;

                    memory->grow(outer_quad_lists_ucell[pqi], outer_quad_neigh_maxes[pqi], 3, "Pair CAC:cell coords expand");
                    if(sector_flag)
                    memory->grow(outer_quad_lists_index[pqi], outer_quad_neigh_maxes[pqi], 3, "Pair CAC:cell indexes expand");
                    else
                    memory->grow(outer_quad_lists_index[pqi], outer_quad_neigh_maxes[pqi], 2, "Pair CAC:cell indexes expand");
                  }
                  outer_quad_lists_ucell[pqi][outer_neigh_index][0] = coords[j][0];
                  outer_quad_lists_ucell[pqi][outer_neigh_index][1] = coords[j][1];
                  outer_quad_lists_ucell[pqi][outer_neigh_index][2] = coords[j][2];
                  outer_quad_lists_index[pqi][outer_neigh_index][0] = j;
                  outer_quad_lists_index[pqi][outer_neigh_index][1] = 0;
                  if(sector_flag) outer_quad_lists_index[pqi][outer_neigh_index][2] = e2quad_index[j];

                  outer_neigh_index++;

                  outer_quad_lists_counts[pqi] = outer_neigh_index;
                  neighbor_weights[iii][0]++;
        }
      }
    }
    
  }
}

////////////////////////////////////////////////////////

void NPairCAC::compute_surface_depths(double &scalex, double &scaley, double &scalez,
  int &countx, int &county, int &countz, int eindex, int etype) {
  double unit_cell_mapped[3];
  int *current_element_scale = atom->element_scale[eindex];
  double rcut;
  rcut = cutneighmax; 
  //flag determines the current element type and corresponding procedure to calculate parameters for 
  //surface penetration depth to be used when computing force density with influences from neighboring
  //elements
 if(etype == 1){
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
    (current_nodal_positions[0][1] - current_nodal_positions[4][1]);
  double dw_z = (current_nodal_positions[0][2] - current_nodal_positions[4][2])*
    (current_nodal_positions[0][2] - current_nodal_positions[4][2]);

  double dw_surf = 2 * rcut / sqrt(dw_x + dw_y + dw_z);
  dw_surf = unit_cell_mapped[2] * (int)(dw_surf / unit_cell_mapped[2]) + unit_cell_mapped[2];
  if (ds_surf > 1) {
    ds_surf = 1;
  }
  if (dt_surf > 1) {
    dt_surf = 1;
  }
  if (dw_surf > 1) {
    dw_surf = 1;
  }

  if (atom->one_layer_flag&&!interface_flags[eindex]) {
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
 }
  if(countx==0) countx=1;
  if(county==0) county=1;
  if(countz==0) countz=1;
}

/////////////////////////////////

void NPairCAC::quadrature_init(int quadrature_rank){

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

}

//compute and store quad points for the element referred to by element_index

int NPairCAC::compute_quad_points(int element_index){

  double unit_cell_mapped[3];
  int **element_scale = atom->element_scale;
  int *element_type = atom->element_type;
  int *current_element_scale = atom->element_scale[element_index];
  //int nodes_per_element;
  double coefficients;
  int signs, signt, signw;
  double **x = atom->x;
  unit_cell_mapped[0] = 2 / double(element_scale[element_index][0]);
  unit_cell_mapped[1] = 2 / double(element_scale[element_index][1]);
  unit_cell_mapped[2] = 2 / double(element_scale[element_index][2]);

  double interior_scale[3];
  double *current_x = x[element_index];
  
  int isurface_counts[3];
  interior_scale[0] = interior_scales[element_index][0];
  interior_scale[1] = interior_scales[element_index][1];
  interior_scale[2] = interior_scales[element_index][2];
  isurface_counts[0] = surface_counts[element_index][0];
  isurface_counts[1] = surface_counts[element_index][1];
  isurface_counts[2] = surface_counts[element_index][2];

  double s, t, w;
  s = t = w = 0;
  double sq, tq, wq;
  int neigh_quad_counter = 0;
  //compute interior quadrature point virtual neighbor lists
  if(quadrature_point_max==0) grow_quad_data();

  //Q8 quadrature point selection
  if (element_type[element_index]==1) {
    if(current_element_scale[0]>2&&current_element_scale[1]>2&&current_element_scale[2]>2){
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

          quadrature_point_data[qi][0]=s;
          quadrature_point_data[qi][1]=t;
          quadrature_point_data[qi][2]=w;
          quadrature_point_data[qi][3]=sq;
          quadrature_point_data[qi][4]=tq;
          quadrature_point_data[qi][5]=wq;
          quadrature_point_data[qi][6]=
          interior_scale[0] * interior_scale[1] * interior_scale[2] *
          quadrature_weights[i] * quadrature_weights[j] * quadrature_weights[k];
          neigh_quad_counter = neigh_quad_counter + 1;
          qi+=1;
          if(qi==quadrature_point_max)
          grow_quad_data();
        }
      }
    }
    }


    //compute surface contributions to element

    int sign[2];
    sign[0] = -1;
    sign[1] = 1;

    // s axis surface contributions
    if(current_element_scale[1]>2&&current_element_scale[2]>2){
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

            quadrature_point_data[qi][0]=s;
            quadrature_point_data[qi][1]=t;
            quadrature_point_data[qi][2]=w;
            quadrature_point_data[qi][3]=sq;
            quadrature_point_data[qi][4]=tq;
            quadrature_point_data[qi][5]=wq;
            quadrature_point_data[qi][6]=
              unit_cell_mapped[0] * interior_scale[1] *
              interior_scale[2] * quadrature_weights[j] * quadrature_weights[k];
            neigh_quad_counter = neigh_quad_counter + 1;
            qi+=1;  
            if(qi==quadrature_point_max)
              grow_quad_data();      
          }
        }
      }
    }
    }
    // t axis contributions
    if(current_element_scale[0]>2&&current_element_scale[2]>2){
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
            
            quadrature_point_data[qi][0]=s;
            quadrature_point_data[qi][1]=t;
            quadrature_point_data[qi][2]=w;
            quadrature_point_data[qi][3]=sq;
            quadrature_point_data[qi][4]=tq;
            quadrature_point_data[qi][5]=wq;
            quadrature_point_data[qi][6]=
              unit_cell_mapped[1] * interior_scale[0] *
              interior_scale[2] * quadrature_weights[j] * quadrature_weights[k];
            neigh_quad_counter = neigh_quad_counter + 1;
                  qi+=1;
            if(qi==quadrature_point_max)
              grow_quad_data(); 
          }
        }
      }
    }
    }
    //w axis contributions
    if(current_element_scale[0]>2&&current_element_scale[1]>2){
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

            quadrature_point_data[qi][0]=s;
            quadrature_point_data[qi][1]=t;
            quadrature_point_data[qi][2]=w;
            quadrature_point_data[qi][3]=sq;
            quadrature_point_data[qi][4]=tq;
            quadrature_point_data[qi][5]=wq;
            quadrature_point_data[qi][6]=
              unit_cell_mapped[2] * interior_scale[1] *
              interior_scale[0] * quadrature_weights[j] * quadrature_weights[k];
            neigh_quad_counter = neigh_quad_counter + 1;
                  qi+=1;
            if(qi==quadrature_point_max)
              grow_quad_data(); 
          }
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
        if(current_element_scale[2]<=2) continue;
        interior_scalez = interior_scale[2];
        surface_countsx = isurface_counts[0];
        surface_countsy = isurface_counts[1];
        unit_mappedx = unit_cell_mapped[0];
        unit_mappedy = unit_cell_mapped[1];
      }
      else if (sc == 4 || sc == 5 || sc == 6 || sc == 7) {
        if(current_element_scale[0]<=2) continue;
        interior_scalez = interior_scale[0];
        surface_countsx = isurface_counts[1];
        surface_countsy = isurface_counts[2];
        unit_mappedx = unit_cell_mapped[1];
        unit_mappedy = unit_cell_mapped[2];
      }
      else if (sc == 8 || sc == 9 || sc == 10 || sc == 11) {
        if(current_element_scale[1]<=2) continue;
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
                        
            quadrature_point_data[qi][0]=s;
            quadrature_point_data[qi][1]=t;
            quadrature_point_data[qi][2]=w;
            quadrature_point_data[qi][3]=sq;
            quadrature_point_data[qi][4]=tq;
            quadrature_point_data[qi][5]=wq;
            quadrature_point_data[qi][6]=
            coefficients = unit_mappedx * unit_mappedy * interior_scalez *
            quadrature_weights[k];
            neigh_quad_counter = neigh_quad_counter + 1;
            qi+=1;
            if(qi==quadrature_point_max)
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
                        
            quadrature_point_data[qi][0]=s;
            quadrature_point_data[qi][1]=t;
            quadrature_point_data[qi][2]=w;
            quadrature_point_data[qi][3]=s;
            quadrature_point_data[qi][4]=t;
            quadrature_point_data[qi][5]=w;
            quadrature_point_data[qi][6]=
            unit_cell_mapped[0] * unit_cell_mapped[1] * unit_cell_mapped[2];	
            neigh_quad_counter = neigh_quad_counter + 1;
            qi+=1;
            if(qi==quadrature_point_max)
              grow_quad_data(); 
          }

        }
      }
    }
  }
  //end of Q8 quadrature point selection
  //atom point selection is trivial; assigned to array for convenience at this point.
  else {
    quadrature_point_data[qi][0]=current_x[0];
    quadrature_point_data[qi][1]=current_x[1];
    quadrature_point_data[qi][2]=current_x[2];
    quadrature_point_data[qi][3]=current_x[0];
    quadrature_point_data[qi][4]=current_x[1];
    quadrature_point_data[qi][5]=current_x[2];
    quadrature_point_data[qi][6]= 1;
    qi+=1;
    if(qi==quadrature_point_max)
    grow_quad_data(); 
    neigh_quad_counter = 1;
  }

  return neigh_quad_counter;
}

//decide if an element is close enough to an element or atom to add as a neighbor

int NPairCAC::CAC_decide_element2element(int origin_element_index, int neighbor_element_index) {
  int found_flag = 0;
  double **eboxes = atom->eboxes;
  int *element_type = atom->element_type;
  int *ebox_ref=atom->ebox_ref;
  double *cutghost = comm->cutghost;
  double **x = atom->x;
  double bounding_boxlo[3], neighbounding_boxlo[3];
  double bounding_boxhi[3], neighbounding_boxhi[3];
  int box_overlap_flag[3];
  double *neigh_ebox;
  double *origin_ebox;
  int neighboring_element_type = element_type[neighbor_element_index];

  origin_ebox = eboxes[ebox_ref[origin_element_index]];
  bounding_boxlo[0] = origin_ebox[0];
  bounding_boxlo[1] = origin_ebox[1];
  bounding_boxlo[2] = origin_ebox[2];
  bounding_boxhi[0] = origin_ebox[3];
  bounding_boxhi[1] = origin_ebox[4];
  bounding_boxhi[2] = origin_ebox[5];
  
  if(neighboring_element_type!=0){
  neigh_ebox = eboxes[ebox_ref[neighbor_element_index]];
  neighbounding_boxlo[0] = neigh_ebox[0]+cutghost[0];
  neighbounding_boxlo[1] = neigh_ebox[1]+cutghost[1];
  neighbounding_boxlo[2] = neigh_ebox[2]+cutghost[2];
  neighbounding_boxhi[0] = neigh_ebox[3]-cutghost[0];
  neighbounding_boxhi[1] = neigh_ebox[4]-cutghost[1];
  neighbounding_boxhi[2] = neigh_ebox[5]-cutghost[2];
  
  box_overlap_flag[2]=box_overlap_flag[1]=box_overlap_flag[0]=0;
    for(int idim=0; idim<domain->dimension; idim++){
      if(bounding_boxlo[idim]>=neighbounding_boxlo[idim]&&bounding_boxlo[idim]<=neighbounding_boxhi[idim])
        box_overlap_flag[idim]=1;
        if(bounding_boxhi[idim]>=neighbounding_boxlo[idim]&&bounding_boxhi[idim]<=neighbounding_boxhi[idim])
        box_overlap_flag[idim]=1;
        if(bounding_boxlo[idim]<=neighbounding_boxlo[idim]&&bounding_boxhi[idim]>=neighbounding_boxhi[idim])
        box_overlap_flag[idim]=1;
      }

    if(box_overlap_flag[0]==1&&box_overlap_flag[1]==1&&box_overlap_flag[2]==1)
      found_flag = 1;
  }
  else{
    box_overlap_flag[2]=box_overlap_flag[1]=box_overlap_flag[0]=0;
    for(int idim=0; idim<domain->dimension; idim++){
      if(x[neighbor_element_index][idim]>=bounding_boxlo[idim]&&x[neighbor_element_index][idim]<=bounding_boxhi[idim])
        box_overlap_flag[idim]=1;
      }

    if(box_overlap_flag[0]==1&&box_overlap_flag[1]==1&&box_overlap_flag[2]==1)
      found_flag = 1;
  }
  
  return found_flag;
}

//decide if an atom or quadrature point is close enough to an element to consider for nonlocal quadrature calculation

int NPairCAC::CAC_decide_quad2element(double *current_quad_point, int neighbor_element_index) {
  int found_flag = 0;
  double **eboxes=atom->eboxes;
  int *ebox_ref=atom->ebox_ref;
  double bounding_boxlo[3];
  double bounding_boxhi[3];
  double *current_ebox;
  
  current_ebox = eboxes[ebox_ref[neighbor_element_index]];
  bounding_boxlo[0] = current_ebox[0];
  bounding_boxlo[1] = current_ebox[1];
  bounding_boxlo[2] = current_ebox[2];
  bounding_boxhi[0] = current_ebox[3];
  bounding_boxhi[1] = current_ebox[4];
  bounding_boxhi[2] = current_ebox[5];
  
  if(current_quad_point[0]>bounding_boxlo[0]&&current_quad_point[0]<bounding_boxhi[0]&&
    current_quad_point[1]>bounding_boxlo[1]&&current_quad_point[1]<bounding_boxhi[1]&&
    current_quad_point[2]>bounding_boxlo[2]&&current_quad_point[2]<bounding_boxhi[2])
  {
    found_flag = 1;
    return found_flag;
  }
  else {
    return found_flag;
  }
  
}

//allocate quadrature based neighbor storage
void NPairCAC::allocate_neigh_list() {
  int *element_type = atom->element_type;
  
  // initialize neighbor list arrays
  if(!neigh_allocated){
    atom->list_container = list_container= (int **) memory->smalloc(atom->nlocal*sizeof(int *), "NPair CAC:list_container");
    memory->create(list_maxes, atom->nlocal, "NPair CAC:list_maxes");
      for (int init = 0; init < atom->nlocal; init++) {
        list_container[init] = memory->create(list_container[init],maxneigh, "NPair CAC:cell_indexes");
      }
    //initialize list_maxes
    for (int init = 0; init < atom->nlocal; init++) list_maxes[init] = maxneigh;
  }

  //grow neighbor list arrays
  if (neigh_allocated&&atom->nlocal > max_atom_count) {
    atom->list_container = list_container = (int **) memory->srealloc(list_container, atom->nlocal*sizeof(int *), "NPair CAC:list_container" );
    memory->grow(list_maxes, atom->nlocal, "NPair CAC:list_maxes");
      for (int init = max_atom_count; init < atom->nlocal; init++) {
        list_container[init] = memory->create(list_container[init],maxneigh, "NPair CAC:cell_indexes");
      }
    //initialize list_maxes
    for (int init = max_atom_count; init < atom->nlocal; init++) list_maxes[init] = maxneigh;
  }
  
  neigh_allocated = 1;
  if(atom->nlocal > max_atom_count)
  max_atom_count = atom->nlocal;
}

//allocate surface counts array
void NPairCAC::allocate_local_arrays() {
  if(!ghost_quad_flag){
    atom->surface_counts = memory->grow(surface_counts, atom->nlocal, 3, "NPairCAC:surface_counts");
    atom->interior_scales = memory->grow(interior_scales, atom->nlocal, 3, "NPairCAC:interior_scales");
    atom->quadrature_counts = memory->grow(quadrature_counts, atom->nlocal, "NPair CAC:quadrature_counts");
    atom->e2quad_index = memory->grow(e2quad_index, atom->nlocal, "NPair CAC:e2quad_index");
    atom->neighbor_weights = memory->grow(neighbor_weights, atom->nlocal,3, "NPair CAC:neighbor_weights");
    nmax = atom->nlocal;
  }
  else{
    atom->surface_counts = memory->grow(surface_counts, atom->nlocal + atom->nghost, 3, "NPairCAC:surface_counts");
    atom->interior_scales = memory->grow(interior_scales, atom->nlocal + atom->nghost, 3, "NPairCAC:interior_scales");
    atom->quadrature_counts = memory->grow(quadrature_counts, atom->nlocal + atom->nghost, "NPair CAC:quadrature_counts");
    atom->e2quad_index = memory->grow(e2quad_index, atom->nlocal + atom->nghost, "NPair CAC:e2quad_index");
    atom->neighbor_weights = memory->grow(neighbor_weights, atom->nlocal + atom->nghost,3, "NPair CAC:neighbor_weights");
    nmax = atom->nlocal + atom->nghost;
  }
}

/* ----------------------------------------------------------------------
   grow quadrature virtual neighbor arrays
------------------------------------------------------------------------- */

void NPairCAC::allocate_quad_neigh_list() {
  int *element_type = atom->element_type;
  int quad = quadrature_node_count;
  int quad_count, maxneigh_quad_inner, maxneigh_quad_outer;
  int sum_limit;
  int *poly_count = atom->poly_count;

  sum_limit = atom->nlocal;
  quad_count = 0;
  maxneigh_quad_inner = atom->max_neigh_inner_init;
  maxneigh_quad_outer = atom->max_neigh_outer_init;
  
  //count the total number of quadrature points
  for(int sum = 0; sum < sum_limit; sum++)
  quad_count += poly_count[sum]*quadrature_counts[sum];
  // initialize quadrature point neighbor list arrays
  if(!quad_allocated){
    atom->inner_quad_lists_ucell = inner_quad_lists_ucell = (double ***) memory->smalloc(sizeof(double **)*quad_count, "NPair CAC:inner_quad_lists_ucell");
    atom->inner_quad_lists_index = inner_quad_lists_index = (int ***) memory->smalloc(sizeof(int **)*quad_count, "NPair CAC:inner_quad_lists_index");
    atom->inner_quad_lists_counts = memory->create(inner_quad_lists_counts, quad_count, "NPair CAC:inner_quad_lists_counts");
    atom->inner_quad_neigh_maxes = memory->create(inner_quad_neigh_maxes, quad_count, "NPair CAC:inner_quad_neigh_maxes");
      for (int neigh_loop = 0; neigh_loop < quad_count; neigh_loop++) {
        memory->create(inner_quad_lists_ucell[neigh_loop], maxneigh_quad_inner, 3, "NPair CAC:inner_quad_lists_ucell");
        if(sector_flag)
        memory->create(inner_quad_lists_index[neigh_loop], maxneigh_quad_inner, 3, "NPair CAC:inner_quad_lists_index");
        else
        memory->create(inner_quad_lists_index[neigh_loop], maxneigh_quad_inner, 2, "NPair CAC:inner_quad_lists_index");
      }
      if (outer_neigh_flag) {
        atom->outer_quad_lists_ucell = outer_quad_lists_ucell = (double ***) memory->smalloc(sizeof(double **)*quad_count, "NPair CAC:outer_quad_lists_ucell");
        atom->outer_quad_lists_index = outer_quad_lists_index = (int ***) memory->smalloc(sizeof(int **)*quad_count, "NPair CAC:outer_quad_lists_index");
        atom->outer_quad_lists_counts = memory->create(outer_quad_lists_counts, quad_count, "NPair CAC:outer_quad_lists_counts");
        atom->outer_quad_neigh_maxes = memory->create(outer_quad_neigh_maxes, quad_count, "NPair CAC:outer_quad_neigh_maxes");
        for (int neigh_loop = 0; neigh_loop < quad_count; neigh_loop++) {
          memory->create(outer_quad_lists_ucell[neigh_loop], maxneigh_quad_outer, 3, "NPair CAC:outer_quad_lists_ucell");
          if(sector_flag)
          memory->create(outer_quad_lists_index[neigh_loop], maxneigh_quad_outer, 3, "NPair CAC:outer_quad_lists_index");
          else
          memory->create(outer_quad_lists_index[neigh_loop], maxneigh_quad_outer, 2, "NPair CAC:outer_quad_lists_index");
        }
      }
    //initialize count maxes
    for (int init = 0; init < quad_count; init++) {
      inner_quad_neigh_maxes[init] = maxneigh_quad_inner;
      if (outer_neigh_flag) {
        outer_quad_neigh_maxes[init] = maxneigh_quad_outer;
      }
    }      
  }

  //grow quadrature point neighbor list arrays if needed
  if (quad_allocated&&quad_count>quad_count_max) {
    atom->inner_quad_lists_ucell = inner_quad_lists_ucell = (double ***) memory->srealloc(inner_quad_lists_ucell, sizeof(double **)*quad_count, "Pair CAC:inner_quad_lists_ucell");
    atom->inner_quad_lists_index = inner_quad_lists_index = (int ***) memory->srealloc(inner_quad_lists_index, sizeof(int **)*quad_count, "Pair CAC:inner_quad_lists_index");
    atom->inner_quad_lists_counts = memory->grow(inner_quad_lists_counts, quad_count, "Pair CAC:inner_quad_lists_counts");
    atom->inner_quad_neigh_maxes = memory->grow(inner_quad_neigh_maxes, quad_count, "Pair CAC:inner_quad_neigh_maxes");
    for(int neigh_loop = quad_count_max; neigh_loop < quad_count; neigh_loop++){
      memory->create(inner_quad_lists_ucell[neigh_loop], maxneigh_quad_inner, 3, "Pair CAC:inner_quad_lists_ucell");
      if(sector_flag)
      memory->create(inner_quad_lists_index[neigh_loop], maxneigh_quad_inner, 3, "Pair CAC:inner_quad_lists_index");
      else
      memory->create(inner_quad_lists_index[neigh_loop], maxneigh_quad_inner, 2, "Pair CAC:inner_quad_lists_index");
    }
    if (outer_neigh_flag) {
      atom->outer_quad_lists_ucell = outer_quad_lists_ucell = (double ***) memory->srealloc(outer_quad_lists_ucell, sizeof(double **)*quad_count, "NPair CAC:outer_quad_lists_ucell");
      atom->outer_quad_lists_index = outer_quad_lists_index = (int ***) memory->srealloc(outer_quad_lists_index, sizeof(int **)*quad_count, "NPair CAC:outer_quad_lists_index");
      atom->outer_quad_lists_counts = memory->grow(outer_quad_lists_counts, quad_count, "NPair CAC:outer_quad_lists_counts");
      atom->outer_quad_neigh_maxes = memory->grow(outer_quad_neigh_maxes, quad_count, "NPair CAC:outer_quad_neigh_maxes");
      for(int neigh_loop = quad_count_max; neigh_loop < quad_count; neigh_loop++){
        memory->create(outer_quad_lists_ucell[neigh_loop], maxneigh_quad_outer, 3, "NPair CAC:outer_quad_lists_ucell");
        if(sector_flag)
        memory->create(outer_quad_lists_index[neigh_loop], maxneigh_quad_outer, 3, "NPair CAC:outer_quad_lists_index");
        else
        memory->create(outer_quad_lists_index[neigh_loop], maxneigh_quad_outer, 2, "NPair CAC:outer_quad_lists_index");
      }
    }

  //initialize count maxes
  for (int init = quad_count_max; init < quad_count; init++) {
    inner_quad_neigh_maxes[init] = maxneigh_quad_inner;
    if (outer_neigh_flag) {
      outer_quad_neigh_maxes[init] = maxneigh_quad_outer;
    }
  }    
  } 
  
  //initialize counts to zero
  for (int init = 0; init < quad_count; init++) {
    inner_quad_lists_counts[init] = 0;
    if (outer_neigh_flag) {
      outer_quad_lists_counts[init] = 0;
    }
  }
  quad_allocated = 1;
  if(quad_count > quad_count_max)
  quad_count_max = quad_count;
}

/* ----------------------------------------------------------------------
  associate a quadrature point sector (region closest to quadrature point
  in mapped undeformed space) with a coordinate; designed for Q8 only now
------------------------------------------------------------------------- */

int NPairCAC::quad_sector_select(double s, double t, double w, int eindex,int pindex){
int *element_scale = atom->element_scale[eindex];
int element_type = atom->element_type[eindex];
int quad_index = 0;
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
quad_index = quad*quad*(signw+1)/2+quad*(signt+1)/2+(signs+1)/2;
return quad_index;
}
//if s surface point
if(surfs!=0&&surft==0&&surfw==0){
  if(surfs==-1) quad_index = quad*quad*quad + quad*quad*s1+quad*(signt+1)/2+(signw+1)/2;
  if(surfs==1) quad_index = quad*quad*quad + quad*quad*n1+quad*quad*s1+quad*(signt+1)/2+(signw+1)/2;
  return quad_index;
}
//if t surface point
if(surfs==0&&surft!=0&&surfw==0){
  if(surft==-1) quad_index = quad*quad*quad + 2*quad*quad*n1+quad*quad*s2+quad*(signs+1)/2+(signw+1)/2;
  if(surft==1) quad_index = quad*quad*quad + 2*quad*quad*n1+quad*quad*n2+quad*quad*s2+quad*(signs+1)/2+(signw+1)/2;
  return quad_index;
}
//if w surface point
if(surfs==0&&surft==0&&surfw!=0){
  if(surfw==-1) quad_index = quad*quad*quad + 2*quad*quad*(n1+n2)+quad*quad*s3+quad*(signs+1)/2+(signt+1)/2;
  if(surfw==1) quad_index = quad*quad*quad + 2*quad*quad*(n1+n2)+quad*quad*n3+quad*quad*s3+quad*(signs+1)/2+(signt+1)/2;
  return quad_index;
}
//if w oriented edge point
if(surfs!=0&&surft!=0&&surfw==0){
quad_index = quad*quad*quad + 2*quad*quad*(n1+n2+n3);
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
quad_index = quad*quad*quad + 2*quad*quad*(n1+n2+n3)+4*quad*n1*n2;
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
quad_index = quad*quad*quad + 2*quad*quad*(n1+n2+n3)+4*quad*n1*n2+4*quad*n2*n3;
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
quad_index = quad*quad*quad + 2*quad*quad*(n1+n2+n3)+4*quad*(n1*n2+n2*n3+n1*n3);
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
return quad_index;
}

/* ----------------------------------------------------------------------
   grow quadrature data arrays
------------------------------------------------------------------------- */

void NPairCAC::grow_quad_data(){
  quadrature_point_max += 1000*atom->maxpoly;	
  atom->quadrature_point_data = memory->grow(quadrature_point_data,quadrature_point_max,7,"pairCAC:quadrature_point_data");
  atom->quadrature_point_max = quadrature_point_max;
}

/* ----------------------------------------------------------------------
  grow custom quadrature int data arrays; requires quad_array to be NULL initialized
  and size must remain constant.
------------------------------------------------------------------------- */

void NPairCAC::allocate_quad_int(int **&int_quad_array, int size){
  if(quadrature_poly_count > quadrature_poly_max){
  quadrature_poly_max = quadrature_poly_count;
  int_quad_array = memory->grow(int_quad_array,quadrature_poly_max,size,"pairCAC:int_quad_array");
  }
}

/* ----------------------------------------------------------------------
  grow custom quadrature double data arrays; requires quad_array to be NULL initialized
  and size must remain constant.
------------------------------------------------------------------------- */

void NPairCAC::allocate_quad_double(double **&double_quad_array, int size){
  if(quadrature_poly_count > quadrature_poly_max){
  quadrature_poly_max = quadrature_poly_count;
  double_quad_array = memory->grow(double_quad_array,quadrature_poly_max,size,"pairCAC:double_quad_array");
  }
}

/* ----------------------------------------------------------------------
   pack the buffer communicating the interface flags that decide the quadrature 
   points for ghost elements
------------------------------------------------------------------------- */

int NPairCAC::pack_forward_comm(int n, int *list, double *buf,
                               int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = interface_flags[j];
  }
  return m;
}

/* ----------------------------------------------------------------------
   pack the buffer communicating the interface flags that decide the quadrature 
   points for ghost elements
------------------------------------------------------------------------- */

void NPairCAC::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    interface_flags[i] = buf[m++];
  }
}

//memory usage due to quadrature point list memory structure
bigint NPairCAC::memory_usage()
{
  int i;
  bigint bytes_used = 0;
  if (neigh_allocated) {
    for (i = 0; i < max_atom_count; i++) {
      bytes_used +=memory->usage(list_container[i],list_maxes[i]);
    }
  }
  
  if(quad_allocated){
    for (i = 0; i < quad_count_max; i++) {
      if(sector_flag){
      bytes_used +=memory->usage(inner_quad_lists_index[i],inner_quad_neigh_maxes[i],3);
      if(outer_neigh_flag)
      bytes_used +=memory->usage(outer_quad_lists_index[i],outer_quad_neigh_maxes[i],3);
      }
      else{
      bytes_used +=memory->usage(inner_quad_lists_index[i],inner_quad_neigh_maxes[i],2);
      if(outer_neigh_flag)
      bytes_used +=memory->usage(outer_quad_lists_index[i],outer_quad_neigh_maxes[i],2);
      }
      bytes_used +=memory->usage(inner_quad_lists_ucell[i],inner_quad_neigh_maxes[i],3);
      if(outer_neigh_flag)
      bytes_used +=memory->usage(outer_quad_lists_ucell[i],outer_quad_neigh_maxes[i],3);
    }
  }
  
  return bytes_used;
}