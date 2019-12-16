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
#include "atom_vec.h"
#include "molecule.h"
#include "domain.h"
#include "my_page.h"
#include "error.h"
#include "memory.h"

#define MAXNEIGH  10
#define EXPAND 50
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NPairCAC::NPairCAC(LAMMPS *lmp) : NPair(lmp) {
  nmax=0;
  maxneigh_quad = MAXNEIGH;
  max_expansion_count = 0;
  bad_bin_flag=0;
  
  surface_counts = NULL;
  interior_scales = NULL;
  old_atom_etype = NULL;
  current_element_quad_points=NULL;
  quad_allocated = 0;
  scan_flags = NULL;
  max_quad_alloc=0;
  surface_counts_max[0] = 1;
  surface_counts_max[1] = 1;
  surface_counts_max[2] = 1;
  surface_counts_max_old[0] = 1;
  surface_counts_max_old[1] = 1;
  surface_counts_max_old[2] = 1;
  memory->create(current_element_quad_points,64,3,"NPair CAC:current_element_quad_points");
  
}

/* ---------------------------------------------------------------------- */

NPairCAC::~NPairCAC() 
{

  if (quad_allocated) {
    for (int init = 0; init < old_atom_count; init++) {
      if (old_atom_etype[init] == 0) {
        memory->destroy(quad_list_container[init][0]);
        memory->sfree(quad_list_container[init]);	
    }
      else {
        for (int neigh_loop = 0; neigh_loop < old_quad_count; neigh_loop++) {
          memory->destroy(quad_list_container[init][neigh_loop]);
    }
    memory->sfree(quad_list_container[init]);		
    }
  }
  memory->sfree(quad_list_container);	
  memory->destroy(current_element_quad_points);
  }
  memory->destroy(surface_counts);
  memory->destroy(interior_scales);
}

/* ----------------------------------------------------------------------
   binned neighbor list construction for all neighbors
   Neighboring is based on the overlap between element bounding boxes + cutoff bound
   and points (quadrature points or pure atoms)
------------------------------------------------------------------------- */

void NPairCAC::build(NeighList *list)
{
  int i, j, k, n, itype, jtype, ibin, iatom;
  double delx, dely, delz, rsq;
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

  cutneighmax = neighbor->cutneighmax;
  
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  int *element_type = atom->element_type;
  int *poly_count = atom->poly_count;
  element_scale = atom->element_scale;
  int current_element_type;
  int neighbor_element_type;

  quadrature_init(2);
  memory->grow(scan_flags,atom->nlocal + atom->nghost,"scan_flags");
  //initialize scan flags to 0
  for (int init=0; init < atom->nlocal + atom->nghost; init++){
    scan_flags[init]=0;
  }

  //compute maximum surface counts for quadrature neighbor list allocation
  

      // initialize or grow surface counts array for quadrature scheme
      // along with interior scaling for the quadrature domain
      if (atom->nlocal  > nmax) {
        allocate_surface_counts();
      }
      surface_counts_max_old[0] = surface_counts_max[0];
      surface_counts_max_old[1] = surface_counts_max[1];
      surface_counts_max_old[2] = surface_counts_max[2];
      for (i = 0; i < atom->nlocal; i++) {
      current_element_scale[0] = element_scale[i][0];
      current_element_scale[1] = element_scale[i][1];
      current_element_scale[2] = element_scale[i][2];	
        if (element_type[i] != 0) {
          
          for (current_poly_counter = 0; current_poly_counter < poly_count[i]; current_poly_counter++) {
            int poly_surface_count[3];
                  current_nodal_positions = nodal_positions[i][current_poly_counter];
            compute_surface_depths(interior_scale[0], interior_scale[1], interior_scale[2],
              poly_surface_count[0], poly_surface_count[1], poly_surface_count[2], 1);
            if(current_poly_counter==0){
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

      // initialize or grow memory for the neighbor list of virtual atoms at quadrature points

      if (atom->nlocal)
      allocate_quad_neigh_list(surface_counts_max[0], surface_counts_max[1], surface_counts_max[2], quadrature_node_count);
      //firstneigh[0] = neighptr;
      if(quadrature_point_count>list->maxatom){
      memory->destroy(list->ilist);
      memory->destroy(list->numneigh);
      memory->create(list->ilist,quadrature_point_count,"NPairCAC:ilist");
      memory->create(list->numneigh,quadrature_point_count,"NPairCAC:numneigh");
      memory->sfree(list->firstneigh);
      list->firstneigh=(int **) memory->smalloc(quadrature_point_count*sizeof(int *),
                                        "NPairCAC:firstneigh");
     
      max_quad_alloc=quadrature_point_count;
      ilist = list->ilist;
      numneigh = list->numneigh;
      firstneigh = list->firstneigh;
                      
      }

  int inum = 0;
  int qi=0;
  int quadrature_count;
 //loop over elements
  for (i = 0; i < nlocal; i++) {
  itype = type[i];
  current_element_type = element_type[i];
  current_element_scale[0] = element_scale[i][0];
  current_element_scale[1] = element_scale[i][1];
  current_element_scale[2] = element_scale[i][2];
  //find the current quadrature points of this element
  if(current_element_type!=0){
    quadrature_count=compute_quad_points(i);
  }
  else{
    quadrature_count=1;
    current_element_quad_points[0][0]=x[i][0];
    current_element_quad_points[0][1]=x[i][1];
    current_element_quad_points[0][2]=x[i][2];
  }
    
// loop over all atoms in surrounding bins in stencil including self
// skip i = j
//loop over quadrature points of elements, by convention an atom is one quadrature point
for (int iquad = 0; iquad < quadrature_count; iquad++) {
  n = 0;
  ibin = quad2bin[qi];
  expansion_count=0;
  int neigh_index;
  
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
      neigh_index=n;
      if (neigh_index == maxneigh_quad + expansion_count*EXPAND) {
          expansion_count += 1;
            memory->grow(quad_list_container[i][iquad], 
            maxneigh_quad + expansion_count*EXPAND, "NPair CAC:cell indexes expand");

      }
      
      current_quad_point[0]=current_element_quad_points[iquad][0];
      current_quad_point[1]=current_element_quad_points[iquad][1];
      current_quad_point[2]=current_element_quad_points[iquad][2];
      neighptr = quad_list_container[i][iquad];
      
      if (neighbor_element_type != 0) {
      if(CAC_decide_quad2element(j)){
       if(scan_flags[j]) continue;
       else scan_flags[j]=1;
       neighptr[n++] = j;
      }
      }
      else if(neighbor_element_type == 0){
      delx = current_quad_point[0] - x[j][0];
      dely = current_quad_point[1] - x[j][1];
      delz = current_quad_point[2] - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
        if (rsq <= cutneighmax*cutneighmax)
       neighptr[n++] = j;
      }

  }
      
  }
  numneigh[qi] = n;
  ilist[qi] = i;
  firstneigh[qi] = neighptr;
  qi++;
  
  //reset scan flags for next quadrture point searching for its neighbors
  for(int reset=0; reset<n; reset++){
    scan_flags[neighptr[reset]]=0;
  }

}

  list->inum = inum;
  list->gnum = 0;
  
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

////////////////////////////////////////////////////////

void NPairCAC::compute_surface_depths(double &scalex, double &scaley, double &scalez,
  int &countx, int &county, int &countz, int flag) {
  int poly = current_poly_counter;
  double unit_cell_mapped[3];
  double rcut;
  rcut = cutneighmax - neighbor->skin; 
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
  if (ds_surf > 1) {
    ds_surf = 1;


  }
  if (dt_surf > 1) {

    dt_surf = 1;


  }
  if (dw_surf > 1) {

    dw_surf = 1;

  }

  if (atom->one_layer_flag) {
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


/////////////////////////////////


void NPairCAC::quadrature_init(int quadrature_rank) {

  if (quadrature_rank == 1) {
    quadrature_node_count = 1;
    
    memory->create(quadrature_abcissae, quadrature_node_count, "pairCAC:quadrature_abcissae");

    quadrature_abcissae[0] = 0;
  }
  if (quadrature_rank == 2) {

    quadrature_node_count = 2;

    memory->create(quadrature_abcissae, quadrature_node_count, "pairCAC:quadrature_abcissae");

    quadrature_abcissae[0] = -0.5773502691896258;
    quadrature_abcissae[1] = 0.5773502691896258;

  }

  if (quadrature_rank == 3)
  {

  }
  if (quadrature_rank == 4)
  {

  }
  if (quadrature_rank == 5)
  {

  }


}


//compute and store quad points for the element referred to by element_index

int NPairCAC::compute_quad_points(int element_index){

  double unit_cell_mapped[3];
  double interior_scale[3];
  int surface_count[3];
  int nodes_per_element;
  double s, t, w;
  s = t = w = 0;
  double quad_position[3];
  double ****nodal_positions = atom->nodal_positions;
  double shape_func;
  int *element_type = atom->element_type;
  int *poly_count = atom->poly_count;
  int *nodes_per_element_list = atom->nodes_per_element_list;
   nodes_per_element = nodes_per_element_list[element_type[element_index]];
  unit_cell_mapped[0] = 2 / double(element_scale[element_index][0]);
  unit_cell_mapped[1] = 2 / double(element_scale[element_index][1]);
  unit_cell_mapped[2] = 2 / double(element_scale[element_index][2]);
  int current_poly_count = poly_count[element_index];
  int quadrature_counter=0;
  int signs, signt, signw;

  //compute quadrature point positions to test for neighboring
  int sign[2];
  sign[0] = -1;
  sign[1] = 1;
  
  surface_count[0]=surface_counts[element_index][0];
  surface_count[1]=surface_counts[element_index][1];
  surface_count[2]=surface_counts[element_index][2];
  interior_scale[0]= interior_scales[element_index][0];
  interior_scale[1]= interior_scales[element_index][1];
  interior_scale[2]= interior_scales[element_index][2];
  for (int poly_counter = 0; poly_counter < current_poly_count; poly_counter++) {
    current_nodal_positions = nodal_positions[element_index][poly_counter];
    //interior contributions
    if(element_scale[element_index][0]>2&&element_scale[element_index][1]>2&&element_scale[element_index][2]>2){
    for (int i = 0; i < quadrature_node_count; i++) {
      for (int j = 0; j < quadrature_node_count; j++) {
        for (int k = 0; k < quadrature_node_count; k++) {

          s = interior_scale[0] * quadrature_abcissae[i];
          t = interior_scale[1] * quadrature_abcissae[j];
          w = interior_scale[2] * quadrature_abcissae[k];
          signs=signt=signw=1;
          if(s<0) signs=-1;
          if(t<0) signt=-1;
          if(w<0) signw=-1;
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


          quad_position[0] = 0;
          quad_position[1] = 0;
          quad_position[2] = 0;
          for (int kkk = 0; kkk < nodes_per_element; kkk++) {
            shape_func = shape_function(s, t, w, 2, kkk + 1);
            quad_position[0] += current_nodal_positions[kkk][0] * shape_func;
            quad_position[1] += current_nodal_positions[kkk][1] * shape_func;
            quad_position[2] += current_nodal_positions[kkk][2] * shape_func;
          }

         current_element_quad_points[quadrature_counter][0]=quad_position[0];
         current_element_quad_points[quadrature_counter][1]=quad_position[1];
         current_element_quad_points[quadrature_counter][2]=quad_position[2];
                 quadrature_counter+=1;

        }
      }
    }
    }
    // s axis surface contributions
  if(element_scale[element_index][1]>2&&element_scale[element_index][2]>2){
  for (int sc = 0; sc < 2; sc++) {
    for (int i = 0; i < surface_count[0]; i++) {
      for (int j = 0; j < quadrature_node_count; j++) {
        for (int k = 0; k < quadrature_node_count; k++) {
          
            s = sign[sc] - i*unit_cell_mapped[0] * sign[sc];
            s = s - 0.5*unit_cell_mapped[0] * sign[sc];
            t = interior_scale[1] * quadrature_abcissae[j];
            w = interior_scale[2] * quadrature_abcissae[k];
            signs=signt=signw=1;
            if(t<0) signt=-1;
            if(w<0) signw=-1;
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

            quad_position[0] = 0;
            quad_position[1] = 0;
            quad_position[2] = 0;
            for (int kkk = 0; kkk < nodes_per_element; kkk++) {
              shape_func = shape_function(s, t, w, 2, kkk + 1);
              quad_position[0] += current_nodal_positions[kkk][0] * shape_func;
              quad_position[1] += current_nodal_positions[kkk][1] * shape_func;
              quad_position[2] += current_nodal_positions[kkk][2] * shape_func;
            }

         current_element_quad_points[quadrature_counter][0]=quad_position[0];
         current_element_quad_points[quadrature_counter][1]=quad_position[1];
         current_element_quad_points[quadrature_counter][2]=quad_position[2];
         quadrature_counter+=1;
          }
        }
      }
    }
  }

  // t axis contributions
  if(element_scale[element_index][0]>2&&element_scale[element_index][2]>2){
  for (int sc = 0; sc < 2; sc++) {
    for (int i = 0; i < surface_count[1]; i++) {
      for (int j = 0; j < quadrature_node_count; j++) {
        for (int k = 0; k < quadrature_node_count; k++) {
            s = interior_scale[0] * quadrature_abcissae[j];
            t = sign[sc] - i*unit_cell_mapped[1] * sign[sc];
            t = t - 0.5*unit_cell_mapped[1] * sign[sc];
            w = interior_scale[2] * quadrature_abcissae[k];
            signs=signt=signw=1;
            if(s<0) signs=-1;
            if(w<0) signw=-1;
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

            quad_position[0] = 0;
            quad_position[1] = 0;
            quad_position[2] = 0;
            for (int kkk = 0; kkk < nodes_per_element; kkk++) {
              shape_func = shape_function(s, t, w, 2, kkk + 1);
              quad_position[0] += current_nodal_positions[kkk][0] * shape_func;
              quad_position[1] += current_nodal_positions[kkk][1] * shape_func;
              quad_position[2] += current_nodal_positions[kkk][2] * shape_func;
            }

         current_element_quad_points[quadrature_counter][0]=quad_position[0];
         current_element_quad_points[quadrature_counter][1]=quad_position[1];
         current_element_quad_points[quadrature_counter][2]=quad_position[2];
         quadrature_counter+=1;

          }
        }
      }
    }
  }

  //w axis surface contributions
  if(element_scale[element_index][0]>2&&element_scale[element_index][1]>2){
  for (int sc = 0; sc < 2; sc++) {
    for (int i = 0; i < surface_count[2]; i++) {
      for (int j = 0; j < quadrature_node_count; j++) {
        for (int k = 0; k < quadrature_node_count; k++) {
            s = interior_scale[0] * quadrature_abcissae[j];
            t = interior_scale[1] * quadrature_abcissae[k];
            signs=signt=signw=1;
            if(s<0) signs=-1;
            if(t<0) signt=-1;
            s = unit_cell_mapped[0] * (int((s+signs) / unit_cell_mapped[0]))-signs;
            t = unit_cell_mapped[1] * (int((t+signt) / unit_cell_mapped[1]))-signt;
            w = sign[sc] - i*unit_cell_mapped[2] * sign[sc];
            w = w - 0.5*unit_cell_mapped[2] * sign[sc];

            if (quadrature_abcissae[j] < 0)
              s = s - 0.5*unit_cell_mapped[0];
            else
              s = s + 0.5*unit_cell_mapped[0];

            if (quadrature_abcissae[k] < 0)
              t = t - 0.5*unit_cell_mapped[1];
            else
              t = t + 0.5*unit_cell_mapped[1];


            quad_position[0] = 0;
            quad_position[1] = 0;
            quad_position[2] = 0;
            for (int kkk = 0; kkk < nodes_per_element; kkk++) {
              shape_func = shape_function(s, t, w, 2, kkk + 1);
              quad_position[0] += current_nodal_positions[kkk][0] * shape_func;
              quad_position[1] += current_nodal_positions[kkk][1] * shape_func;
              quad_position[2] += current_nodal_positions[kkk][2] * shape_func;
            }

         current_element_quad_points[quadrature_counter][0]=quad_position[0];
         current_element_quad_points[quadrature_counter][1]=quad_position[1];
         current_element_quad_points[quadrature_counter][2]=quad_position[2];
         quadrature_counter+=1;

          }
        }
      }
    }
  }
  int surface_countx;
  int surface_county;

  //compute edge contributions

  for (int sc = 0; sc < 12; sc++) {
    if (sc == 0 || sc == 1 || sc == 2 || sc == 3) {
      if(element_scale[element_index][2]<=2) continue;
      surface_countx = surface_count[0];
      surface_county = surface_count[1];
    }
    else if (sc == 4 || sc == 5 || sc == 6 || sc == 7) {
      if(element_scale[element_index][0]<=2) continue;
      surface_countx = surface_count[1];
      surface_county = surface_count[2];
    }
    else if (sc == 8 || sc == 9 || sc == 10 || sc == 11) {
      if(element_scale[element_index][1]<=2) continue;
      surface_countx = surface_count[0];
      surface_county = surface_count[2];
    }

    
    for (int i = 0; i < surface_countx; i++) {
      for (int j = 0; j < surface_county; j++) {
        
        for (int k = 0; k < quadrature_node_count; k++) {
          
            if (sc == 0) {

              s = -1 + (i + 0.5)*unit_cell_mapped[0];
              t = -1 + (j + 0.5)*unit_cell_mapped[1];
              w = interior_scale[2] * quadrature_abcissae[k];
              signs=signt=signw=1;
                  if(w<0) signw=-1;
              w = unit_cell_mapped[2] * (int((w+signw) / unit_cell_mapped[2]))-signw;
              if (quadrature_abcissae[k] < 0)
                w = w - 0.5*unit_cell_mapped[2];
              else
                w = w + 0.5*unit_cell_mapped[2];
            }
            else if (sc == 1) {
              s = 1 - (i + 0.5)*unit_cell_mapped[0];
              t = -1 + (j + 0.5)*unit_cell_mapped[1];
              w = interior_scale[2] * quadrature_abcissae[k];
              signs=signt=signw=1;
                  if(w<0) signw=-1;
              w = unit_cell_mapped[2] * (int((w+signw) / unit_cell_mapped[2]))-signw;
              if (quadrature_abcissae[k] < 0)
                w = w - 0.5*unit_cell_mapped[2];
              else
                w = w + 0.5*unit_cell_mapped[2];
            }
            else if (sc == 2) {
              s = -1 + (i + 0.5)*unit_cell_mapped[0];
              t = 1 - (j + 0.5)*unit_cell_mapped[1];
              w = interior_scale[2] * quadrature_abcissae[k];
              signs=signt=signw=1;
                  if(w<0) signw=-1;
              w = unit_cell_mapped[2] * (int((w+signw) / unit_cell_mapped[2]))-signw;
              if (quadrature_abcissae[k] < 0)
                w = w - 0.5*unit_cell_mapped[2];
              else
                w = w + 0.5*unit_cell_mapped[2];
            }
            else if (sc == 3) {
              s = 1 - (i + 0.5)*unit_cell_mapped[0];
              t = 1 - (j + 0.5)*unit_cell_mapped[1];
              w = interior_scale[2] * quadrature_abcissae[k];
              signs=signt=signw=1;
                  if(w<0) signw=-1;
              w = unit_cell_mapped[2] * (int((w+signw) / unit_cell_mapped[2]))-signw;
              if (quadrature_abcissae[k] < 0)
                w = w - 0.5*unit_cell_mapped[2];
              else
                w = w + 0.5*unit_cell_mapped[2];
            }
            else if (sc == 4) {
              s = interior_scale[0] * quadrature_abcissae[k];
              signs=signt=signw=1;
                  if(s<0) signs=-1;
              s = unit_cell_mapped[0] * (int((s+signs) / unit_cell_mapped[0]))-signs;
              t = -1 + (i + 0.5)*unit_cell_mapped[1];
              w = -1 + (j + 0.5)*unit_cell_mapped[2];
              if (quadrature_abcissae[k] < 0)
                s = s - 0.5*unit_cell_mapped[0];
              else
                s = s + 0.5*unit_cell_mapped[0];

            }
            else if (sc == 5) {
              s = interior_scale[0] * quadrature_abcissae[k];
              signs=signt=signw=1;
                  if(s<0) signs=-1;
              s = unit_cell_mapped[0] * (int((s+signs) / unit_cell_mapped[0]))-signs;
              t = 1 - (i + 0.5)*unit_cell_mapped[1];
              w = -1 + (j + 0.5)*unit_cell_mapped[2];
              if (quadrature_abcissae[k] < 0)
                s = s - 0.5*unit_cell_mapped[0];
              else
                s = s + 0.5*unit_cell_mapped[0];
            }
            else if (sc == 6) {
              s = interior_scale[0] * quadrature_abcissae[k];
              signs=signt=signw=1;
                  if(s<0) signs=-1;
              s = unit_cell_mapped[0] * (int((s+signs) / unit_cell_mapped[0]))-signs;
              t = -1 + (i + 0.5)*unit_cell_mapped[1];
              w = 1 - (j + 0.5)*unit_cell_mapped[2];
              if (quadrature_abcissae[k] < 0)
                s = s - 0.5*unit_cell_mapped[0];
              else
                s = s + 0.5*unit_cell_mapped[0];
            }
            else if (sc == 7) {
              s = interior_scale[0] * quadrature_abcissae[k];
              signs=signt=signw=1;
                  if(s<0) signs=-1;
              s = unit_cell_mapped[0] * (int((s+signs) / unit_cell_mapped[0]))-signs;
              t = 1 - (i + 0.5)*unit_cell_mapped[1];
              w = 1 - (j + 0.5)*unit_cell_mapped[2];
              if (quadrature_abcissae[k] < 0)
                s = s - 0.5*unit_cell_mapped[0];
              else
                s = s + 0.5*unit_cell_mapped[0];
            }
            else if (sc == 8) {
              s = -1 + (i + 0.5)*unit_cell_mapped[0];
              t = interior_scale[1] * quadrature_abcissae[k];
              signs=signt=signw=1;
                  if(t<0) signt=-1;
              t = unit_cell_mapped[1] * (int((t+signt) / unit_cell_mapped[1]))-signt;
              w = -1 + (j + 0.5)*unit_cell_mapped[2];
              if (quadrature_abcissae[k] < 0)
                t = t - 0.5*unit_cell_mapped[1];
              else
                t = t + 0.5*unit_cell_mapped[1];

            }
            else if (sc == 9) {
              s = 1 - (i + 0.5)*unit_cell_mapped[0];
              t = interior_scale[1] * quadrature_abcissae[k];
              signs=signt=signw=1;
                  if(t<0) signt=-1;
              t = unit_cell_mapped[1] * (int((t+signt) / unit_cell_mapped[1]))-signt;
              w = -1 + (j + 0.5)*unit_cell_mapped[2];
              if (quadrature_abcissae[k] < 0)
                t = t - 0.5*unit_cell_mapped[1];
              else
                t = t + 0.5*unit_cell_mapped[1];
            }
            else if (sc == 10) {
              s = -1 + (i + 0.5)*unit_cell_mapped[0];
              t = interior_scale[1] * quadrature_abcissae[k];
              signs=signt=signw=1;
                  if(t<0) signt=-1;
              t = unit_cell_mapped[1] * (int((t+signt) / unit_cell_mapped[1]))-signt;
              w = 1 - (j + 0.5)*unit_cell_mapped[2];
              if (quadrature_abcissae[k] < 0)
                t = t - 0.5*unit_cell_mapped[1];
              else
                t = t + 0.5*unit_cell_mapped[1];
            }
            else if (sc == 11) {
              s = 1 - (i + 0.5)*unit_cell_mapped[0];
              t = interior_scale[1] * quadrature_abcissae[k];
              signs=signt=signw=1;
                  if(t<0) signt=-1;
              t = unit_cell_mapped[1] * (int((t+signt) / unit_cell_mapped[1]))-signt;
              w = 1 - (j + 0.5)*unit_cell_mapped[2];
              if (quadrature_abcissae[k] < 0)
                t = t - 0.5*unit_cell_mapped[1];
              else
                t = t + 0.5*unit_cell_mapped[1];
            }

            quad_position[0] = 0;
            quad_position[1] = 0;
            quad_position[2] = 0;
            for (int kkk = 0; kkk < nodes_per_element; kkk++) {
              shape_func = shape_function(s, t, w, 2, kkk + 1);
              quad_position[0] += current_nodal_positions[kkk][0] * shape_func;
              quad_position[1] += current_nodal_positions[kkk][1] * shape_func;
              quad_position[2] += current_nodal_positions[kkk][2] * shape_func;
            }

         current_element_quad_points[quadrature_counter][0]=quad_position[0];
         current_element_quad_points[quadrature_counter][1]=quad_position[1];
         current_element_quad_points[quadrature_counter][2]=quad_position[2];
         quadrature_counter+=1;
          }

        }

      }
    }
  
  //compute corner contributions

  for (int sc = 0; sc < 8; sc++) {
    for (int i = 0; i < surface_count[0]; i++) {
      for (int j = 0; j < surface_count[1]; j++) {
        for (int k = 0; k < surface_count[2]; k++) {
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

            quad_position[0] = 0;
            quad_position[1] = 0;
            quad_position[2] = 0;
            for (int kkk = 0; kkk < nodes_per_element; kkk++) {
              shape_func = shape_function(s, t, w, 2, kkk + 1);
              quad_position[0] += current_nodal_positions[kkk][0] * shape_func;
              quad_position[1] += current_nodal_positions[kkk][1] * shape_func;
              quad_position[2] += current_nodal_positions[kkk][2] * shape_func;
            }

         current_element_quad_points[quadrature_counter][0]=quad_position[0];
         current_element_quad_points[quadrature_counter][1]=quad_position[1];
         current_element_quad_points[quadrature_counter][2]=quad_position[2];
         quadrature_counter+=1;

          }
        }

      }
    }
  }

  return quadrature_counter;

}

//decide if an atom or quadrature point is close enough to an element to consider for nonlocal quadrature calculation

int NPairCAC::CAC_decide_quad2element(int neighbor_element_index) {
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

void NPairCAC::allocate_quad_neigh_list(int n1,int n2,int n3,int quad) {
  int *element_type = atom->element_type;
  int *poly_count = atom->poly_count;
  int max_quad_count = quad*quad*quad + 2 * n1*quad*quad + 2 * n2*quad*quad +
    +2 * n3*quad*quad + 4 * n1*n2*quad + 4 * n3*n2*quad + 4 * n1*n3*quad
    + 8 * n1*n2*n3;
    //grow array used to store the set of quadrature points for each element at a time
    memory->grow(current_element_quad_points,max_quad_count*atom->maxpoly,3,"NPair CAC:current_element_quad_points");
 int c1,c2,c3;
 int current_quad_count;
 //number of atoms and quadrature points, i.e. atoms are counted as quadrature points.
  quadrature_point_count=0;
  
  maxneigh_quad += max_expansion_count*EXPAND;
  
  max_expansion_count = 0;
  
  // initialize quadrature point neighbor list vectors
  if (quad_allocated) {
    for (int init = 0; init < old_atom_count; init++) {

      if (old_atom_etype[init] == 0) {
        memory->destroy(quad_list_container[init][0]);
        memory->sfree(quad_list_container[init]); 
        
      }
      else {
        
        for (int neigh_loop = 0; neigh_loop < old_quad_count; neigh_loop++) {
          memory->destroy(quad_list_container[init][neigh_loop]);
        }
        memory->sfree(quad_list_container[init]);
      
      }
    }

    memory->sfree(quad_list_container);
  
  }
  
  
    quad_list_container= (int ***) memory->smalloc(atom->nlocal*sizeof(int **), "NPair CAC:quad_list_container");
    for (int init = 0; init < atom->nlocal; init++) {
      
      if (element_type[init] == 0) {
        quadrature_point_count+=1;
        quad_list_container[init]= (int **) memory->smalloc(sizeof(int *), "NPair CAC:list2ucell");
        quad_list_container[init][0]=memory->create(quad_list_container[init][0], maxneigh_quad, "NPair CAC:cell_indexes");
        
      }
      else {
        c1=surface_counts[init][0];
        c2=surface_counts[init][1];
        c3=surface_counts[init][2];
        current_quad_count = quad*quad*quad + 2 * c1*quad*quad + 2 * c2*quad*quad +
        +2 * c3*quad*quad + 4 * c1*c2*quad + 4 * c3*c2*quad + 4 * c1*c3*quad
         + 8 * c1*c2*c3;
        quadrature_point_count+=current_quad_count*poly_count[init];
        
        quad_list_container[init]= (int **) memory->smalloc(max_quad_count*atom->maxpoly*sizeof(int *), "NPair CAC:list2ucell");
        for (int neigh_loop = 0; neigh_loop < max_quad_count*atom->maxpoly; neigh_loop++) {

          quad_list_container[init][neigh_loop]=memory->create(quad_list_container[init][neigh_loop], maxneigh_quad, "NPair CAC:cell_indexes");


        }
      }
    }

  
  quad_allocated = 1;
  old_atom_count = atom->nlocal;
  old_quad_count = max_quad_count*atom->maxpoly;
  memory->grow(old_atom_etype, atom->nlocal, "NPair CAC:old_element_type_map");
  for (int init = 0; init < atom->nlocal; init++) {
    old_atom_etype[init]= element_type[init];
  }
}

//allocate surface counts array
void NPairCAC::allocate_surface_counts() {
  memory->grow(surface_counts, atom->nlocal , 3, "NPairCAC:surface_counts");
  memory->grow(interior_scales, atom->nlocal , 3, "NPairCAC:interior_scales");
  nmax = atom->nlocal;
}

//memory usage due to quadrature point list memory structure
bigint NPairCAC::memory_usage()
{
    int quad= quadrature_node_count;
    int n1=surface_counts_max[0];
    int n2=surface_counts_max[1];
    int n3=surface_counts_max[2];
    int max_quad_count = quad*quad*quad + 2 * n1*quad*quad + 2 * n2*quad*quad +
    +2 * n3*quad*quad + 4 * n1*n2*quad + 4 * n3*n2*quad + 4 * n1*n3*quad
    + 8 * n1*n2*n3;
  bigint bytes_used = 0;
    if (quad_allocated) {
    for (int init = 0; init < old_atom_count; init++) {

      if (old_atom_etype[init] == 0) {
        bytes_used +=memory->usage(quad_list_container[init][0], maxneigh_quad);
        
      }
      else {

        for (int neigh_loop = 0; neigh_loop < old_quad_count; neigh_loop++) {
          bytes_used +=memory->usage(quad_list_container[init][neigh_loop], maxneigh_quad);
        }
        
      
      }
    }

    bytes_used +=memory->usage(current_element_quad_points,max_quad_count);

  }
    
   
  return bytes_used;

}