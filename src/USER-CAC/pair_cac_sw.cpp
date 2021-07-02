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


#include "pair_cac_sw.h"
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "domain.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "utils.h"
#include "tokenizer.h"
#include "potential_file_reader.h"

//#include "math_extra.h"
#define PLANE_EPSILON  1e-6 //error tolerance for surface flux calculation
#define MAXNEIGHOUT  110
#define MAXNEIGHIN  10
#define MAXLINE 1024
#define DELTA 4
#define EXPAND 10
using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairCACSW::PairCACSW(LAMMPS *lmp) : PairCAC(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  type_map = 1;
  flux_max = add_ncluster = 0;
  outer_neighflag = 1;
  flux_enable = 2;
  nelements = 0;
  elements = NULL;
  nparams = maxparam = 0;
  params = NULL;
  elem2param = NULL;
  map = NULL;
  cluster_neighbors = NULL;
  cluster_neighbor_counts = NULL;
  add_cluster_neighbors = NULL;
  add_cluster_neighbor_counts = NULL;
}

/* ---------------------------------------------------------------------- */

PairCACSW::~PairCACSW() {

  if (copymode) return;

  if (elements)
    for (int i = 0; i < nelements; i++) delete[] elements[i];
    delete[] elements;
    memory->destroy(params);
    memory->destroy(elem2param);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    delete[] map;
  }

  memory->destroy(inner_neighbor_coords);
  memory->destroy(outer_neighbor_coords);
  memory->destroy(inner_neighbor_types);
  memory->destroy(outer_neighbor_types);
  memory->destroy(cluster_neighbor_counts);
  if(local_inner_max)
  for(int i=0; i < local_inner_max +1; i++)
    memory->destroy(cluster_neighbors[i]);
  memory->destroy(cluster_neighbors);

  if(atom->cac_flux_flag){
  memory->destroy(add_cluster_neighbor_counts);
  for(int i=0; i < flux_max; i++)
    memory->destroy(add_cluster_neighbors[i]);
  memory->destroy(add_cluster_neighbors);
  }
}

/* ---------------------------------------------------------------------- */

void PairCACSW::allocate()
{
  allocated = 1;
  int n = atom->ntypes;
  max_nodes_per_element = atom->nodes_per_element;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");

  map = new int[n + 1];

  memory->create(mass_matrix,max_nodes_per_element, max_nodes_per_element,"pairCAC:mass_matrix");
  memory->create(mass_copy, max_nodes_per_element, max_nodes_per_element,"pairCAC:copy_mass_matrix");
  memory->create(force_column, max_nodes_per_element,3,"pairCAC:force_residue");
  memory->create(current_force_column, max_nodes_per_element,"pairCAC:current_force_residue");
  memory->create(current_nodal_forces, max_nodes_per_element,"pairCAC:current_nodal_force");
  memory->create(pivot, max_nodes_per_element+1,"pairCAC:pivots");
  quadrature_init(2);
}

/* ----------------------------------------------------------------------
 set coeffs for one or more type pairs
 ------------------------------------------------------------------------- */

void PairCACSW::coeff(int narg, char **arg) {
  int i, j, n;

  if (!allocated) allocate();

  if (narg != 3 + atom->ntypes)
    error->all(FLERR, "Incorrect args for pair coefficients");

  // insure I,J args are * *

  if (strcmp(arg[0], "*") != 0 || strcmp(arg[1], "*") != 0)
    error->all(FLERR, "Incorrect args for pair coefficients");

  // read args that map atom types to elements in potential file
  // map[i] = which element the Ith atom type is, -1 if NULL
  // nelements = # of unique elements
  // elements = list of element names

  if (elements) {
    for (i = 0; i < nelements; i++) delete[] elements[i];
    delete[] elements;
  }
  elements = new char*[atom->ntypes];
  for (i = 0; i < atom->ntypes; i++) elements[i] = NULL;

  nelements = 0;
  for (i = 3; i < narg; i++) {
    if (strcmp(arg[i], "NULL") == 0) {
      map[i - 2] = -1;
      continue;
    }
    for (j = 0; j < nelements; j++)
      if (strcmp(arg[i], elements[j]) == 0) break;
    map[i - 2] = j;
    if (j == nelements) {
      n = strlen(arg[i]) + 1;
      elements[j] = new char[n];
      strcpy(elements[j], arg[i]);
      nelements++;
    }
  }

  // read potential file and initialize potential parameters

  read_file(arg[2]);
  setup_params();

  // clear setflag since coeff() called once with I,J = * *

  n = atom->ntypes;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements

  int count = 0;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      if (map[i] >= 0 && map[j] >= 0) {
        setflag[i][j] = 1;
        count++;
      }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
  type_maps = map;
}

/* ----------------------------------------------------------------------
 init for one type pair i,j and corresponding j,i
 ------------------------------------------------------------------------- */

double PairCACSW::init_one(int i, int j) {

  if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");

  if (outer_neighflag)
  return 2*cutmax;
  else
  return cutmax;
}

/* ---------------------------------------------------------------------- */


void PairCACSW::init_style()
{
  PairCAC::init_style();
  atom->max_neigh_inner_init = maxneigh_quad_inner = MAXNEIGHIN;
  atom->max_neigh_outer_init = maxneigh_quad_outer = MAXNEIGHOUT;
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style cac/sw requires atom IDs");
  atom->outer_neigh_flag=1;
}

/* ---------------------------------------------------------------------- */

void PairCACSW::read_file(char *file)
{
  memory->sfree(params);
  params = nullptr;
  nparams = maxparam = 0;

  // open file on proc 0

  if (comm->me == 0) {
    PotentialFileReader reader(lmp, file, "sw", unit_convert_flag);
    char * line;

    // transparently convert units for supported conversions

    int unit_convert = reader.get_unit_convert();
    double conversion_factor = utils::get_conversion_factor(utils::ENERGY,
                                                            unit_convert);

    while((line = reader.next_line(NPARAMS_PER_LINE))) {
      try {
        ValueTokenizer values(line);

        std::string iname = values.next_string();
        std::string jname = values.next_string();
        std::string kname = values.next_string();

        // ielement,jelement,kelement = 1st args
        // if all 3 args are in element list, then parse this line
        // else skip to next entry in file
        int ielement, jelement, kelement;

        for (ielement = 0; ielement < nelements; ielement++)
          if (iname == elements[ielement]) break;
        if (ielement == nelements) continue;
        for (jelement = 0; jelement < nelements; jelement++)
          if (jname == elements[jelement]) break;
        if (jelement == nelements) continue;
        for (kelement = 0; kelement < nelements; kelement++)
          if (kname == elements[kelement]) break;
        if (kelement == nelements) continue;

        // load up parameter settings and error check their values

        if (nparams == maxparam) {
          maxparam += DELTA;
          params = (Param *) memory->srealloc(params,maxparam*sizeof(Param),
                                              "pair:params");

          // make certain all addional allocated storage is initialized
          // to avoid false positives when checking with valgrind

          memset(params + nparams, 0, DELTA*sizeof(Param));
        }

        params[nparams].ielement = ielement;
        params[nparams].jelement = jelement;
        params[nparams].kelement = kelement;
        params[nparams].epsilon  = values.next_double();
        params[nparams].sigma    = values.next_double();
        params[nparams].littlea  = values.next_double();
        params[nparams].lambda   = values.next_double();
        params[nparams].gamma    = values.next_double();
        params[nparams].costheta = values.next_double();
        params[nparams].biga     = values.next_double();
        params[nparams].bigb     = values.next_double();
        params[nparams].powerp   = values.next_double();
        params[nparams].powerq   = values.next_double();
        params[nparams].tol      = values.next_double();
      } catch (TokenizerException & e) {
        error->one(FLERR, e.what());
      }

      if (unit_convert) {
        params[nparams].epsilon *= conversion_factor;
      }

      if (params[nparams].epsilon < 0.0 || params[nparams].sigma < 0.0 ||
          params[nparams].littlea < 0.0 || params[nparams].lambda < 0.0 ||
          params[nparams].gamma < 0.0 || params[nparams].biga < 0.0 ||
          params[nparams].bigb < 0.0 || params[nparams].powerp < 0.0 ||
          params[nparams].powerq < 0.0 || params[nparams].tol < 0.0)
        error->one(FLERR,"Illegal Stillinger-Weber parameter");

      nparams++;
    }
  }

  MPI_Bcast(&nparams, 1, MPI_INT, 0, world);
  MPI_Bcast(&maxparam, 1, MPI_INT, 0, world);

  if(comm->me != 0) {
    params = (Param *) memory->srealloc(params,maxparam*sizeof(Param), "pair:params");
  }

  MPI_Bcast(params, maxparam*sizeof(Param), MPI_BYTE, 0, world);
}

/* ---------------------------------------------------------------------- */

void PairCACSW::setup_params()
{
  int i, j, k, m, n;
  double rtmp;

  // set elem2param for all triplet combinations
  // must be a single exact match to lines read from file
  // do not allow for ACB in place of ABC

  memory->destroy(elem2param);
  memory->create(elem2param, nelements, nelements, nelements, "pair:elem2param");

  for (i = 0; i < nelements; i++)
    for (j = 0; j < nelements; j++)
      for (k = 0; k < nelements; k++) {
        n = -1;
        for (m = 0; m < nparams; m++) {
          if (i == params[m].ielement && j == params[m].jelement &&
            k == params[m].kelement) {
            if (n >= 0) error->all(FLERR, "Potential file has duplicate entry");
            n = m;
          }
        }
        if (n < 0) error->all(FLERR, "Potential file is missing an entry");
        elem2param[i][j][k] = n;
      }


  // compute parameter values derived from inputs

  // set cutsq using shortcut to reduce neighbor list for accelerated
  // calculations. cut must remain unchanged as it is a potential parameter
  // (cut = a*sigma)

  for (m = 0; m < nparams; m++) {
    params[m].cut = params[m].sigma*params[m].littlea;

    rtmp = params[m].cut;
    if (params[m].tol > 0.0) {
      if (params[m].tol > 0.01) params[m].tol = 0.01;
      if (params[m].gamma < 1.0)
        rtmp = rtmp +
        params[m].gamma * params[m].sigma / log(params[m].tol);
      else rtmp = rtmp +
        params[m].sigma / log(params[m].tol);
    }
    params[m].cutsq = rtmp * rtmp;

    params[m].sigma_gamma = params[m].sigma*params[m].gamma;
    params[m].lambda_epsilon = params[m].lambda*params[m].epsilon;
    params[m].lambda_epsilon2 = 2.0*params[m].lambda*params[m].epsilon;
    params[m].c1 = params[m].biga*params[m].epsilon *
      params[m].powerp*params[m].bigb *
      pow(params[m].sigma, params[m].powerp);
    params[m].c2 = params[m].biga*params[m].epsilon*params[m].powerq *
      pow(params[m].sigma, params[m].powerq);
    params[m].c3 = params[m].biga*params[m].epsilon*params[m].bigb *
      pow(params[m].sigma, params[m].powerp + 1.0);
    params[m].c4 = params[m].biga*params[m].epsilon *
      pow(params[m].sigma, params[m].powerq + 1.0);
    params[m].c5 = params[m].biga*params[m].epsilon*params[m].bigb *
      pow(params[m].sigma, params[m].powerp);
    params[m].c6 = params[m].biga*params[m].epsilon *
      pow(params[m].sigma, params[m].powerq);
  }

  // set cutmax to max of all params

  cutmax = 0.0;
  for (m = 0; m < nparams; m++) {
    rtmp = sqrt(params[m].cutsq);
    if (rtmp > cutmax) cutmax = rtmp;
  }
  cut_global_s=cutmax;
}

/* ---------------------------------------------------------------------- */

void PairCACSW::twobody(Param *param, double rsq, double &fforce,
  int eflag, double &eng)
{
  double r, rinvsq, rp, rq, rainv, rainvsq, expsrainv;

  r = sqrt(rsq);
  rinvsq = 1.0 / rsq;
  rp = pow(r, -param->powerp);
  rq = pow(r, -param->powerq);
  rainv = 1.0 / (r - param->cut);
  rainvsq = rainv*rainv*r;
  expsrainv = exp(param->sigma * rainv);
  fforce = (param->c1*rp - param->c2*rq +
    (param->c3*rp - param->c4*rq) * rainvsq) * expsrainv * rinvsq;
  if (eflag) eng = (param->c5*rp - param->c6*rq) * expsrainv;
}

/* ---------------------------------------------------------------------- */

void PairCACSW::threebody(Param *paramij, Param *paramik, Param *paramijk,
  double rsq1, double rsq2,
  double *delr1, double *delr2,
  double *fj, double *fk, int eflag, double &eng)
{
  double r1, rinvsq1, rainv1, gsrainv1, gsrainvsq1, expgsrainv1;
  double r2, rinvsq2, rainv2, gsrainv2, gsrainvsq2, expgsrainv2;
  double rinv12, cs, delcs, delcssq, facexp, facrad, frad1, frad2;
  double facang, facang12, csfacang, csfac1, csfac2;

  r1 = sqrt(rsq1);
  rinvsq1 = 1.0 / rsq1;
  rainv1 = 1.0 / (r1 - paramij->cut);
  gsrainv1 = paramij->sigma_gamma * rainv1;
  gsrainvsq1 = gsrainv1*rainv1 / r1;
  expgsrainv1 = exp(gsrainv1);

  r2 = sqrt(rsq2);
  rinvsq2 = 1.0 / rsq2;
  rainv2 = 1.0 / (r2 - paramik->cut);
  gsrainv2 = paramik->sigma_gamma * rainv2;
  gsrainvsq2 = gsrainv2*rainv2 / r2;
  expgsrainv2 = exp(gsrainv2);

  rinv12 = 1.0 / (r1*r2);
  cs = (delr1[0] * delr2[0] + delr1[1] * delr2[1] + delr1[2] * delr2[2]) * rinv12;
  delcs = cs - paramijk->costheta;
  delcssq = delcs*delcs;

  facexp = expgsrainv1*expgsrainv2;

  // facrad = sqrt(paramij->lambda_epsilon*paramik->lambda_epsilon) *
  //          facexp*delcssq;

  facrad = paramijk->lambda_epsilon * facexp*delcssq;
  frad1 = facrad*gsrainvsq1;
  frad2 = facrad*gsrainvsq2;
  facang = paramijk->lambda_epsilon2 * facexp*delcs;
  facang12 = rinv12*facang;
  csfacang = cs*facang;
  csfac1 = rinvsq1*csfacang;

  fj[0] = delr1[0] * (frad1 + csfac1) - delr2[0] * facang12;
  fj[1] = delr1[1] * (frad1 + csfac1) - delr2[1] * facang12;
  fj[2] = delr1[2] * (frad1 + csfac1) - delr2[2] * facang12;

  csfac2 = rinvsq2*csfacang;

  fk[0] = delr2[0] * (frad2 + csfac2) - delr1[0] * facang12;
  fk[1] = delr2[1] * (frad2 + csfac2) - delr1[1] * facang12;
  fk[2] = delr2[2] * (frad2 + csfac2) - delr1[2] * facang12;

  if (eflag) eng = facrad;
}

//-----------------------------------------------------------------------

void PairCACSW::force_densities( int iii, double s,double t, double w, double coefficients,
  double &force_densityx,double &force_densityy,double &force_densityz){

  double delx,dely,delz;
  double cutshortsq = cutmax*cutmax;
  double fpair, flux_interaction[3];
  int *type = atom->type;
  double distancesq;
  double scan_position[3], scan_position2[3];
  int current_type = poly_counter;

  int nodes_per_element;
  int *nodes_count_list = atom->nodes_per_element_list;
  int origin_type = type_array[poly_counter];

  int listtype;
  int scan_type, scan_type2, scan_type3;
  int element_index;
  int poly_index, short_scan, short_scan2, short_scan3, short_scan_index;
  int neigh_max_inner = inner_quad_lists_counts[pqi];
  int neigh_max_outer = outer_quad_lists_counts[pqi];
  int neigh_max_add, all_neigh;
  int itype, jtype, ktype, ijparam, ikparam, ijkparam;
  double energy_contribution;
  double rsq, rsq1, rsq2;
  double delr1[3], delr2[3], fj[3], fk[3];
  double ****nodal_positions = atom->nodal_positions;
  int **node_types = atom->node_types;
  origin_type = map[type_array[poly_counter]];
  double inner_scan_position[3];
  int **inner_quad_indices = inner_quad_lists_index[pqi];
  int **outer_quad_indices = outer_quad_lists_index[pqi];
  double cut_add = atom->cut_add;
  energy_contribution = 0;

  //allocate arrays that store neighbor information around just this quadrature point
  if(neigh_max_inner>local_inner_max){
  memory->grow(cluster_neighbor_counts, neigh_max_inner+EXPAND+1, "pair_cac_sw:cluster_neighbor_counts");
  cluster_neighbors = (int **) memory->srealloc(cluster_neighbors,sizeof(int *)*(neigh_max_inner+EXPAND+1), "pair_cac_sw:cluster_neighbors");
  //initialize sizing of cluster neighbors using neigh_max_inner
  for(int cinit = 0; cinit < neigh_max_inner + EXPAND + 1; cinit++){
    if(cinit>=local_inner_max+1||local_inner_max==0) cluster_neighbors[cinit] = NULL;
    memory->grow(cluster_neighbors[cinit], neigh_max_inner + neigh_max_outer+EXPAND, "pair_cac_sw:inner_neighbor_types");
  }
  }

  if(quad_flux_flag){
    neigh_max_add = add_quad_lists_counts[pqi];
    all_neigh = neigh_max_inner+neigh_max_outer+neigh_max_add;
    add_ncluster = 0;
    if(neigh_max_outer+neigh_max_add + 1 > flux_max){
      memory->grow(add_cluster_neighbor_counts, neigh_max_outer+neigh_max_add+EXPAND+1, "pair_cac_sw:cluster_neighbor_counts");
      add_cluster_neighbors = (int **) memory->srealloc(add_cluster_neighbors,
                              sizeof(int *)*(neigh_max_outer+neigh_max_add+EXPAND+1), "pair_cac_sw:add_cluster_neighbors");
      //initialize sizing of cluster neighbors using neigh_max_inner
      for(int cinit = 0; cinit < neigh_max_outer+neigh_max_add+EXPAND+1; cinit++){
        if(cinit>=flux_max) add_cluster_neighbors[cinit] = NULL;
          memory->grow(add_cluster_neighbors[cinit], all_neigh + EXPAND, "pair_cac_tersoff:inner_neighbor_types");
      }
      flux_max = neigh_max_outer+neigh_max_add + EXPAND + 1;
    }

    for(int cinit = 0; cinit < flux_max; cinit++)
      add_cluster_neighbor_counts[cinit] = 0;
  }
  allocate_quad_memory();
  //set virtual neighbor types, etc.
  init_quad_arrays();
  //compute virtual neighbor positions at the current timestep
  interpolation(iii, s, t, w);
  
  //initialize cluster counts
  for(int cinit = 0; cinit < local_inner_max + 1; cinit++)
    cluster_neighbor_counts[cinit] = 0;

  //two body contribution
  for (int l = 0; l < neigh_max_inner; l++) {
    scan_type = inner_neighbor_types[l];
    scan_position[0] = inner_neighbor_coords[l][0];
    scan_position[1] = inner_neighbor_coords[l][1];
    scan_position[2] = inner_neighbor_coords[l][2];

    delx = current_position[0] - scan_position[0];
    dely = current_position[1] - scan_position[1];
    delz = current_position[2] - scan_position[2];
    distancesq = delx*delx + dely*dely + delz*delz;

    //check to add particle to smaller neighbor list of this virtual atom
    if(distancesq < cutshortsq)
      cluster_neighbors[0][cluster_neighbor_counts[0]++] = l;
    
    //extend cluster for flux calculation
    if(quad_flux_flag){
      if(distancesq < (cutmax+cut_add)*(cutmax+cut_add))
      add_cluster_neighbors[0][add_cluster_neighbor_counts[0]++] = l;
    }

    ijparam = elem2param[origin_type][scan_type][scan_type];
    if (distancesq >= params[ijparam].cutsq) continue;

    twobody(&params[ijparam], distancesq, fpair, quad_eflag, energy_contribution);
        quadrature_energy += energy_contribution/2;

    force_densityx += delx*fpair;
    force_densityy += dely*fpair;
    force_densityz += delz*fpair;
  
    if(atom->CAC_virial){
    virial_density[0] += 0.5*delx*delx*fpair;
    virial_density[1] += 0.5*dely*dely*fpair;
    virial_density[2] += 0.5*delz*delz*fpair;
    virial_density[3] += 0.5*delx*dely*fpair;
    virial_density[4] += 0.5*delx*delz*fpair;
    virial_density[5] += 0.5*dely*delz*fpair;
    }

    //cac flux contribution due to current quadrature point and neighbor pair interactions
    if(quad_flux_flag)
      current_quad_flux(l,-delx*fpair,-dely*fpair,-delz*fpair);
  }

  //search lists for more add cluster sites for flux calculation
  if(quad_flux_flag){
    //possible virtual atoms in the outer neigh band
    for (int l = 0; l < neigh_max_outer; l++) {
      scan_type = outer_neighbor_types[l];
      scan_position[0] = outer_neighbor_coords[l][0];
      scan_position[1] = outer_neighbor_coords[l][1];
      scan_position[2] = outer_neighbor_coords[l][2];

      delx = current_position[0] - scan_position[0];
      dely = current_position[1] - scan_position[1];
      delz = current_position[2] - scan_position[2];
      distancesq = delx*delx + dely*dely + delz*delz;
    
      if(distancesq < (cutmax+cut_add)*(cutmax+cut_add))
      add_cluster_neighbors[0][add_cluster_neighbor_counts[0]++] = l + neigh_max_inner;
    }

    //possible virtual atoms in the additional band
    if(cut_add > cutmax)
    for (int l = 0; l < neigh_max_add; l++) {
      scan_type = add_neighbor_types[l];
      scan_position[0] = add_neighbor_coords[l][0];
      scan_position[1] = add_neighbor_coords[l][1];
      scan_position[2] = add_neighbor_coords[l][2];

      delx = current_position[0] - scan_position[0];
      dely = current_position[1] - scan_position[1];
      delz = current_position[2] - scan_position[2];
      distancesq = delx*delx + dely*dely + delz*delz;
    
      if(distancesq < (cutmax+cut_add)*(cutmax+cut_add))
      add_cluster_neighbors[0][add_cluster_neighbor_counts[0]++] = l + neigh_max_inner + neigh_max_outer;
    }
  }
  
  //construct the list of cluster neighbors using the short cutoff
  for (int l = 0; l < cluster_neighbor_counts[0]; l++) {
    short_scan = cluster_neighbors[0][l];
    scan_position[0] = inner_neighbor_coords[short_scan][0];
    scan_position[1] = inner_neighbor_coords[short_scan][1];
    scan_position[2] = inner_neighbor_coords[short_scan][2];
    
    for (int j = 0; j < neigh_max_inner; j++) {
      if(short_scan==j) continue;
      scan_position2[0] = inner_neighbor_coords[j][0];
      scan_position2[1] = inner_neighbor_coords[j][1];
      scan_position2[2] = inner_neighbor_coords[j][2];
      delr1[0] = scan_position2[0] - scan_position[0];
      delr1[1] = scan_position2[1] - scan_position[1];
      delr1[2] = scan_position2[2] - scan_position[2];
      distancesq = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
      if(distancesq < cutshortsq)
        cluster_neighbors[l+1][cluster_neighbor_counts[l+1]++] = j;
    }
    for (int j = 0; j < neigh_max_outer; j++) {
      scan_position2[0] = outer_neighbor_coords[j][0];
      scan_position2[1] = outer_neighbor_coords[j][1];
      scan_position2[2] = outer_neighbor_coords[j][2];
      delr1[0] = scan_position2[0] - scan_position[0];
      delr1[1] = scan_position2[1] - scan_position[1];
      delr1[2] = scan_position2[2] - scan_position[2];
      distancesq = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
      if(distancesq < cutshortsq)
        cluster_neighbors[l+1][cluster_neighbor_counts[l+1]++] = j + neigh_max_inner;
    }
  }

  //construct additional cluster neighbors for flux calculation
  if(quad_flux_flag)
  for (int l = 0; l < add_cluster_neighbor_counts[0]; l++) {
    short_scan_index = short_scan = add_cluster_neighbors[0][l];
    if(short_scan < neigh_max_inner){
      scan_position[0] = inner_neighbor_coords[short_scan][0];
      scan_position[1] = inner_neighbor_coords[short_scan][1];
      scan_position[2] = inner_neighbor_coords[short_scan][2];
    }
    else if(short_scan < neigh_max_inner + neigh_max_outer){
      short_scan_index -= neigh_max_inner;
      scan_position[0] = outer_neighbor_coords[short_scan_index][0];
      scan_position[1] = outer_neighbor_coords[short_scan_index][1];
      scan_position[2] = outer_neighbor_coords[short_scan_index][2];
    }
    else{
      short_scan_index -= neigh_max_inner + neigh_max_outer;
      scan_position[0] = add_neighbor_coords[short_scan_index][0];
      scan_position[1] = add_neighbor_coords[short_scan_index][1];
      scan_position[2] = add_neighbor_coords[short_scan_index][2];
    }
    
    for (int j = 0; j < neigh_max_inner; j++) {
      if(short_scan==j) continue;
      scan_position2[0] = inner_neighbor_coords[j][0];
      scan_position2[1] = inner_neighbor_coords[j][1];
      scan_position2[2] = inner_neighbor_coords[j][2];
      delr1[0] = scan_position2[0] - scan_position[0];
      delr1[1] = scan_position2[1] - scan_position[1];
      delr1[2] = scan_position2[2] - scan_position[2];
      distancesq = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
      if(distancesq < cutshortsq)
        add_cluster_neighbors[l+1][add_cluster_neighbor_counts[l+1]++] = j;
    }
    for (int j = 0; j < neigh_max_outer; j++) {
      if(short_scan-neigh_max_inner==j) continue;
      scan_position2[0] = outer_neighbor_coords[j][0];
      scan_position2[1] = outer_neighbor_coords[j][1];
      scan_position2[2] = outer_neighbor_coords[j][2];
      delr1[0] = scan_position2[0] - scan_position[0];
      delr1[1] = scan_position2[1] - scan_position[1];
      delr1[2] = scan_position2[2] - scan_position[2];
      distancesq = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
      if(distancesq < cutshortsq)
        add_cluster_neighbors[l+1][add_cluster_neighbor_counts[l+1]++] = j + neigh_max_inner;
    }
    for (int j = 0; j < neigh_max_add; j++) {
      if(short_scan-neigh_max_outer-neigh_max_inner==j) continue;
      scan_position2[0] = add_neighbor_coords[j][0];
      scan_position2[1] = add_neighbor_coords[j][1];
      scan_position2[2] = add_neighbor_coords[j][2];
      delr1[0] = scan_position2[0] - scan_position[0];
      delr1[1] = scan_position2[1] - scan_position[1];
      delr1[2] = scan_position2[2] - scan_position[2];
      distancesq = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];
      if(distancesq < cutshortsq)
        add_cluster_neighbors[l+1][add_cluster_neighbor_counts[l+1]++] = j + neigh_max_inner + neigh_max_outer;
    }
  }

  int loop_limit = cluster_neighbor_counts[0] - 1;
  if(quad_flux_flag) loop_limit = cluster_neighbor_counts[0];

  //ith three body contributions
  for (int l = 0; l < loop_limit; l++) {
    short_scan = cluster_neighbors[0][l];
    scan_type = inner_neighbor_types[short_scan];
    scan_position[0] = inner_neighbor_coords[short_scan][0];
    scan_position[1] = inner_neighbor_coords[short_scan][1];
    scan_position[2] = inner_neighbor_coords[short_scan][2];
    delr1[0] = scan_position[0] - current_position[0];
    delr1[1] = scan_position[1] - current_position[1];
    delr1[2] = scan_position[2] - current_position[2];

    ijparam = elem2param[origin_type][scan_type][scan_type];
    rsq1 = delr1[0] * delr1[0] + delr1[1] * delr1[1] + delr1[2] * delr1[2];
    if (rsq1 >= params[ijparam].cutsq) continue;
    if(quad_flux_flag) flux_interaction[0] = flux_interaction[1] = flux_interaction[2] = 0;

    for (int k = l + 1; k < cluster_neighbor_counts[0]; k++) {
      short_scan2 = cluster_neighbors[0][k];
      scan_type2 = inner_neighbor_types[short_scan2];
      scan_position2[0] = inner_neighbor_coords[short_scan2][0];
      scan_position2[1] = inner_neighbor_coords[short_scan2][1];
      scan_position2[2] = inner_neighbor_coords[short_scan2][2];

      ikparam = elem2param[origin_type][scan_type2][scan_type2];
      ijkparam = elem2param[origin_type][scan_type][scan_type2];

      delr2[0] = scan_position2[0] - current_position[0];
      delr2[1] = scan_position2[1] - current_position[1];
      delr2[2] = scan_position2[2] - current_position[2];
      rsq2 = delr2[0] * delr2[0] + delr2[1] * delr2[1] + delr2[2] * delr2[2];
      if (rsq2 >= params[ikparam].cutsq) continue;

      threebody(&params[ijparam], &params[ikparam], &params[ijkparam],
        rsq1, rsq2, delr1, delr2, fj, fk, quad_eflag, energy_contribution);
        quadrature_energy += energy_contribution/3;

      force_densityx -= fj[0] + fk[0];
      force_densityy -= fj[1] + fk[1];
      force_densityz -= fj[2] + fk[2];

      if(atom->CAC_virial){
        virial_density[0] += THIRD*(delr1[0]*fj[0] + delr2[0]*fk[0]);
        virial_density[1] += THIRD*(delr1[1]*fj[1] + delr2[1]*fk[1]);
        virial_density[2] += THIRD*(delr1[2]*fj[2] + delr2[2]*fk[2]);
        virial_density[3] += THIRD*(delr1[0]*fj[1] + delr2[0]*fk[1]);
        virial_density[4] += THIRD*(delr1[0]*fj[2] + delr2[0]*fk[2]);
        virial_density[5] += THIRD*(delr1[1]*fj[2] + delr2[1]*fk[2]);
      }

      if(quad_flux_flag){
        flux_interaction[0] += fj[0];
        flux_interaction[1] += fj[1];
        flux_interaction[2] += fj[2];
      }
    }

    if(quad_flux_flag)
      for(int k = 0; k < l; k++){
        short_scan2 = cluster_neighbors[0][k];
        scan_type2 = inner_neighbor_types[short_scan2];
        scan_position2[0] = inner_neighbor_coords[short_scan2][0];
        scan_position2[1] = inner_neighbor_coords[short_scan2][1];
        scan_position2[2] = inner_neighbor_coords[short_scan2][2];

        ikparam = elem2param[origin_type][scan_type2][scan_type2];
        ijkparam = elem2param[origin_type][scan_type][scan_type2];

        delr2[0] = scan_position2[0] - current_position[0];
        delr2[1] = scan_position2[1] - current_position[1];
        delr2[2] = scan_position2[2] - current_position[2];
        rsq2 = delr2[0] * delr2[0] + delr2[1] * delr2[1] + delr2[2] * delr2[2];
        if (rsq2 >= params[ikparam].cutsq) continue;
        ijkparam = elem2param[origin_type][scan_type2][scan_type];
        threebody(&params[ikparam], &params[ijparam], &params[ijkparam],
        rsq2, rsq1, delr2, delr1, fk, fj, quad_eflag, energy_contribution);
        flux_interaction[0] += fj[0];
        flux_interaction[1] += fj[1];
        flux_interaction[2] += fj[2];
      }
    //cac flux contribution due to current quadrature point and neighbor pair interactions
    if(quad_flux_flag)
      current_quad_flux(short_scan,flux_interaction[0],flux_interaction[1],flux_interaction[2]);
  }

  //jk three body contributions to i
  for (int l = 0; l < cluster_neighbor_counts[0]; l++) {
    short_scan = cluster_neighbors[0][l];
    scan_type = inner_neighbor_types[short_scan];
    scan_position[0] = inner_neighbor_coords[short_scan][0];
    scan_position[1] = inner_neighbor_coords[short_scan][1];
    scan_position[2] = inner_neighbor_coords[short_scan][2];

    delr1[0] = current_position[0] - scan_position[0];
    delr1[1] = current_position[1] - scan_position[1];
    delr1[2] = current_position[2] - scan_position[2];

    ijparam = elem2param[scan_type][origin_type][origin_type];

    rsq1 = delr1[0] * delr1[0] + delr1[1] * delr1[1] + delr1[2] * delr1[2];
    if (rsq1 >= params[ijparam].cutsq) continue;
    if(quad_flux_flag) flux_interaction[0] = flux_interaction[1] = flux_interaction[2] = 0;

    for (int k = 0; k < cluster_neighbor_counts[l+1]; k++) {
      short_scan2 = cluster_neighbors[l+1][k];
      if(short_scan2 >= neigh_max_inner){ 
        short_scan2 -= neigh_max_inner;
        scan_type2 = outer_neighbor_types[short_scan2];
        scan_position2[0] = outer_neighbor_coords[short_scan2][0];
        scan_position2[1] = outer_neighbor_coords[short_scan2][1];
        scan_position2[2] = outer_neighbor_coords[short_scan2][2];
      }
      else{
        scan_type2 = inner_neighbor_types[short_scan2];
        scan_position2[0] = inner_neighbor_coords[short_scan2][0];
        scan_position2[1] = inner_neighbor_coords[short_scan2][1];
        scan_position2[2] = inner_neighbor_coords[short_scan2][2];
      }

      ikparam = elem2param[scan_type][scan_type2][scan_type2];
      ijkparam = elem2param[scan_type][origin_type][scan_type2];

      delr2[0] = scan_position2[0] - scan_position[0];
      delr2[1] = scan_position2[1] - scan_position[1];
      delr2[2] = scan_position2[2] - scan_position[2];
      rsq2 = delr2[0] * delr2[0] + delr2[1] * delr2[1] + delr2[2] * delr2[2];
      if (rsq2 >= params[ikparam].cutsq) continue;

      threebody(&params[ijparam], &params[ikparam], &params[ijkparam],
        rsq1, rsq2, delr1, delr2, fj, fk, quad_eflag, energy_contribution);
        quadrature_energy += energy_contribution/3;

      force_densityx += fj[0];
      force_densityy += fj[1];
      force_densityz += fj[2];
      if(atom->CAC_virial){
        virial_density[0] += THIRD*(delr1[0]*fj[0] + delr2[0]*fk[0]);
        virial_density[1] += THIRD*(delr1[1]*fj[1] + delr2[1]*fk[1]);
        virial_density[2] += THIRD*(delr1[2]*fj[2] + delr2[2]*fk[2]);
        virial_density[3] += THIRD*(delr1[0]*fj[1] + delr2[0]*fk[1]);
        virial_density[4] += THIRD*(delr1[0]*fj[2] + delr2[0]*fk[2]);
        virial_density[5] += THIRD*(delr1[1]*fj[2] + delr2[1]*fk[2]);
      }
      if(quad_flux_flag){
        flux_interaction[0] += fj[0];
        flux_interaction[1] += fj[1];
        flux_interaction[2] += fj[2];
      }
    }
    //cac flux contribution due to current quadrature point and neighbor pair interactions
    if(quad_flux_flag)
      current_quad_flux(short_scan,-flux_interaction[0],-flux_interaction[1],-flux_interaction[2]);
  }
  //end of scanning loop
 
  //additional cac flux contributions due to neighbors interacting with neighbors
  //  in the vicinity of this quadrature point
  if (quad_flux_flag) {
    //compute_intersections();
    quad_neigh_flux();
  }
}

/* ---------------------------------------------------------------------- 
 Compute the cac flux density due to virtual neighbors around a quadrature point
---------------------------------------------------------------------- */

void PairCACSW::quad_neigh_flux(){
  int all_neigh, is, isl, normal_flag, sign, scan_type1, scan_type2, scan_type3, index, jindex;
  int short_scan_index, short_scan_index2, short_scan_index3, *current_cluster, current_ncluster;
  int dim1, dim2, intersection_flag, intersection_count, icontrib, short_scan, short_scan2, short_scan3;
  int ijparam, ikparam, ijkparam;
  double m, fpair, interaction_forceij[3], interaction_forceji[3], delxa[3], energy_contribution;
  double scan_position1[3], scan_position2[3], scan_position3[3], delx, dely, delz, distancesq;  
  double intersection_point[3], proj, lparam, planecoord, plane_limits[2][2];
  double rsq, rsq1, rsq2;
  double delr1[3], delr2[3], fj[3], fk[3];
  double vix, viy, viz, vjx, vjy, vjz, interactionx, interactiony, interactionz;
  double *box_center = atom->box_center;
  double *box_size = atom->box_size;
  int neigh_max_inner = inner_quad_lists_counts[pqi];
  int neigh_max_outer;
  int neigh_add = add_quad_lists_counts[pqi];
  neigh_max_outer = outer_quad_lists_counts[pqi];
  all_neigh = cluster_neighbor_counts[0] + add_cluster_neighbor_counts[0];
  //icontrib = 0;

  //determine which of the 6 planes of the atom box are intersected by a given i-j pair
  for(int ineigh=0; ineigh < all_neigh; ineigh++){
    if(ineigh<cluster_neighbor_counts[0]){
      short_scan_index = short_scan = cluster_neighbors[0][ineigh];
      current_ncluster = cluster_neighbor_counts[ineigh+1];
      current_cluster = cluster_neighbors[ineigh+1];
      scan_type1 = inner_neighbor_types[short_scan];
      scan_position1[0] = inner_neighbor_coords[short_scan][0];
      scan_position1[1] = inner_neighbor_coords[short_scan][1];
      scan_position1[2] = inner_neighbor_coords[short_scan][2];
      vix = inner_neighbor_velocities[short_scan][0];
      viy = inner_neighbor_velocities[short_scan][1];
      viz = inner_neighbor_velocities[short_scan][2];
    }
    else{
      short_scan_index = short_scan = add_cluster_neighbors[0][ineigh-cluster_neighbor_counts[0]];
      current_ncluster = add_cluster_neighbor_counts[ineigh-cluster_neighbor_counts[0]+1];
      current_cluster = add_cluster_neighbors[ineigh-cluster_neighbor_counts[0]+1];
      if(short_scan >= neigh_max_inner+neigh_max_outer){
        short_scan_index -= neigh_max_inner+neigh_max_outer;
        scan_type1 = add_neighbor_types[short_scan_index];
        scan_position1[0] = add_neighbor_coords[short_scan_index][0];
        scan_position1[1] = add_neighbor_coords[short_scan_index][1];
        scan_position1[2] = add_neighbor_coords[short_scan_index][2];
        vix = add_neighbor_velocities[short_scan_index][0];
        viy = add_neighbor_velocities[short_scan_index][1];
        viz = add_neighbor_velocities[short_scan_index][2];
      }
      else if(short_scan >= neigh_max_inner){
        short_scan_index -= neigh_max_inner;
        scan_type1 = outer_neighbor_types[short_scan_index];
        scan_position1[0] = outer_neighbor_coords[short_scan_index][0];
        scan_position1[1] = outer_neighbor_coords[short_scan_index][1];
        scan_position1[2] = outer_neighbor_coords[short_scan_index][2];
        vix = outer_neighbor_velocities[short_scan_index][0];
        viy = outer_neighbor_velocities[short_scan_index][1];
        viz = outer_neighbor_velocities[short_scan_index][2];
      }
      else{
        scan_type1 = inner_neighbor_types[short_scan_index];
        scan_position1[0] = inner_neighbor_coords[short_scan_index][0];
        scan_position1[1] = inner_neighbor_coords[short_scan_index][1];
        scan_position1[2] = inner_neighbor_coords[short_scan_index][2];
        vix = inner_neighbor_velocities[short_scan_index][0];
        viy = inner_neighbor_velocities[short_scan_index][1];
        viz = inner_neighbor_velocities[short_scan_index][2];
      }
    }
    for(int jneigh=0; jneigh < current_ncluster; jneigh++){
      intersection_count = 0;
      short_scan_index2 = short_scan2 = current_cluster[jneigh];
      if(short_scan2 >= neigh_max_inner+neigh_max_outer){
        short_scan_index2 -= neigh_max_inner+neigh_max_outer;
        scan_type2 = add_neighbor_types[short_scan_index2];
        scan_position2[0] = add_neighbor_coords[short_scan_index2][0];
        scan_position2[1] = add_neighbor_coords[short_scan_index2][1];
        scan_position2[2] = add_neighbor_coords[short_scan_index2][2];
        vjx = add_neighbor_velocities[short_scan_index2][0];
        vjy = add_neighbor_velocities[short_scan_index2][1];
        vjz = add_neighbor_velocities[short_scan_index2][2];
      }
      else if(short_scan2 >= neigh_max_inner){
        short_scan_index2 -= neigh_max_inner;
        scan_type2 = outer_neighbor_types[short_scan_index2];
        scan_position2[0] = outer_neighbor_coords[short_scan_index2][0];
        scan_position2[1] = outer_neighbor_coords[short_scan_index2][1];
        scan_position2[2] = outer_neighbor_coords[short_scan_index2][2];
        vjx = outer_neighbor_velocities[short_scan_index2][0];
        vjy = outer_neighbor_velocities[short_scan_index2][1];
        vjz = outer_neighbor_velocities[short_scan_index2][2];
      }
      else{
        scan_type2 = inner_neighbor_types[short_scan_index2];
        scan_position2[0] = inner_neighbor_coords[short_scan_index2][0];
        scan_position2[1] = inner_neighbor_coords[short_scan_index2][1];
        scan_position2[2] = inner_neighbor_coords[short_scan_index2][2];
        vjx = inner_neighbor_velocities[short_scan_index2][0];
        vjy = inner_neighbor_velocities[short_scan_index2][1];
        vjz = inner_neighbor_velocities[short_scan_index2][2];
      }
      
      delxa[0] = delx = scan_position1[0] - scan_position2[0];
      delxa[1] = dely = scan_position1[1] - scan_position2[1];
      delxa[2] = delz = scan_position1[2] - scan_position2[2];
      delr1[0] = - delx;
      delr1[1] = - dely;
      delr1[2] = - delz;
      rsq1 = delx*delx + dely*dely + delz*delz;

      //compute pair interaction
      //two body part

      ijparam = elem2param[scan_type1][scan_type2][scan_type2];
      if (rsq1 >= params[ijparam].cutsq) continue;
      interaction_forceij[0] = interaction_forceij[1] = interaction_forceij[2] = 0;
      
      if(short_scan_index2>short_scan_index){
      twobody(&params[ijparam], rsq1, fpair, 0, energy_contribution);
      interaction_forceij[0] -= delx*fpair;
      interaction_forceij[1] -= dely*fpair;
      interaction_forceij[2] -= delz*fpair;
      }

      //three body term
      //ith three body contributions
      for (int kneigh = 0; kneigh < current_ncluster; kneigh++) {
        if (jneigh==kneigh) continue;
        short_scan_index3 = short_scan3 = current_cluster[kneigh];
        if(short_scan3 >= neigh_max_inner+neigh_max_outer){
          short_scan_index3 -= neigh_max_inner+neigh_max_outer;
          scan_type3 = add_neighbor_types[short_scan_index3];
          scan_position3[0] = add_neighbor_coords[short_scan_index3][0];
          scan_position3[1] = add_neighbor_coords[short_scan_index3][1];
          scan_position3[2] = add_neighbor_coords[short_scan_index3][2];
        }
        else if(short_scan3 >= neigh_max_inner){
          short_scan_index3 -= neigh_max_inner;
          scan_type3 = outer_neighbor_types[short_scan_index3];
          scan_position3[0] = outer_neighbor_coords[short_scan_index3][0];
          scan_position3[1] = outer_neighbor_coords[short_scan_index3][1];
          scan_position3[2] = outer_neighbor_coords[short_scan_index3][2];
        }
        else{
          scan_type3 = inner_neighbor_types[short_scan_index3];
          scan_position3[0] = inner_neighbor_coords[short_scan_index3][0];
          scan_position3[1] = inner_neighbor_coords[short_scan_index3][1];
          scan_position3[2] = inner_neighbor_coords[short_scan_index3][2];
        }
        delr2[0] = scan_position1[0] - scan_position3[0];
        delr2[1] = scan_position1[1] - scan_position3[1];
        delr2[2] = scan_position1[2] - scan_position3[2];
        ikparam = elem2param[scan_type1][scan_type3][scan_type3];
        rsq2 = delr2[0] * delr2[0] + delr2[1] * delr2[1] + delr2[2] * delr2[2];
        if (rsq2 >= params[ikparam].cutsq) continue;
        if(kneigh > jneigh){
        ijkparam = elem2param[scan_type1][scan_type2][scan_type3];
        threebody(&params[ijparam], &params[ikparam], &params[ijkparam],
          rsq1, rsq2, delr1, delr2, fj, fk, 0, energy_contribution);

        interaction_forceij[0] += fj[0];
        interaction_forceij[1] += fj[1];
        interaction_forceij[2] += fj[2];
        }
        else{
        ijkparam = elem2param[scan_type1][scan_type3][scan_type2];
        threebody(&params[ikparam], &params[ijparam], &params[ijkparam],
          rsq2, rsq1, delr2, delr1, fk, fj, 0, energy_contribution);
        interaction_forceij[0] += fj[0];
        interaction_forceij[1] += fj[1];
        interaction_forceij[2] += fj[2];
        }
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
          flux_density[4*isl] += (interaction_forceij[0]*vjx + 
          interaction_forceij[1]*vjy+interaction_forceij[2]*vjz)*sign;
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
