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

#include "pair_cac_tersoff.h"
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "atom.h"
#include "force.h"
#include "comm.h"
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
#define MAXNEIGHOUT  110
#define MAXNEIGHIN  10
#define MAXLINE 1024
#define DELTA 4
#define EXPAND 10
using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairCACTersoff::PairCACTersoff(LAMMPS *lmp) : PairCAC(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;

  nelements = 0;
  elements = NULL;
  nparams = maxparam = 0;
  params = NULL;
  elem2param = NULL;
  map = NULL;
  outer_neighflag = 1;
  nmax = 0;
  maxshort = 10;
  neighshort = NULL;

  inner_neighbor_coords = NULL;
  outer_neighbor_coords = NULL;
  inner_neighbor_types = NULL;
  outer_neighbor_types = NULL;
  cluster_neighbors = NULL;
  cluster_neighbor_counts = NULL;
}

/* ---------------------------------------------------------------------- */

PairCACTersoff::~PairCACTersoff() {

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
    memory->destroy(neighshort);
    memory->destroy(inner_neighbor_coords);
    memory->destroy(outer_neighbor_coords);
    memory->destroy(inner_neighbor_types);
    memory->destroy(outer_neighbor_types);
  }
}

/* ---------------------------------------------------------------------- */

void PairCACTersoff::allocate()
{
  allocated = 1;
  int n = atom->ntypes;
  max_nodes_per_element = atom->nodes_per_element;

  memory->create(setflag, n + 1, n + 1, "pairCAC:setflag");
  memory->create(cutsq, n + 1, n + 1, "pairCAC:cutsq");
  memory->create(neighshort,maxshort,"pairCAC:neighshort");
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

void PairCACTersoff::coeff(int narg, char **arg) {

  int i,j,n;

  if (!allocated) allocate();

  if (narg != 3 + atom->ntypes)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // insure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // read args that map atom types to elements in potential file
  // map[i] = which element the Ith atom type is, -1 if NULL
  // nelements = # of unique elements
  // elements = list of element names

  if (elements) {
    for (i = 0; i < nelements; i++) delete [] elements[i];
    delete [] elements;
  }
  elements = new char*[atom->ntypes];
  for (i = 0; i < atom->ntypes; i++) elements[i] = NULL;

  nelements = 0;
  for (i = 3; i < narg; i++) {
    if (strcmp(arg[i],"NULL") == 0) {
      map[i-2] = -1;
      continue;
    }
    for (j = 0; j < nelements; j++)
      if (strcmp(arg[i],elements[j]) == 0) break;
    map[i-2] = j;
    if (j == nelements) {
      n = strlen(arg[i]) + 1;
      elements[j] = new char[n];
      strcpy(elements[j],arg[i]);
      nelements++;
    }
  }

  // read potential file and initialize potential parameters

  read_file(arg[2]);
  setup_params();

  // clear setflag since coeff() called once with I,J = * *

  n = atom->ntypes;
  for (i = 1; i <= n; i++)
    for (j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements

  int count = 0;
  for (i = 1; i <= n; i++)
    for (j = i; j <= n; j++)
      if (map[i] >= 0 && map[j] >= 0) {
        setflag[i][j] = 1;
        count++;
      }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
 init for one type pair i,j and corresponding j,i
 ------------------------------------------------------------------------- */

double PairCACTersoff::init_one(int i, int j) {

  if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");

  if (outer_neighflag)
  return 2*cutmax;
  else
  return cutmax;
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairCACTersoff::init_style()
{
  check_existence_flags();
  atom->max_neigh_inner_init = maxneigh_quad_inner = MAXNEIGHIN;
  atom->max_neigh_outer_init = maxneigh_quad_outer = MAXNEIGHOUT;

  atom->outer_neigh_flag=1;

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->cac = 1;
}

/* ---------------------------------------------------------------------- */

void PairCACTersoff::read_file(char *file)
{
  memory->sfree(params);
  params = nullptr;
  nparams = maxparam = 0;

  // open file on proc 0

  if (comm->me == 0) {
    PotentialFileReader reader(lmp, file, "tersoff", unit_convert_flag);
    char *line;

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

        params[nparams].ielement  = ielement;
        params[nparams].jelement  = jelement;
        params[nparams].kelement  = kelement;
        params[nparams].powerm    = values.next_double();
        params[nparams].gamma     = values.next_double();
        params[nparams].lam3      = values.next_double();
        params[nparams].c         = values.next_double();
        params[nparams].d         = values.next_double();
        params[nparams].h         = values.next_double();
        params[nparams].powern    = values.next_double();
        params[nparams].beta      = values.next_double();
        params[nparams].lam2      = values.next_double();
        params[nparams].bigb      = values.next_double();
        params[nparams].bigr      = values.next_double();
        params[nparams].bigd      = values.next_double();
        params[nparams].lam1      = values.next_double();
        params[nparams].biga      = values.next_double();
        params[nparams].powermint = int(params[nparams].powerm);

        if (unit_convert) {
          params[nparams].biga *= conversion_factor;
          params[nparams].bigb *= conversion_factor;
        }
      } catch (TokenizerException & e) {
        error->one(FLERR, e.what());
      }

      // currently only allow m exponent of 1 or 3
      if (params[nparams].c < 0.0 ||
          params[nparams].d < 0.0 ||
          params[nparams].powern < 0.0 ||
          params[nparams].beta < 0.0 ||
          params[nparams].lam2 < 0.0 ||
          params[nparams].bigb < 0.0 ||
          params[nparams].bigr < 0.0 ||
          params[nparams].bigd < 0.0 ||
          params[nparams].bigd > params[nparams].bigr ||
          params[nparams].lam1 < 0.0 ||
          params[nparams].biga < 0.0 ||
          params[nparams].powerm - params[nparams].powermint != 0.0 ||
          (params[nparams].powermint != 3 &&
          params[nparams].powermint != 1) ||
          params[nparams].gamma < 0.0)
        error->one(FLERR,"Illegal Tersoff parameter");

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

//-------------------------------------

void PairCACTersoff::setup_params()
{
  int i,j,k,m,n;

  // set elem2param for all element triplet combinations
  // must be a single exact match to lines read from file
  // do not allow for ACB in place of ABC

  memory->destroy(elem2param);
  memory->create(elem2param,nelements,nelements,nelements,"pair:elem2param");

  for (i = 0; i < nelements; i++)
    for (j = 0; j < nelements; j++)
      for (k = 0; k < nelements; k++) {
        n = -1;
        for (m = 0; m < nparams; m++) {
          if (i == params[m].ielement && j == params[m].jelement &&
              k == params[m].kelement) {
            if (n >= 0) error->all(FLERR,"Potential file has duplicate entry");
            n = m;
          }
        }
        if (n < 0) error->all(FLERR,"Potential file is missing an entry");
        elem2param[i][j][k] = n;
      }


  // compute parameter values derived from inputs

  for (m = 0; m < nparams; m++) {
    params[m].cut = params[m].bigr + params[m].bigd;
    params[m].cutsq = params[m].cut*params[m].cut;

    params[m].c1 = pow(2.0*params[m].powern*1.0e-16,-1.0/params[m].powern);
    params[m].c2 = pow(2.0*params[m].powern*1.0e-8,-1.0/params[m].powern);
    params[m].c3 = 1.0/params[m].c2;
    params[m].c4 = 1.0/params[m].c1;
  }

  // set cutmax to max of all params

  cutmax = 0.0;
  for (m = 0; m < nparams; m++)
    if (params[m].cut > cutmax) cutmax = params[m].cut;
  cut_global_s = cutmax;
}

/* ---------------------------------------------------------------------- */

void PairCACTersoff::repulsive(Param *param, double rsq, double &fforce,
                            int eflag, double &eng)
{
  double r,tmp_fc,tmp_fc_d,tmp_exp;

  r = sqrt(rsq);
  tmp_fc = ters_fc(r,param);
  tmp_fc_d = ters_fc_d(r,param);
  tmp_exp = exp(-param->lam1 * r);
  fforce = -param->biga * tmp_exp * (tmp_fc_d - tmp_fc*param->lam1) / r;
  if (eflag) eng = tmp_fc * param->biga * tmp_exp;
}

/* ---------------------------------------------------------------------- */

double PairCACTersoff::zeta(Param *param, double rsqij, double rsqik,
                         double *delrij, double *delrik)
{
  double rij,rik,costheta,arg,ex_delr;

  rij = sqrt(rsqij);
  rik = sqrt(rsqik);
  costheta = (delrij[0]*delrik[0] + delrij[1]*delrik[1] +
              delrij[2]*delrik[2]) / (rij*rik);

  if (param->powermint == 3) arg = pow(param->lam3 * (rij-rik),3.0);
  else arg = param->lam3 * (rij-rik);

  if (arg > 69.0776) ex_delr = 1.e30;
  else if (arg < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(arg);

  return ters_fc(rik,param) * ters_gijk(costheta,param) * ex_delr;
}

/* ---------------------------------------------------------------------- */

void PairCACTersoff::force_zeta(Param *param, double rsq, double zeta_ij,
                             double &fforce, double &prefactor,
                             int eflag, double &eng)
{
  double r,fa,fa_d,bij;

  r = sqrt(rsq);
  fa = ters_fa(r,param);
  fa_d = ters_fa_d(r,param);
  bij = ters_bij(zeta_ij,param);
  fforce = 0.5*bij*fa_d / r;
  prefactor = -0.5*fa * ters_bij_d(zeta_ij,param);
  if (eflag) eng = 0.5*bij*fa;
}

/* ----------------------------------------------------------------------
   attractive term
   use param_ij cutoff for rij test
   use param_ijk cutoff for rik test
------------------------------------------------------------------------- */

void PairCACTersoff::attractive(Param *param, double prefactor,
                             double rsqij, double rsqik,
                             double *delrij, double *delrik,
                             double *fi, double *fj, double *fk)
{
  double rij_hat[3],rik_hat[3];
  double rij,rijinv,rik,rikinv;

  rij = sqrt(rsqij);
  rijinv = 1.0/rij;
  vec3_scale(rijinv,delrij,rij_hat);

  rik = sqrt(rsqik);
  rikinv = 1.0/rik;
  vec3_scale(rikinv,delrik,rik_hat);

  ters_zetaterm_d(prefactor,rij_hat,rij,rik_hat,rik,fi,fj,fk,param);
}

/* ---------------------------------------------------------------------- */

double PairCACTersoff::ters_fc(double r, Param *param)
{
  double ters_R = param->bigr;
  double ters_D = param->bigd;

  if (r < ters_R-ters_D) return 1.0;
  if (r > ters_R+ters_D) return 0.0;
  return 0.5*(1.0 - sin(MY_PI2*(r - ters_R)/ters_D));
}

/* ---------------------------------------------------------------------- */

double PairCACTersoff::ters_fc_d(double r, Param *param)
{
  double ters_R = param->bigr;
  double ters_D = param->bigd;

  if (r < ters_R-ters_D) return 0.0;
  if (r > ters_R+ters_D) return 0.0;
  return -(MY_PI4/ters_D) * cos(MY_PI2*(r - ters_R)/ters_D);
}

/* ---------------------------------------------------------------------- */

double PairCACTersoff::ters_fa(double r, Param *param)
{
  if (r > param->bigr + param->bigd) return 0.0;
  return -param->bigb * exp(-param->lam2 * r) * ters_fc(r,param);
}

/* ---------------------------------------------------------------------- */

double PairCACTersoff::ters_fa_d(double r, Param *param)
{
  if (r > param->bigr + param->bigd) return 0.0;
  return param->bigb * exp(-param->lam2 * r) *
    (param->lam2 * ters_fc(r,param) - ters_fc_d(r,param));
}

/* ---------------------------------------------------------------------- */

double PairCACTersoff::ters_bij(double zeta, Param *param)
{
  double tmp = param->beta * zeta;
  if (tmp > param->c1) return 1.0/sqrt(tmp);
  if (tmp > param->c2)
    return (1.0 - pow(tmp,-param->powern) / (2.0*param->powern))/sqrt(tmp);
  if (tmp < param->c4) return 1.0;
  if (tmp < param->c3)
    return 1.0 - pow(tmp,param->powern)/(2.0*param->powern);
  return pow(1.0 + pow(tmp,param->powern), -1.0/(2.0*param->powern));
}

/* ---------------------------------------------------------------------- */

double PairCACTersoff::ters_bij_d(double zeta, Param *param)
{
  double tmp = param->beta * zeta;
  if (tmp > param->c1) return param->beta * -0.5*pow(tmp,-1.5);
  if (tmp > param->c2)
    return param->beta * (-0.5*pow(tmp,-1.5) *
                          // error in negligible 2nd term fixed 9/30/2015
                          // (1.0 - 0.5*(1.0 +  1.0/(2.0*param->powern)) *
                          (1.0 - (1.0 +  1.0/(2.0*param->powern)) *
                           pow(tmp,-param->powern)));
  if (tmp < param->c4) return 0.0;
  if (tmp < param->c3)
    return -0.5*param->beta * pow(tmp,param->powern-1.0);

  double tmp_n = pow(tmp,param->powern);
  return -0.5 * pow(1.0+tmp_n, -1.0-(1.0/(2.0*param->powern)))*tmp_n / zeta;
}

/* ---------------------------------------------------------------------- */

void PairCACTersoff::ters_zetaterm_d(double prefactor,
                                  double *rij_hat, double rij,
                                  double *rik_hat, double rik,
                                  double *dri, double *drj, double *drk,
                                  Param *param)
{
  double gijk,gijk_d,ex_delr,ex_delr_d,fc,dfc,cos_theta,tmp;
  double dcosdri[3],dcosdrj[3],dcosdrk[3];

  fc = ters_fc(rik,param);
  dfc = ters_fc_d(rik,param);
  if (param->powermint == 3) tmp = pow(param->lam3 * (rij-rik),3.0);
  else tmp = param->lam3 * (rij-rik);

  if (tmp > 69.0776) ex_delr = 1.e30;
  else if (tmp < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(tmp);

  if (param->powermint == 3)
    ex_delr_d = 3.0*pow(param->lam3,3.0) * pow(rij-rik,2.0)*ex_delr;
  else ex_delr_d = param->lam3 * ex_delr;

  cos_theta = vec3_dot(rij_hat,rik_hat);
  gijk = ters_gijk(cos_theta,param);
  gijk_d = ters_gijk_d(cos_theta,param);
  costheta_d(rij_hat,rij,rik_hat,rik,dcosdri,dcosdrj,dcosdrk);

  // compute the derivative wrt Ri
  // dri = -dfc*gijk*ex_delr*rik_hat;
  // dri += fc*gijk_d*ex_delr*dcosdri;
  // dri += fc*gijk*ex_delr_d*(rik_hat - rij_hat);

  vec3_scale(-dfc*gijk*ex_delr,rik_hat,dri);
  vec3_scaleadd(fc*gijk_d*ex_delr,dcosdri,dri,dri);
  vec3_scaleadd(fc*gijk*ex_delr_d,rik_hat,dri,dri);
  vec3_scaleadd(-fc*gijk*ex_delr_d,rij_hat,dri,dri);
  vec3_scale(prefactor,dri,dri);

  // compute the derivative wrt Rj
  // drj = fc*gijk_d*ex_delr*dcosdrj;
  // drj += fc*gijk*ex_delr_d*rij_hat;

  vec3_scale(fc*gijk_d*ex_delr,dcosdrj,drj);
  vec3_scaleadd(fc*gijk*ex_delr_d,rij_hat,drj,drj);
  vec3_scale(prefactor,drj,drj);

  // compute the derivative wrt Rk
  // drk = dfc*gijk*ex_delr*rik_hat;
  // drk += fc*gijk_d*ex_delr*dcosdrk;
  // drk += -fc*gijk*ex_delr_d*rik_hat;

  vec3_scale(dfc*gijk*ex_delr,rik_hat,drk);
  vec3_scaleadd(fc*gijk_d*ex_delr,dcosdrk,drk,drk);
  vec3_scaleadd(-fc*gijk*ex_delr_d,rik_hat,drk,drk);
  vec3_scale(prefactor,drk,drk);
}

/* ---------------------------------------------------------------------- */

void PairCACTersoff::costheta_d(double *rij_hat, double rij,
                             double *rik_hat, double rik,
                             double *dri, double *drj, double *drk)
{
  // first element is devative wrt Ri, second wrt Rj, third wrt Rk

  double cos_theta = vec3_dot(rij_hat,rik_hat);

  vec3_scaleadd(-cos_theta,rij_hat,rik_hat,drj);
  vec3_scale(1.0/rij,drj,drj);
  vec3_scaleadd(-cos_theta,rik_hat,rij_hat,drk);
  vec3_scale(1.0/rik,drk,drk);
  vec3_add(drj,drk,dri);
  vec3_scale(-1.0,dri,dri);
}

//-----------------------------------------------------------------------

void PairCACTersoff::force_densities( int iii, double s,double t, double w, double coefficients,
  double &force_densityx,double &force_densityy,double &force_densityz){

double delx,dely,delz;
double shape_func;
double shape_func2;
double scanning_unit_cell[3];
double fpair;
int *type = atom->type;
double unit_cell[3];
double distancesq;
double current_position[3];
double scan_position[3];
double rcut;
int current_type = poly_counter;

int nodes_per_element;
int *nodes_count_list = atom->nodes_per_element_list;
double cutshortsq = cutmax*cutmax;

unit_cell[0] = s;
unit_cell[1] = t;
unit_cell[2] = w;

//scan the surrounding unit cell locations in a cartesian grid
//of isoparametric space until the cutoff is exceeded
//for each grid scan

int distanceflag=0;
  current_position[0]=0;
  current_position[1]=0;
  current_position[2]=0;

  if (!atomic_flag) {
    nodes_per_element = nodes_count_list[current_element_type];
    for (int kkk = 0; kkk < nodes_per_element; kkk++) {
      shape_func = shape_function(unit_cell[0], unit_cell[1], unit_cell[2], 2, kkk + 1);
      current_position[0] += current_nodal_positions[kkk][0] * shape_func;
      current_position[1] += current_nodal_positions[kkk][1] * shape_func;
      current_position[2] += current_nodal_positions[kkk][2] * shape_func;
    }
  }
  else {
    current_position[0] = s;
    current_position[1] = t;
    current_position[2] = w;
  }

  rcut = cut_global_s;
  int origin_type = type_array[poly_counter];

  int listtype;
  int scan_type, scan_type2, scan_type3;
  int element_index;
  int poly_index;
  int neigh_max_inner = inner_quad_lists_counts[pqi];
  int neigh_max_outer = outer_quad_lists_counts[pqi];
  int itype, jtype, ktype, ijparam, ikparam, ijkparam;
  double energy_contribution;
  energy_contribution = 0;

  if(neigh_max_inner>local_inner_max){
  memory->grow(inner_neighbor_types, neigh_max_inner+EXPAND, "pair_cac_tersoff:inner_neighbor_types");
  memory->grow(inner_neighbor_coords, neigh_max_inner+EXPAND, 3, "pair_cac_tersoff:inner_neighbor_coords");
  memory->grow(neighshort, neigh_max_inner+EXPAND, "pair_cac_tersoff:neighshort");
  memory->grow(cluster_neighbor_counts, neigh_max_inner+EXPAND+1, "pair_cac_tersoff:cluster_neighbor_counts");
  cluster_neighbors = (int **) memory->srealloc(cluster_neighbors,sizeof(int *)*(neigh_max_inner+EXPAND+1), "pair_cac_tersoff:cluster_neighbors");
  //initialize sizing of cluster neighbors using neigh_max_inner
  for(int cinit = 0; cinit < neigh_max_inner + EXPAND + 1; cinit++){
    if(cinit>=local_inner_max+1||local_inner_max==0) cluster_neighbors[cinit] = NULL;
    memory->grow(cluster_neighbors[cinit], neigh_max_inner + neigh_max_outer+EXPAND, "pair_cac_tersoff:inner_neighbor_types");
  }
  local_inner_max=neigh_max_inner+EXPAND;
  }
  if(neigh_max_outer>local_outer_max){
  memory->grow(outer_neighbor_coords, neigh_max_outer+EXPAND, 3, "Pair_CAC_sw:outer_neighbor_coords");
  memory->grow(outer_neighbor_types, neigh_max_outer+EXPAND, "Pair_CAC_sw:outer_neighbor_types");
  local_outer_max=neigh_max_outer+EXPAND;
  }
  for(int cinit = 0; cinit < local_inner_max + 1; cinit++)
    cluster_neighbor_counts[cinit] = 0;

  tagint itag, jtag;
  double rsq, rsq1, rsq2, rsq3;
  double delr1[3], delr2[3], ndelr1[3], delr3[3], fi[3], fj[3], fk[3];
  double ****nodal_positions = atom->nodal_positions;
  int **node_types = atom->node_types;
  origin_type = map[type_array[poly_counter]];
  double scan_position2[3], scan_position3[3];
  int **inner_quad_indices = inner_quad_lists_index[pqi];
  int **outer_quad_indices = outer_quad_lists_index[pqi];

  //precompute virtual neighbor atom locations
  for (int l = 0; l < neigh_max_inner; l++) {
    //listtype = quad_list_container[iii].inner_list2ucell[neigh_quad_counter].cell_indexes[l][0];
    element_index = inner_quad_indices[l][0];
    poly_index = inner_quad_indices[l][1];
    inner_neighbor_types[l] = map[node_types[element_index][poly_index]];
  }
  for (int l = 0; l < neigh_max_outer; l++) {
    //listtype = quad_list_container[iii].inner_list2ucell[neigh_quad_counter].cell_indexes[l][0];
    element_index = outer_quad_indices[l][0];
    poly_index = outer_quad_indices[l][1];
    outer_neighbor_types[l] = map[node_types[element_index][poly_index]];
  }
  //compute virtual neighbor positions at the current timestep
  interpolation(iii);

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
    
    //check to add particle to neighshort list
    if(distancesq < cutshortsq)
      cluster_neighbors[0][cluster_neighbor_counts[0]++] = l;

    ijparam = elem2param[origin_type][scan_type][scan_type];
    if (distancesq >= params[ijparam].cutsq) continue;

    repulsive(&params[ijparam], distancesq, fpair, quad_eflag, energy_contribution);
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
  }

  int short_scan, short_scan2, short_scan3;
  double zeta_ij,prefactor,r_scan,fa_scan;

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
  
  //three body contributions
  for (int l = 0; l < cluster_neighbor_counts[0]; l++) {
    short_scan = cluster_neighbors[0][l];
    scan_type = inner_neighbor_types[short_scan];
    scan_position[0] = inner_neighbor_coords[short_scan][0];
    scan_position[1] = inner_neighbor_coords[short_scan][1];
    scan_position[2] = inner_neighbor_coords[short_scan][2];
    delr1[0] = delx = scan_position[0] - current_position[0];
    delr1[1] = dely = scan_position[1] - current_position[1];
    delr1[2] = delz = scan_position[2] - current_position[2];
    ndelr1[0] = -delr1[0];
    ndelr1[1] = -delr1[1];
    ndelr1[2] = -delr1[2];
    ijparam = elem2param[origin_type][scan_type][scan_type];
    rsq1 = delx * delx + dely * dely + delz * delz;
    if (rsq1 >= params[ijparam].cutsq) continue;
    //first tally pairwise forces due to zeta
    zeta_ij = 0;

    for (int j = 0; j < cluster_neighbor_counts[0]; j++) {
      short_scan2 = cluster_neighbors[0][j];
      if (short_scan==short_scan2) continue;
      scan_type2 = inner_neighbor_types[short_scan2];
      scan_position2[0] = inner_neighbor_coords[short_scan2][0];
      scan_position2[1] = inner_neighbor_coords[short_scan2][1];
      scan_position2[2] = inner_neighbor_coords[short_scan2][2];
      ijkparam = elem2param[origin_type][scan_type][scan_type2];

      delr2[0] = scan_position2[0] - current_position[0];
      delr2[1] = scan_position2[1] - current_position[1];
      delr2[2] = scan_position2[2] - current_position[2];
      rsq2 = delr2[0] * delr2[0] + delr2[1] * delr2[1] + delr2[2] * delr2[2];
      if (rsq2 >= params[ijkparam].cutsq) continue;
      zeta_ij += zeta(&params[ijkparam],rsq1,rsq2,delr1,delr2);
    }
    force_zeta(&params[ijparam],rsq1,zeta_ij,fpair,prefactor,quad_eflag, energy_contribution);
    quadrature_energy += energy_contribution/2;
    force_densityx += delr1[0]*fpair;
    force_densityy += delr1[1]*fpair;
    force_densityz += delr1[2]*fpair;
    
    if(atom->CAC_virial){
    virial_density[0] += 0.5*delx*delx*fpair;
    virial_density[1] += 0.5*dely*dely*fpair;
    virial_density[2] += 0.5*delz*delz*fpair;
    virial_density[3] += 0.5*delx*dely*fpair;
    virial_density[4] += 0.5*delx*delz*fpair;
    virial_density[5] += 0.5*dely*delz*fpair;
    }
    
    //tally three body term for ijl and ilj permutations of triplets that contribute to i
    for (int j = 0; j < cluster_neighbor_counts[0]; j++) {
      short_scan2 = cluster_neighbors[0][j];
      if (short_scan==short_scan2) continue;
      scan_type2 = inner_neighbor_types[short_scan2];
      scan_position2[0] = inner_neighbor_coords[short_scan2][0];
      scan_position2[1] = inner_neighbor_coords[short_scan2][1];
      scan_position2[2] = inner_neighbor_coords[short_scan2][2];
      ijkparam = elem2param[origin_type][scan_type][scan_type2];

      delr2[0] = scan_position2[0] - current_position[0];
      delr2[1] = scan_position2[1] - current_position[1];
      delr2[2] = scan_position2[2] - current_position[2];
      rsq2 = delr2[0] * delr2[0] + delr2[1] * delr2[1] + delr2[2] * delr2[2];
      if (rsq2 >= params[ijkparam].cutsq) continue;
      attractive(&params[ijkparam],prefactor,
                  rsq1,rsq2,delr1,delr2,fi,fj,fk);

      force_densityx += fi[0];
      force_densityy += fi[1];
      force_densityz += fi[2];
      if(atom->CAC_virial){
        virial_density[0] += THIRD*(delr1[0]*fj[0] + delr2[0]*fk[0]);
        virial_density[1] += THIRD*(delr1[1]*fj[1] + delr2[1]*fk[1]);
        virial_density[2] += THIRD*(delr1[2]*fj[2] + delr2[2]*fk[2]);
        virial_density[3] += THIRD*(delr1[0]*fj[1] + delr2[0]*fk[1]);
        virial_density[4] += THIRD*(delr1[0]*fj[2] + delr2[0]*fk[2]);
        virial_density[5] += THIRD*(delr1[1]*fj[2] + delr2[1]*fk[2]);
      }
    }
    
    zeta_ij = 0;
    ijparam = elem2param[scan_type][origin_type][origin_type];
    //add pair ji force contribution to i due to zeta from atoms in the short list
    for (int j = 0; j < cluster_neighbor_counts[l+1]; j++) {
      short_scan2 = cluster_neighbors[l+1][j];
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
      delr2[0] = scan_position2[0] - scan_position[0];
      delr2[1] = scan_position2[1] - scan_position[1];
      delr2[2] = scan_position2[2] - scan_position[2];
      ijkparam = elem2param[scan_type][origin_type][scan_type2];
      rsq2 = delr2[0] * delr2[0] + delr2[1] * delr2[1] + delr2[2] * delr2[2];
      if (rsq2 >= params[ijkparam].cutsq) continue;
      zeta_ij += zeta(&params[ijkparam],rsq1,rsq2,ndelr1,delr2);
    }
    
    force_zeta(&params[ijparam],rsq1,zeta_ij,fpair,prefactor,quad_eflag, energy_contribution);
    quadrature_energy += energy_contribution/2;
    force_densityx += delr1[0]*fpair;
    force_densityy += delr1[1]*fpair;
    force_densityz += delr1[2]*fpair;
    
    //tally three body term for jil and lij permutations of triplets that contribute to i
    for (int j = 0; j < cluster_neighbor_counts[l+1]; j++) {
      short_scan2 = cluster_neighbors[l+1][j];
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
      delr2[0] = scan_position2[0] - scan_position[0];
      delr2[1] = scan_position2[1] - scan_position[1];
      delr2[2] = scan_position2[2] - scan_position[2];
      ijkparam = elem2param[scan_type][origin_type][scan_type2];
      rsq2 = delr2[0] * delr2[0] + delr2[1] * delr2[1] + delr2[2] * delr2[2];
      if (rsq2 >= params[ijkparam].cutsq) continue;
      attractive(&params[ijkparam],prefactor,
                  rsq1,rsq2,ndelr1,delr2,fj,fi,fk);
      force_densityx += fi[0];
      force_densityy += fi[1];
      force_densityz += fi[2];
      if(atom->CAC_virial){
        virial_density[0] += THIRD*(delr1[0]*fj[0] + delr2[0]*fk[0]);
        virial_density[1] += THIRD*(delr1[1]*fj[1] + delr2[1]*fk[1]);
        virial_density[2] += THIRD*(delr1[2]*fj[2] + delr2[2]*fk[2]);
        virial_density[3] += THIRD*(delr1[0]*fj[1] + delr2[0]*fk[1]);
        virial_density[4] += THIRD*(delr1[0]*fj[2] + delr2[0]*fk[2]);
        virial_density[5] += THIRD*(delr1[1]*fj[2] + delr2[1]*fk[2]);
      }
    }

    //tally three body term for jli and lji permutations of triplets that contribute to i
    //compute zeta_jl and zeta_lj
    
    for (int j = 0; j < cluster_neighbor_counts[l+1]; j++) {
      zeta_ij = 0;
      short_scan2 = cluster_neighbors[l+1][j];
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
      delr2[0] = scan_position2[0] - scan_position[0];
      delr2[1] = scan_position2[1] - scan_position[1];
      delr2[2] = scan_position2[2] - scan_position[2];
      ijparam = elem2param[scan_type][scan_type2][scan_type2];
      rsq2 = delr2[0] * delr2[0] + delr2[1] * delr2[1] + delr2[2] * delr2[2];
      if (rsq2 >= params[ijparam].cutsq) continue;

      for (int k = 0; k < cluster_neighbor_counts[l+1]; k++) {
        short_scan3 = cluster_neighbors[l+1][k];
        if(short_scan3 >= neigh_max_inner){ 
          short_scan3 -= neigh_max_inner;
          if(short_scan2==short_scan3) continue;
          scan_type3 = outer_neighbor_types[short_scan3];
          scan_position3[0] = outer_neighbor_coords[short_scan3][0];
          scan_position3[1] = outer_neighbor_coords[short_scan3][1];
          scan_position3[2] = outer_neighbor_coords[short_scan3][2];
        }
        else{
          if(short_scan2==short_scan3) continue;
          scan_type3 = inner_neighbor_types[short_scan3];
          scan_position3[0] = inner_neighbor_coords[short_scan3][0];
          scan_position3[1] = inner_neighbor_coords[short_scan3][1];
          scan_position3[2] = inner_neighbor_coords[short_scan3][2];
        }
        delr3[0] = scan_position3[0] - scan_position[0];
        delr3[1] = scan_position3[1] - scan_position[1];
        delr3[2] = scan_position3[2] - scan_position[2];
        rsq3 = delr3[0] * delr3[0] + delr3[1] * delr3[1] + delr3[2] * delr3[2];
        ijkparam = elem2param[scan_type][scan_type2][scan_type3];
        if (rsq3 >= params[ijkparam].cutsq) continue;
        zeta_ij += zeta(&params[ijkparam],rsq2,rsq3,delr2,delr3);
      }
      //add contribution to zeta_ij from the virtual atom at the quadrature point
      delr3[0] = current_position[0] - scan_position[0];
      delr3[1] = current_position[1] - scan_position[1];
      delr3[2] = current_position[2] - scan_position[2];
      rsq3 = delr3[0] * delr3[0] + delr3[1] * delr3[1] + delr3[2] * delr3[2];
      ijkparam = elem2param[scan_type][scan_type2][origin_type];
      if (rsq3 >= params[ijkparam].cutsq) continue;
      zeta_ij += zeta(&params[ijkparam],rsq2,rsq3,delr2,delr3);
      //compute prefactor
      r_scan = sqrt(rsq2);
      fa_scan = ters_fa(r_scan,&params[ijparam]);
      prefactor = -0.5*fa_scan * ters_bij_d(zeta_ij,&params[ijparam]);

      //compute three body contribution to quadrature site
      attractive(&params[ijkparam],prefactor,
                  rsq2,rsq3,delr2,delr3,fj,fk,fi);
      force_densityx += fi[0];
      force_densityy += fi[1];
      force_densityz += fi[2];
      if(atom->CAC_virial){
        virial_density[0] += THIRD*(delr1[0]*fj[0] + delr2[0]*fk[0]);
        virial_density[1] += THIRD*(delr1[1]*fj[1] + delr2[1]*fk[1]);
        virial_density[2] += THIRD*(delr1[2]*fj[2] + delr2[2]*fk[2]);
        virial_density[3] += THIRD*(delr1[0]*fj[1] + delr2[0]*fk[1]);
        virial_density[4] += THIRD*(delr1[0]*fj[2] + delr2[0]*fk[2]);
        virial_density[5] += THIRD*(delr1[1]*fj[2] + delr2[1]*fk[2]);
      }
    }
    
  }

//end of scanning loop
}
//------------------------------------------------------------------------
