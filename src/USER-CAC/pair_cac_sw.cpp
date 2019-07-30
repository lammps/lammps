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
#include "pair_cac_sw.h"
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

//#include "math_extra.h"
#define MAXNEIGH1  110
#define MAXNEIGH2  10
#define MAXLINE 1024
#define DELTA 4
using namespace LAMMPS_NS;
using namespace MathConst;
using namespace std;

/* ---------------------------------------------------------------------- */

PairCACSW::PairCACSW(LAMMPS *lmp) : PairCAC(lmp)
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
  
  inner_neighbor_coords = NULL;
  outer_neighbor_coords = NULL;
  inner_neighbor_types = NULL;
  outer_neighbor_types = NULL;
  interior_scales = NULL;
  surface_counts = NULL;

  surface_counts_max[0] = 0;
  surface_counts_max[1] = 0;
  surface_counts_max[2] = 0;
  surface_counts_max_old[0] = 0;
  surface_counts_max_old[1] = 0;
  surface_counts_max_old[2] = 0;
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
	memory->destroy(inner_neighbor_coords);
	memory->destroy(outer_neighbor_coords);
	memory->destroy(inner_neighbor_types);
	memory->destroy(outer_neighbor_types);
    memory->destroy(mass_matrix);
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
  memory->create(surf_set, 6, 2, "pairCAC:surf_set");
  memory->create(dof_set, 6, 4, "pairCAC:surf_set");
  memory->create(sort_surf_set, 6, 2, "pairCAC:surf_set");
  memory->create(sort_dof_set, 6, 4, "pairCAC:surf_set");
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
  check_existence_flags();
	maxneigh_quad_inner = MAXNEIGH2;
	maxneigh_quad_outer = MAXNEIGH1;
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style cac/sw requires atom IDs");

  atom->outer_neigh_flag=1;
  // need a full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  //neighbor->requests[irequest]->full = 1;
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

  for (int si = 0; si < 6; si++) {
	  sort_dof_set[si][0] = dof_set[si][0];
	  sort_dof_set[si][1] = dof_set[si][1];
	  sort_dof_set[si][2] = dof_set[si][2];
	  sort_dof_set[si][3] = dof_set[si][3];
	  sort_surf_set[si][0] = surf_set[si][0];
	  sort_surf_set[si][1] = surf_set[si][1];
  }

}

////////////////////////
void PairCACSW::read_file(char *file)
{
	int params_per_line = 14;
	char **words = new char*[params_per_line + 1];

	memory->sfree(params);
	params = NULL;
	nparams = maxparam = 0;

	// open file on proc 0

	FILE *fp;
	if (comm->me == 0) {
		fp = force->open_potential(file);
		if (fp == NULL) {
			char str[128];
			sprintf(str, "Cannot open Stillinger-Weber potential file %s", file);
			error->one(FLERR, str);
		}
	}

	// read each set of params from potential file
	// one set of params can span multiple lines
	// store params if all 3 element tags are in element list

	int n, nwords, ielement, jelement, kelement;
	char line[MAXLINE], *ptr;
	int eof = 0;

	while (1) {
		if (comm->me == 0) {
			ptr = fgets(line, MAXLINE, fp);
			if (ptr == NULL) {
				eof = 1;
				fclose(fp);
			}
			else n = strlen(line) + 1;
		}
		MPI_Bcast(&eof, 1, MPI_INT, 0, world);
		if (eof) break;
		MPI_Bcast(&n, 1, MPI_INT, 0, world);
		MPI_Bcast(line, n, MPI_CHAR, 0, world);

		// strip comment, skip line if blank

		if ((ptr = strchr(line, '#'))) *ptr = '\0';
		nwords = atom->count_words(line);
		if (nwords == 0) continue;

		// concatenate additional lines until have params_per_line words

		while (nwords < params_per_line) {
			n = strlen(line);
			if (comm->me == 0) {
				ptr = fgets(&line[n], MAXLINE - n, fp);
				if (ptr == NULL) {
					eof = 1;
					fclose(fp);
				}
				else n = strlen(line) + 1;
			}
			MPI_Bcast(&eof, 1, MPI_INT, 0, world);
			if (eof) break;
			MPI_Bcast(&n, 1, MPI_INT, 0, world);
			MPI_Bcast(line, n, MPI_CHAR, 0, world);
			if ((ptr = strchr(line, '#'))) *ptr = '\0';
			nwords = atom->count_words(line);
		}

		if (nwords != params_per_line)
			error->all(FLERR, "Incorrect format in Stillinger-Weber potential file");

		// words = ptrs to all words in line

		nwords = 0;
		words[nwords++] = strtok(line, " \t\n\r\f");
		while ((words[nwords++] = strtok(NULL, " \t\n\r\f"))) continue;

		// ielement,jelement,kelement = 1st args
		// if all 3 args are in element list, then parse this line
		// else skip to next entry in file

		for (ielement = 0; ielement < nelements; ielement++)
			if (strcmp(words[0], elements[ielement]) == 0) break;
		if (ielement == nelements) continue;
		for (jelement = 0; jelement < nelements; jelement++)
			if (strcmp(words[1], elements[jelement]) == 0) break;
		if (jelement == nelements) continue;
		for (kelement = 0; kelement < nelements; kelement++)
			if (strcmp(words[2], elements[kelement]) == 0) break;
		if (kelement == nelements) continue;

		// load up parameter settings and error check their values

		if (nparams == maxparam) {
			maxparam += DELTA;
			params = (Param *)memory->srealloc(params, maxparam * sizeof(Param),
				"pair:params");
		}

		params[nparams].ielement = ielement;
		params[nparams].jelement = jelement;
		params[nparams].kelement = kelement;
		params[nparams].epsilon = atof(words[3]);
		params[nparams].sigma = atof(words[4]);
		params[nparams].littlea = atof(words[5]);
		params[nparams].lambda = atof(words[6]);
		params[nparams].gamma = atof(words[7]);
		params[nparams].costheta = atof(words[8]);
		params[nparams].biga = atof(words[9]);
		params[nparams].bigb = atof(words[10]);
		params[nparams].powerp = atof(words[11]);
		params[nparams].powerq = atof(words[12]);
		params[nparams].tol = atof(words[13]);

		if (params[nparams].epsilon < 0.0 || params[nparams].sigma < 0.0 ||
			params[nparams].littlea < 0.0 || params[nparams].lambda < 0.0 ||
			params[nparams].gamma < 0.0 || params[nparams].biga < 0.0 ||
			params[nparams].bigb < 0.0 || params[nparams].powerp < 0.0 ||
			params[nparams].powerq < 0.0 || params[nparams].tol < 0.0)
			error->all(FLERR, "Illegal Stillinger-Weber parameter");

		nparams++;
	}

	delete[] words;
}

//-------------------------------------

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

double r2inv;
double r6inv;
double shape_func;
double shape_func2;
double unit_cell_mapped[3];
double scanning_unit_cell[3];

double forcelj,factor_lj,fpair;
int *type = atom->type;
double unit_cell[3];
double distancesq;
double current_position[3];
double scan_position[3];
double rcut;
int current_type = poly_counter;

int nodes_per_element;
int *nodes_count_list = atom->nodes_per_element_list;	

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

int distanceflag=0;
    current_position[0]=0;
    current_position[1]=0;
    current_position[2]=0;

	if (!atomic_flag) {
		nodes_per_element = nodes_count_list[current_element_type];
		for (int kkk = 0; kkk < nodes_per_element; kkk++) {
			shape_func = shape_function(unit_cell[0], unit_cell[1], unit_cell[2], 2, kkk + 1);
			current_position[0] += current_nodal_positions[kkk][poly_counter][0] * shape_func;
			current_position[1] += current_nodal_positions[kkk][poly_counter][1] * shape_func;
			current_position[2] += current_nodal_positions[kkk][poly_counter][2] * shape_func;
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
	int scan_type, scan_type2;
	int listindex;
	int poly_index;
	int element_index;
	int *ilist, *jlist, *numneigh, **firstneigh;
	int neigh_max_inner = inner_quad_lists_counts[iii][neigh_quad_counter];
	int neigh_max_outer = outer_quad_lists_counts[iii][neigh_quad_counter];
	int itype, jtype, ktype, ijparam, ikparam, ijkparam;
	double energy_contribution;
	energy_contribution = 0;

	if(neigh_max_inner>local_inner_max){
	memory->grow(inner_neighbor_types, neigh_max_inner, "Pair_CAC_sw:inner_neighbor_types");
	memory->grow(inner_neighbor_coords, neigh_max_inner, 3, "Pair_CAC_sw:inner_neighbor_coords");
	local_inner_max=neigh_max_inner;
	}
	if(neigh_max_outer>local_outer_max){
	memory->grow(outer_neighbor_coords, neigh_max_outer, 3, "Pair_CAC_sw:outer_neighbor_coords");
	memory->grow(outer_neighbor_types, neigh_max_outer, "Pair_CAC_sw:outer_neighbor_types");
	local_outer_max=neigh_max_outer;
	}

	tagint itag, jtag;
	double rsq, rsq1, rsq2;
	double delr1[3], delr2[3], fj[3], fk[3];
	ilist = list->ilist;
	numneigh = list->numneigh;
	firstneigh = list->firstneigh;
	jlist = firstneigh[iii];
	double ****nodal_positions = atom->nodal_positions;
	int **node_types = atom->node_types;
	origin_type = map[type_array[poly_counter]];
	double inner_scan_position[3];
	//precompute virtual neighbor atom locations
	for (int l = 0; l < neigh_max_inner; l++) {
		scanning_unit_cell[0] = inner_quad_lists_ucell[iii][neigh_quad_counter][l][0];
		scanning_unit_cell[1] = inner_quad_lists_ucell[iii][neigh_quad_counter][l][1];
		scanning_unit_cell[2] = inner_quad_lists_ucell[iii][neigh_quad_counter][l][2];
		//listtype = quad_list_container[iii].inner_list2ucell[neigh_quad_counter].cell_indexes[l][0];
		listindex = inner_quad_lists_index[iii][neigh_quad_counter][l][0];
		poly_index = inner_quad_lists_index[iii][neigh_quad_counter][l][1];
		element_index = listindex;
		element_index &= NEIGHMASK;
		inner_neighbor_types[l] = map[node_types[element_index][poly_index]];
		neigh_list_cord(inner_neighbor_coords[l][0], inner_neighbor_coords[l][1], inner_neighbor_coords[l][2],
			element_index, poly_index, scanning_unit_cell[0], scanning_unit_cell[1], scanning_unit_cell[2]);
	}
	for (int l = 0; l < neigh_max_outer; l++) {
        scanning_unit_cell[0] = outer_quad_lists_ucell[iii][neigh_quad_counter][l][0];
		scanning_unit_cell[1] = outer_quad_lists_ucell[iii][neigh_quad_counter][l][1];
		scanning_unit_cell[2] = outer_quad_lists_ucell[iii][neigh_quad_counter][l][2];
		//listtype = quad_list_container[iii].inner_list2ucell[neigh_quad_counter].cell_indexes[l][0];
		listindex = outer_quad_lists_index[iii][neigh_quad_counter][l][0];
		poly_index = outer_quad_lists_index[iii][neigh_quad_counter][l][1];
		element_index = listindex;
		element_index &= NEIGHMASK;
		outer_neighbor_types[l] = map[node_types[element_index][poly_index]];
		neigh_list_cord(outer_neighbor_coords[l][0], outer_neighbor_coords[l][1], outer_neighbor_coords[l][2],
			element_index, poly_index, scanning_unit_cell[0], scanning_unit_cell[1], scanning_unit_cell[2]);
	}
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

		ijparam = elem2param[origin_type][scan_type][scan_type];
		if (distancesq >= params[ijparam].cutsq) continue;

		twobody(&params[ijparam], distancesq, fpair, quad_eflag, energy_contribution);
        quadrature_energy += energy_contribution/2;

		force_densityx += delx*fpair;
		force_densityy += dely*fpair;
		force_densityz += delz*fpair;
		if(atom->CAC_virial){
		virial_density[0] +=delx*delx*fpair;
		virial_density[1] +=dely*dely*fpair;
		virial_density[2] +=delz*delz*fpair;
		virial_density[3] +=delx*dely*fpair;
		virial_density[4] +=delx*delz*fpair;
		virial_density[5] +=dely*delz*fpair;
		}


	}

	//ith three body contributions
	for (int l = 0; l < neigh_max_inner - 1; l++) {
		scan_type = inner_neighbor_types[l];
		scan_position[0] = inner_neighbor_coords[l][0];
		scan_position[1] = inner_neighbor_coords[l][1];
		scan_position[2] = inner_neighbor_coords[l][2];
		delr1[0] = scan_position[0] - current_position[0];
		delr1[1] = scan_position[1] - current_position[1];
		delr1[2] = scan_position[2] - current_position[2];
		distancesq = delx*delx + dely*dely + delz*delz;

		ijparam = elem2param[origin_type][scan_type][scan_type];

		rsq1 = delr1[0] * delr1[0] + delr1[1] * delr1[1] + delr1[2] * delr1[2];
		if (rsq1 >= params[ijparam].cutsq) continue;



		for (int k = l + 1; k < neigh_max_inner; k++) {
			scan_type2 = inner_neighbor_types[k];
			scan_position[0] = inner_neighbor_coords[k][0];
			scan_position[1] = inner_neighbor_coords[k][1];
			scan_position[2] = inner_neighbor_coords[k][2];

			ikparam = elem2param[origin_type][scan_type2][scan_type2];
			ijkparam = elem2param[origin_type][scan_type][scan_type2];

			delr2[0] = scan_position[0] - current_position[0];
			delr2[1] = scan_position[1] - current_position[1];
			delr2[2] = scan_position[2] - current_position[2];
			rsq2 = delr2[0] * delr2[0] + delr2[1] * delr2[1] + delr2[2] * delr2[2];
			if (rsq2 >= params[ikparam].cutsq) continue;


			threebody(&params[ijparam], &params[ikparam], &params[ijkparam],
				rsq1, rsq2, delr1, delr2, fj, fk, quad_eflag, energy_contribution);
            quadrature_energy += energy_contribution/3;

			force_densityx -= fj[0] + fk[0];
			force_densityy -= fj[1] + fk[1];
			force_densityz -= fj[2] + fk[2];
      if(atom->CAC_virial){
	   	virial_density[0] += delr1[0]*fj[0] + delr2[0]*fk[0];
		  virial_density[1] += delr1[1]*fj[1] + delr2[1]*fk[1];
		  virial_density[2] += delr1[2]*fj[2] + delr2[2]*fk[2];
		  virial_density[3] += delr1[0]*fj[1] + delr2[0]*fk[1];
		  virial_density[4] += delr1[0]*fj[2] + delr2[0]*fk[2];
		  virial_density[5] += delr1[1]*fj[2] + delr2[1]*fk[2];
		  }

		}

	}
	//jk three body contributions to i
	for (int l = 0; l < neigh_max_inner; l++) {
		scan_type = inner_neighbor_types[l];
		inner_scan_position[0] = inner_neighbor_coords[l][0];
		inner_scan_position[1] = inner_neighbor_coords[l][1];
		inner_scan_position[2] = inner_neighbor_coords[l][2];

		delr1[0] = current_position[0] - inner_scan_position[0];
		delr1[1] = current_position[1] - inner_scan_position[1];
		delr1[2] = current_position[2] - inner_scan_position[2];
		distancesq = delx*delx + dely*dely + delz*delz;

		ijparam = elem2param[origin_type][scan_type][scan_type];

		rsq1 = delr1[0] * delr1[0] + delr1[1] * delr1[1] + delr1[2] * delr1[2];
		if (rsq1 >= params[ijparam].cutsq) continue;



		for (int k = 0; k < neigh_max_inner; k++) {
			//add ji as well as ij contributions
			if (k == l) continue;
			scan_type2 = inner_neighbor_types[k];
			scan_position[0] = inner_neighbor_coords[k][0];
			scan_position[1] = inner_neighbor_coords[k][1];
			scan_position[2] = inner_neighbor_coords[k][2];

			ikparam = elem2param[origin_type][scan_type2][scan_type2];
			ijkparam = elem2param[origin_type][scan_type][scan_type2];

			delr2[0] = scan_position[0] - inner_scan_position[0];
			delr2[1] = scan_position[1] - inner_scan_position[1];
			delr2[2] = scan_position[2] - inner_scan_position[2];
			rsq2 = delr2[0] * delr2[0] + delr2[1] * delr2[1] + delr2[2] * delr2[2];
			if (rsq2 >= params[ikparam].cutsq) continue;


			threebody(&params[ijparam], &params[ikparam], &params[ijkparam],
				rsq1, rsq2, delr1, delr2, fj, fk, quad_eflag, energy_contribution);
            quadrature_energy += energy_contribution/3;

			force_densityx += fj[0];
			force_densityy += fj[1];
			force_densityz += fj[2];
			if(atom->CAC_virial){
	   	virial_density[0] += delr1[0]*fj[0] + delr2[0]*fk[0];
		  virial_density[1] += delr1[1]*fj[1] + delr2[1]*fk[1];
		  virial_density[2] += delr1[2]*fj[2] + delr2[2]*fk[2];
		  virial_density[3] += delr1[0]*fj[1] + delr2[0]*fk[1];
		  virial_density[4] += delr1[0]*fj[2] + delr2[0]*fk[2];
		  virial_density[5] += delr1[1]*fj[2] + delr2[1]*fk[2];
		  }
		}
		for (int k = 0; k < neigh_max_outer; k++) {
			//add contributions that come from outer neighbor band (farther particle triplets connected to i)

			scan_type2 = outer_neighbor_types[k];
			scan_position[0] = outer_neighbor_coords[k][0];
			scan_position[1] = outer_neighbor_coords[k][1];
			scan_position[2] = outer_neighbor_coords[k][2];

			ikparam = elem2param[origin_type][scan_type2][scan_type2];
			ijkparam = elem2param[origin_type][scan_type][scan_type2];

			delr2[0] = scan_position[0] - inner_scan_position[0];
			delr2[1] = scan_position[1] - inner_scan_position[1];
			delr2[2] = scan_position[2] - inner_scan_position[2];
			rsq2 = delr2[0] * delr2[0] + delr2[1] * delr2[1] + delr2[2] * delr2[2];
			if (rsq2 >= params[ikparam].cutsq) continue;


			threebody(&params[ijparam], &params[ikparam], &params[ijkparam],
				rsq1, rsq2, delr1, delr2, fj, fk, quad_eflag, energy_contribution);
            quadrature_energy += energy_contribution/3;

			force_densityx += fj[0];
			force_densityy += fj[1];
			force_densityz += fj[2];
			if(atom->CAC_virial){
	   	virial_density[0] += delr1[0]*fj[0] + delr2[0]*fk[0];
		  virial_density[1] += delr1[1]*fj[1] + delr2[1]*fk[1];
		  virial_density[2] += delr1[2]*fj[2] + delr2[2]*fk[2];
		  virial_density[3] += delr1[0]*fj[1] + delr2[0]*fk[1];
		  virial_density[4] += delr1[0]*fj[2] + delr2[0]*fk[2];
		  virial_density[5] += delr1[1]*fj[2] + delr2[1]*fk[2];
		  }
		}

	}

	if(atom->CAC_virial){
	   	virial_density[0] *= 0.3333333;
		  virial_density[1] *= 0.3333333;
		  virial_density[2] *= 0.3333333;
		  virial_density[3] *= 0.3333333;
		  virial_density[4] *= 0.3333333;
		  virial_density[5] *= 0.3333333;
		  }
//end of scanning loop
}
//------------------------------------------------------------------------
