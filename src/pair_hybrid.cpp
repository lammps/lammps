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

/* ----------------------------------------------------------------------
   Contributing authors: 
     James Fischer (High Performance Technologies, Inc)
     Vincent Natoli (Stone Ridge Technology)
     David Richie (Stone Ridge Technology)
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "ctype.h"
#include "pair_hybrid.h"
#include "atom.h"
#include "force.h"
#include "pair.h"
#include "neighbor.h"
#include "update.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define NEIGHEXTRA 10000

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

PairHybrid::PairHybrid(LAMMPS *lmp) : Pair(lmp)
{
  nstyles = 0;
}

/* ---------------------------------------------------------------------- */

PairHybrid::~PairHybrid()
{
  if (nstyles) {
    for (int m = 0; m < nstyles; m++) delete styles[m];
    delete [] styles;
    for (int m = 0; m < nstyles; m++) delete [] keywords[m];
    delete [] keywords;
  }

  if (allocated) {
    memory->destroy_2d_int_array(setflag);
    memory->destroy_2d_int_array(map);
    memory->destroy_2d_double_array(cutsq);

    delete [] nnlist;
    delete [] maxneigh;
    for (int m = 0; m < nstyles; m++) memory->sfree(nlist[m]);
    delete [] nlist;

    delete [] nnlist_full;
    delete [] maxneigh_full;
    for (int m = 0; m < nstyles; m++) memory->sfree(nlist_full[m]);
    delete [] nlist_full;

    for (int m = 0; m < nstyles; m++) delete [] firstneigh[m];
    delete [] firstneigh;
    for (int m = 0; m < nstyles; m++) delete [] numneigh[m];
    delete [] numneigh;

    for (int m = 0; m < nstyles; m++) delete [] firstneigh_full[m];
    delete [] firstneigh_full;
    for (int m = 0; m < nstyles; m++) delete [] numneigh_full[m];
    delete [] numneigh_full;
  }
}

/* ---------------------------------------------------------------------- */

void PairHybrid::compute(int eflag, int vflag)
{
  int i,j,k,m,n,jfull,nneigh;
  int *neighs,*mapi;
  double **f_original;

  // save ptrs to original neighbor lists

  int **firstneigh_original = neighbor->firstneigh;
  int *numneigh_original = neighbor->numneigh;
  int **firstneigh_full_original = neighbor->firstneigh_full;
  int *numneigh_full_original = neighbor->numneigh_full;

  // if this is re-neighbor step, create sub-style lists

  if (neighbor->ago == 0) {
    int *type = atom->type;
    int nlocal = atom->nlocal;
    int nall = atom->nlocal + atom->nghost;

    // realloc per-atom per-style firstneigh/numneigh half/full if necessary

    if (nlocal > maxlocal) {
      maxlocal = nlocal;
      if (neigh_half_every) {
	for (m = 0; m < nstyles; m++) {
	  delete [] firstneigh[m];
	  delete [] numneigh[m];
	}
	for (m = 0; m < nstyles; m++) {
	  firstneigh[m] = new int*[maxlocal];
	  numneigh[m] = new int[maxlocal];
	}
      }
      if (neigh_full_every) {
	for (m = 0; m < nstyles; m++) {
	  delete [] firstneigh_full[m];
	  delete [] numneigh_full[m];
	}
	for (m = 0; m < nstyles; m++) {
	  firstneigh_full[m] = new int*[maxlocal];
	  numneigh_full[m] = new int[maxlocal];
	}
      }
    }

    // nnlist[] = length of each sub-style half list
    // nnlist_full[] = length of each sub-style full list
    // count from half and/or full list depending on what sub-styles use

    for (m = 0; m < nstyles; m++) nnlist[m] = nnlist_full[m] = 0;

    if (neigh_half_every && neigh_full_every) {
      for (i = 0; i < nlocal; i++) {
	mapi = map[type[i]];
	neighs = firstneigh_original[i];
	nneigh = numneigh_original[i];
	for (k = 0; k < nneigh; k++) {
	  j = neighs[k];
	  if (j >= nall) j %= nall;
	  m = mapi[type[j]];
	  if (styles[m] && styles[m]->neigh_half_every) nnlist[m]++;
	}
	neighs = firstneigh_full_original[i];
	nneigh = numneigh_full_original[i];
	for (k = 0; k < nneigh; k++) {
	  j = neighs[k];
	  if (j >= nall) j %= nall;
	  m = mapi[type[j]];
	  if (styles[m] && styles[m]->neigh_full_every) nnlist_full[m]++;
	}
      }

    } else if (neigh_half_every) {
      for (i = 0; i < nlocal; i++) {
	mapi = map[type[i]];
	neighs = firstneigh_original[i];
	nneigh = numneigh_original[i];
	for (k = 0; k < nneigh; k++) {
	  j = neighs[k];
	  if (j >= nall) j %= nall;
	  nnlist[mapi[type[j]]]++;
	}
      }

    } else if (neigh_full_every) {
      for (i = 0; i < nlocal; i++) {
	mapi = map[type[i]];
	neighs = firstneigh_full_original[i];
	nneigh = numneigh_full_original[i];
	for (k = 0; k < nneigh; k++) {
	  j = neighs[k];
	  if (j >= nall) j %= nall;
	  nnlist_full[mapi[type[j]]]++;
	}
      }
    }

    // realloc sub-style nlist and nlist_full if necessary

    if (neigh_half_every) {
      for (m = 0; m < nstyles; m++) {
	if (nnlist[m] > maxneigh[m]) {
	  memory->sfree(nlist[m]);
	  maxneigh[m] = nnlist[m] + NEIGHEXTRA;
	  nlist[m] = (int *)
	    memory->smalloc(maxneigh[m]*sizeof(int),"pair_hybrid:nlist");
	}
	nnlist[m] = 0;
      }
    }
    if (neigh_full_every) {
      for (m = 0; m < nstyles; m++) {
	if (nnlist_full[m] > maxneigh_full[m]) {
	  memory->sfree(nlist_full[m]);
	  maxneigh_full[m] = nnlist_full[m] + NEIGHEXTRA;
	  nlist_full[m] = (int *)
	    memory->smalloc(maxneigh_full[m]*sizeof(int),
			    "pair_hybrid:nlist_full");
	}
	nnlist_full[m] = 0;
      }
    }

    // load sub-style half/full list with values from original lists
    // load from half and/or full list depending on what sub-styles use

    if (neigh_half_every && neigh_full_every) {
      for (i = 0; i < nlocal; i++) {
	for (m = 0; m < nstyles; m++) {
	  firstneigh[m][i] = &nlist[m][nnlist[m]];
	  numneigh[m][i] = nnlist[m];
	  firstneigh_full[m][i] = &nlist_full[m][nnlist_full[m]];
	  numneigh_full[m][i] = nnlist_full[m];
	}
	mapi = map[type[i]];
	neighs = firstneigh_original[i];
	nneigh = numneigh_original[i];
	for (k = 0; k < nneigh; k++) {
	  j = jfull = neighs[k];
	  if (j >= nall) j %= nall;
	  m = mapi[type[j]];
	  if (styles[m] && styles[m]->neigh_half_every)
	    nlist[m][nnlist[m]++] = jfull;
	}
	neighs = firstneigh_full_original[i];
	nneigh = numneigh_full_original[i];
	for (k = 0; k < nneigh; k++) {
	  j = jfull = neighs[k];
	  if (j >= nall) j %= nall;
	  m = mapi[type[j]];
	  if (styles[m] && styles[m]->neigh_full_every)
	    nlist_full[m][nnlist_full[m]++] = jfull;
	}
	for (m = 0; m < nstyles; m++) {
	  numneigh[m][i] = nnlist[m] - numneigh[m][i];
	  numneigh_full[m][i] = nnlist_full[m] - numneigh_full[m][i];
	}
      }

    } else if (neigh_half_every) {
      for (i = 0; i < nlocal; i++) {
	for (m = 0; m < nstyles; m++) {
	  firstneigh[m][i] = &nlist[m][nnlist[m]];
	  numneigh[m][i] = nnlist[m];
	}
	mapi = map[type[i]];
	neighs = firstneigh_original[i];
	nneigh = numneigh_original[i];
	for (k = 0; k < nneigh; k++) {
	  j = jfull = neighs[k];
	  if (j >= nall) j %= nall;
	  m = mapi[type[j]];
	  nlist[m][nnlist[m]++] = jfull;
	}
	for (m = 0; m < nstyles; m++)
	  numneigh[m][i] = nnlist[m] - numneigh[m][i];
      }

    } else if (neigh_full_every) {
      for (i = 0; i < nlocal; i++) {
	for (m = 0; m < nstyles; m++) {
	  firstneigh_full[m][i] = &nlist_full[m][nnlist_full[m]];
	  numneigh_full[m][i] = nnlist_full[m];
	}
	mapi = map[type[i]];
	neighs = firstneigh_full_original[i];
	nneigh = numneigh_full_original[i];
	for (k = 0; k < nneigh; k++) {
	  j = jfull = neighs[k];
	  if (j >= nall) j %= nall;
	  m = mapi[type[j]];
	  nlist_full[m][nnlist_full[m]++] = jfull;
	}
	for (m = 0; m < nstyles; m++)
	  numneigh_full[m][i] = nnlist_full[m] - numneigh_full[m][i];
      }
    }
  }
  
  // call each sub-style's compute function
  // set neighbor->firstneigh/numneigh to sub-style lists before call
  //   set half or full or both depending on what sub-style uses
  // for vflag = 1:
  //   sub-style accumulates in its virial[6]
  //   sum sub-style virial[6] to hybrid's virial[6]
  // for vflag = 2:
  //   set atom->f to update->f_pair so sub-style will sum its f to f_pair
  //   call sub-style compute() with vflag % 2 to prevent sub-style
  //     from calling virial_compute()
  //   reset atom->f to stored f_original
  //   call hybrid virial_compute() which will use update->f_pair
  // accumulate sub-style energy,virial in hybrid's energy,virial

  eng_vdwl = eng_coul = 0.0;
  if (vflag) for (n = 0; n < 6; n++) virial[n] = 0.0;

  if (vflag == 2) {
    f_original = atom->f;
    atom->f = update->f_pair;
  }

  for (m = 0; m < nstyles; m++) {
    if (styles[m] == NULL) continue;
    if (styles[m]->neigh_half_every) {
      neighbor->firstneigh = firstneigh[m];
      neighbor->numneigh = numneigh[m];
    } 
    if (styles[m]->neigh_full_every) {
      neighbor->firstneigh_full = firstneigh_full[m];
      neighbor->numneigh_full = numneigh_full[m];
    }
    styles[m]->compute(eflag,vflag % 2);
    if (eflag) {
      eng_vdwl += styles[m]->eng_vdwl;
      eng_coul += styles[m]->eng_coul;
    }
    if (vflag == 1) for (n = 0; n < 6; n++) virial[n] += styles[m]->virial[n];
  }

  if (vflag == 2) {
    atom->f = f_original;
    virial_compute();
  }

  // restore ptrs to original neighbor lists

  neighbor->firstneigh = firstneigh_original;
  neighbor->numneigh = numneigh_original;
  neighbor->firstneigh_full = firstneigh_full_original;
  neighbor->numneigh_full = numneigh_full_original;
}

/* ----------------------------------------------------------------------
   allocate all arrays 
------------------------------------------------------------------------- */

void PairHybrid::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  map = memory->create_2d_int_array(n+1,n+1,"pair:map");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      map[i][j] = -1;

  setflag = memory->create_2d_int_array(n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  cutsq = memory->create_2d_double_array(n+1,n+1,"pair:cutsq");

  nnlist = new int[nstyles];
  maxneigh = new int[nstyles];
  nlist = new int*[nstyles];
  for (int m = 0; m < nstyles; m++) maxneigh[m] = 0;
  for (int m = 0; m < nstyles; m++) nlist[m] = NULL;

  nnlist_full = new int[nstyles];
  maxneigh_full = new int[nstyles];
  nlist_full = new int*[nstyles];
  for (int m = 0; m < nstyles; m++) maxneigh_full[m] = 0;
  for (int m = 0; m < nstyles; m++) nlist_full[m] = NULL;

  maxlocal = 0;
  firstneigh = new int**[nstyles];
  numneigh = new int*[nstyles];
  for (int m = 0; m < nstyles; m++) firstneigh[m] = NULL;
  for (int m = 0; m < nstyles; m++) numneigh[m] = NULL;
  firstneigh_full = new int**[nstyles];
  numneigh_full = new int*[nstyles];
  for (int m = 0; m < nstyles; m++) firstneigh_full[m] = NULL;
  for (int m = 0; m < nstyles; m++) numneigh_full[m] = NULL;
}

/* ----------------------------------------------------------------------
   create one pair style for each arg in list
------------------------------------------------------------------------- */

void PairHybrid::settings(int narg, char **arg)
{
  int i,m,istyle;

  if (narg < 1) error->all("Illegal pair_style command");

  // delete old lists, since cannot just change settings

  if (nstyles) {
    for (m = 0; m < nstyles; m++) delete styles[m];
    delete [] styles;
    for (m = 0; m < nstyles; m++) delete [] keywords[m];
    delete [] keywords;
  }

  if (allocated) {
    memory->destroy_2d_int_array(setflag);
    memory->destroy_2d_int_array(map);
    memory->destroy_2d_double_array(cutsq);
  }
  allocated = 0;

  // count sub-styles by skipping numeric args
  // one exception is 1st arg of style "table", which is non-numeric word

  nstyles = i = 0;
  while (i < narg) {
    if (strcmp(arg[i],"table") == 0) i++;
    i++;
    while (i < narg && !isalpha(arg[i][0])) i++;
    nstyles++;
  }

  // allocate list of sub-styles

  styles = new Pair*[nstyles];
  keywords = new char*[nstyles];

  // allocate each sub-style and call its settings() with subset of args
  // define subset of sub-styles by skipping numeric args
  // one exception is 1st arg of style "table", which is non-numeric word

  nstyles = i = 0;
  while (i < narg) {
    for (m = 0; m < nstyles; m++)
      if (strcmp(arg[i],keywords[m]) == 0) 
	error->all("Pair style hybrid cannot use same pair style twice");
    if (strcmp(arg[i],"hybrid") == 0) 
      error->all("Pair style hybrid cannot have hybrid as an argument");
    styles[nstyles] = force->new_pair(arg[i]);
    keywords[nstyles] = new char[strlen(arg[i])+1];
    strcpy(keywords[nstyles],arg[i]);
    istyle = i;
    if (strcmp(arg[i],"table") == 0) i++;
    i++;
    while (i < narg && !isalpha(arg[i][0])) i++;
    if (styles[nstyles]) styles[nstyles]->settings(i-istyle-1,&arg[istyle+1]);
    nstyles++;
  }

  // set comm_forward, comm_reverse to max of any sub-style

  for (m = 0; m < nstyles; m++) {
    if (styles[m]) comm_forward = MAX(comm_forward,styles[m]->comm_forward);
    if (styles[m]) comm_reverse = MAX(comm_reverse,styles[m]->comm_reverse);
  }

  // neigh_every = 1 if any sub-style = 1

  neigh_half_every = neigh_full_every = 0;
  for (m = 0; m < nstyles; m++) {
    if (styles[m] && styles[m]->neigh_half_every) neigh_half_every = 1;
    if (styles[m] && styles[m]->neigh_full_every) neigh_full_every = 1;
  }

  // single_enable = 0 if any sub-style = 0

  for (m = 0; m < nstyles; m++)
    if (styles[m] && styles[m]->single_enable == 0) single_enable = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairHybrid::coeff(int narg, char **arg)
{
  if (narg < 3) error->all("Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  // 3rd arg = pair style name

  int m;
  for (m = 0; m < nstyles; m++)
    if (strcmp(arg[2],keywords[m]) == 0) break;
  if (m == nstyles) error->all("Pair coeff for hybrid has invalid style");

  // move 1st/2nd args to 2nd/3rd args

  sprintf(arg[2],"%s",arg[1]);
  sprintf(arg[1],"%s",arg[0]);

  // invoke sub-style coeff() starting with 1st arg

  if (styles[m]) styles[m]->coeff(narg-1,&arg[1]);

  // if pair style only allows one pair coeff call (with * * and type mapping)
  // then unset any old setflag/map assigned to that style first
  // in case pair coeff for this sub-style is being called for 2nd time

  if (styles[m] && styles[m]->one_coeff)
    for (int i = 1; i <= atom->ntypes; i++)
      for (int j = i; j <= atom->ntypes; j++)
	if (map[i][j] == m) {
	  map[i][j] = -1;
	  setflag[i][j] = 0;
	}

  // set hybrid map & setflag only if substyle set its setflag
  // if sub-style is NULL for "none", still set setflag

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      if (styles[m] == NULL || styles[m]->setflag[i][j]) {
	map[i][j] = m;
	setflag[i][j] = 1;
	count++;
      }
    }
  }

  if (count == 0) error->all("Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairHybrid::init_one(int i, int j)
{
  // if i,j is set explicity, call its sub-style
  // if i,j is not set and i,i sub-style = j,j sub-style
  // then set map[i][j] to this sub-style and call sub-style for init/mixing
  // else i,j has not been set by user
  // check for special case = style none

  double cut = 0.0;
  if (setflag[i][j]) {
    if (styles[map[i][j]]) {
      cut = styles[map[i][j]]->init_one(i,j);
      styles[map[i][j]]->cutsq[i][j] = 
	styles[map[i][j]]->cutsq[j][i] = cut*cut;
      if (tail_flag) {
	etail_ij = styles[map[i][j]]->etail_ij;
	ptail_ij = styles[map[i][j]]->ptail_ij;
      }
    }
  } else if (map[i][i] == map[j][j]) {
    map[i][j] = map[i][i];
    if (styles[map[i][j]]) {
      cut = styles[map[i][j]]->init_one(i,j);
      styles[map[i][j]]->cutsq[i][j] = 
	styles[map[i][j]]->cutsq[j][i] = cut*cut;
      if (tail_flag) {
	etail_ij = styles[map[i][j]]->etail_ij;
	ptail_ij = styles[map[i][j]]->ptail_ij;
      }
    }
  } else error->all("All pair coeffs are not set");

  map[j][i] = map[i][j];

  return cut;
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairHybrid::init_style()
{
  for (int m = 0; m < nstyles; m++)
    if (styles[m]) styles[m]->init_style();
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairHybrid::write_restart(FILE *fp)
{
  fwrite(&nstyles,sizeof(int),1,fp);

  // each sub-style writes its settings

  int n;
  for (int m = 0; m < nstyles; m++) {
    n = strlen(keywords[m]) + 1;
    fwrite(&n,sizeof(int),1,fp);
    fwrite(keywords[m],sizeof(char),n,fp);
    if (styles[m]) styles[m]->write_restart_settings(fp);
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairHybrid::read_restart(FILE *fp)
{
  allocate();

  int me = comm->me;
  if (me == 0) fread(&nstyles,sizeof(int),1,fp);
  MPI_Bcast(&nstyles,1,MPI_INT,0,world);

  styles = new Pair*[nstyles];
  keywords = new char*[nstyles];
  
  // each sub-style is created via new_pair() and reads its settings

  int n;
  for (int m = 0; m < nstyles; m++) {
    if (me == 0) fread(&n,sizeof(int),1,fp);
    MPI_Bcast(&n,1,MPI_INT,0,world);
    keywords[m] = new char[n];
    if (me == 0) fread(keywords[m],sizeof(char),n,fp);
    MPI_Bcast(keywords[m],n,MPI_CHAR,0,world);
    styles[m] = force->new_pair(keywords[m]);
    if (styles[m]) styles[m]->read_restart_settings(fp);
  }
}

/* ---------------------------------------------------------------------- */

void PairHybrid::single(int i, int j, int itype, int jtype,
			double rsq, double factor_coul, double factor_lj,
			int eflag, One &one)
{
  if (map[itype][jtype] == -1)
    error->one("Invoked pair single on pair style none");

  styles[map[itype][jtype]]->
    single(i,j,itype,jtype,rsq,factor_coul,factor_lj,eflag,one);
}

/* ---------------------------------------------------------------------- */

void PairHybrid::single_embed(int i, int itype, double &phi)
{
  if (map[itype][itype] == -1)
    error->one("Invoked pair single on pair style none");
  
  styles[map[itype][itype]]->single_embed(i,itype,phi);
}

/* ----------------------------------------------------------------------
   modify parameters of the pair style
   call modify_params of PairHybrid
   also pass command args to each sub-style of hybrid
------------------------------------------------------------------------- */

void PairHybrid::modify_params(int narg, char **arg)
{
  Pair::modify_params(narg,arg);
  for (int m = 0; m < nstyles; m++)
    if (styles[m]) styles[m]->modify_params(narg,arg);
}

/* ----------------------------------------------------------------------
   memory usage of sub-style firstneigh, numneigh, neighbor list
   add in memory usage of each sub-style itself
------------------------------------------------------------------------- */

int PairHybrid::memory_usage()
{
  int bytes = nstyles*maxlocal * (sizeof(int *) + sizeof(int));
  for (int m = 0; m < nstyles; m++) bytes += maxneigh[m] * sizeof(int);
  for (int m = 0; m < nstyles; m++)
    if (styles[m]) bytes += styles[m]->memory_usage();
  return bytes;
}
