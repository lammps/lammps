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

#include "math.h"
#include "string.h"
#include "ctype.h"
#include "pair_hybrid.h"
#include "atom.h"
#include "force.h"
#include "pair.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "update.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

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
    memory->destroy_2d_double_array(cutsq);
    memory->destroy_2d_int_array(nmap);
    memory->destroy_3d_int_array(map);
  }
}

/* ----------------------------------------------------------------------
  call each sub-style's compute function
  accumulate sub-style energy/virial in hybrid's energy/virial
  for vflag = 1:
    each sub-style computes own virial[6]
    sum sub-style virial[6] to hybrid's virial[6]
  for vflag = 2:
    call sub-style compute() with vflag % 2
      to prevent sub-style from calling virial_compute()
    hybrid calls virial_compute() on final accumulated f
------------------------------------------------------------------------- */

void PairHybrid::compute(int eflag, int vflag)
{
  int m,n;

  eng_vdwl = eng_coul = 0.0;
  if (vflag) for (n = 0; n < 6; n++) virial[n] = 0.0;

  for (m = 0; m < nstyles; m++) {
    styles[m]->compute(eflag,vflag % 2);
    if (eflag) {
      eng_vdwl += styles[m]->eng_vdwl;
      eng_coul += styles[m]->eng_coul;
    }
    if (vflag == 1) for (n = 0; n < 6; n++) virial[n] += styles[m]->virial[n];
  }

  if (vflag == 2) virial_compute();
}

/* ---------------------------------------------------------------------- */

void PairHybrid::compute_inner()
{
  for (int m = 0; m < nstyles; m++)
    if (styles[m]->respa_enable) styles[m]->compute_inner();
}

/* ---------------------------------------------------------------------- */

void PairHybrid::compute_middle()
{
  for (int m = 0; m < nstyles; m++)
    if (styles[m]->respa_enable) styles[m]->compute_middle();
}

/* ---------------------------------------------------------------------- */

void PairHybrid::compute_outer(int eflag, int vflag)
{
  for (int m = 0; m < nstyles; m++)
    if (styles[m]->respa_enable) styles[m]->compute_outer(eflag,vflag);
}

/* ----------------------------------------------------------------------
   allocate all arrays 
------------------------------------------------------------------------- */

void PairHybrid::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  setflag = memory->create_2d_int_array(n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  cutsq = memory->create_2d_double_array(n+1,n+1,"pair:cutsq");

  nmap = memory->create_2d_int_array(n+1,n+1,"pair:nmap");
  map = memory->create_3d_int_array(n+1,n+1,nstyles,"pair:map");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      nmap[i][j] = 0;
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
    memory->destroy_2d_double_array(cutsq);
    memory->destroy_2d_int_array(nmap);
    memory->destroy_3d_int_array(map);
  }
  allocated = 0;

  // count sub-styles by skipping numeric args
  // one exception is 1st arg of style "table", which is non-numeric word

  nstyles = 0;
  i = 0;
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
  // define subset of args for a sub-style by skipping numeric args
  // one exception is 1st arg of style "table", which is non-numeric word

  nstyles = 0;
  i = 0;
  while (i < narg) {
    for (m = 0; m < nstyles; m++)
      if (strcmp(arg[i],keywords[m]) == 0) 
	error->all("Pair style hybrid cannot use same pair style twice");
    if (strcmp(arg[i],"hybrid") == 0) 
      error->all("Pair style hybrid cannot have hybrid as an argument");
    if (strcmp(arg[i],"none") == 0) 
      error->all("Pair style hybrid cannot have none as an argument");
    styles[nstyles] = force->new_pair(arg[i]);
    keywords[nstyles] = new char[strlen(arg[i])+1];
    strcpy(keywords[nstyles],arg[i]);
    istyle = i;
    if (strcmp(arg[i],"table") == 0) i++;
    i++;
    while (i < narg && !isalpha(arg[i][0])) i++;
    styles[nstyles]->settings(i-istyle-1,&arg[istyle+1]);
    nstyles++;
  }

  // set comm_forward, comm_reverse to max of any sub-style

  for (m = 0; m < nstyles; m++) {
    if (styles[m]) comm_forward = MAX(comm_forward,styles[m]->comm_forward);
    if (styles[m]) comm_reverse = MAX(comm_reverse,styles[m]->comm_reverse);
  }

  // single_enable = 0 if any sub-style = 0
  // respa_enable = 1 if any sub-style is set

  for (m = 0; m < nstyles; m++)
    if (styles[m]->single_enable == 0) single_enable = 0;
  for (m = 0; m < nstyles; m++)
    if (styles[m]->respa_enable) respa_enable = 1;
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

  // 3rd arg = pair sub-style name
  // allow for "none" as valid sub-style name

  int m;
  for (m = 0; m < nstyles; m++)
    if (strcmp(arg[2],keywords[m]) == 0) break;

  int none = 0;
  if (m == nstyles) {
    if (strcmp(arg[2],"none") == 0) none = 1;
    else error->all("Pair coeff for hybrid has invalid style");
  }

  // move 1st/2nd args to 2nd/3rd args

  sprintf(arg[2],"%s",arg[1]);
  sprintf(arg[1],"%s",arg[0]);

  // invoke sub-style coeff() starting with 1st arg

  if (!none) styles[m]->coeff(narg-1,&arg[1]);

  // if sub-style only allows one pair coeff call (with * * and type mapping)
  // then unset setflag/map assigned to that style before setting it below
  // in case pair coeff for this sub-style is being called for 2nd time

  if (!none && styles[m]->one_coeff)
    for (int i = 1; i <= atom->ntypes; i++)
      for (int j = i; j <= atom->ntypes; j++)
	if (nmap[i][j] && map[i][j][0] == m) {
	  setflag[i][j] = 0;
	  nmap[i][j] = 0;
	}

  // set setflag and which type pairs map to which sub-style
  // if sub-style is none: set hybrid subflag, wipe out map
  // else: set hybrid setflag & map only if substyle setflag is set
  //   previous mappings are wiped out

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      if (none) {
	setflag[i][j] = 1;
	nmap[i][j] = 0;
	count++;
      } else if (styles[m]->setflag[i][j]) {
	setflag[i][j] = 1;
	nmap[i][j] = 1;
	map[i][j][0] = m;
	count++;
      }
    }
  }

  if (count == 0) error->all("Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairHybrid::init_style()
{
  int i,m,itype,jtype,used,istyle,skip;

  // error if a sub-style is not used

  int ntypes = atom->ntypes;

  for (istyle = 0; istyle < nstyles; istyle++) {
    used = 0;
    for (itype = 1; itype <= ntypes; itype++)
      for (jtype = itype; jtype <= ntypes; jtype++)
	for (m = 0; m < nmap[itype][jtype]; m++)
	  if (map[itype][jtype][m] == istyle) used = 1;
    if (used == 0) error->all("Pair hybrid sub-style is not used");
  }

  // each sub-style makes its neighbor list request(s)

  for (istyle = 0; istyle < nstyles; istyle++) styles[istyle]->init_style();

  // create skip lists for each neigh request

  for (i = 0; i < neighbor->nrequest; i++) {

    // only relevant for half, full, gran, respaouter lists

    if (neighbor->requests[i]->granhistory) continue;
    if (neighbor->requests[i]->respamiddle) continue;
    if (neighbor->requests[i]->respainner) continue;
    if (neighbor->requests[i]->half_from_full) continue;

    // find associated sub-style

    for (istyle = 0; istyle < nstyles; istyle++)
      if (styles[istyle] == neighbor->requests[i]->requestor) break;

    // allocate iskip and ijskip
    // initialize so as to skip all pair types
    // set ijskip = 0 if type pair matches any entry in sub-style map
    // set ijskip = 0 if mixing will assign type pair to this sub-style
    //   will occur if both I,I and J,J are assigned to single sub-style
    //   and sub-style for both I,I and J,J match istyle
    // set iskip = 1 only if all ijskip for itype are 1

    int *iskip = new int[ntypes+1];
    int **ijskip = memory->create_2d_int_array(ntypes+1,ntypes+1,
					       "pair_hybrid:ijskip");

    for (itype = 1; itype <= ntypes; itype++)
      for (jtype = 1; jtype <= ntypes; jtype++)
	ijskip[itype][jtype] = 1;

    for (itype = 1; itype <= ntypes; itype++)
      for (jtype = itype; jtype <= ntypes; jtype++) {
	for (m = 0; m < nmap[itype][jtype]; m++)
	  if (map[itype][jtype][m] == istyle) 
	    ijskip[itype][jtype] = ijskip[jtype][itype] = 0;
	if (nmap[itype][itype] == 1 && map[itype][itype][0] == istyle &&
	    nmap[jtype][jtype] == 1 && map[jtype][jtype][0] == istyle)
	  ijskip[itype][jtype] = ijskip[jtype][itype] = 0;
      }

    for (itype = 1; itype <= ntypes; itype++) {
      iskip[itype] = 1;
      for (jtype = 1; jtype <= ntypes; jtype++)
	if (ijskip[itype][jtype] == 0) iskip[itype] = 0;
    }

    // if any skipping occurs
    // set request->skip and copy iskip and ijskip into request
    // else delete iskip and ijskip

    skip = 0;
    for (itype = 1; itype <= ntypes; itype++)
      for (jtype = 1; jtype <= ntypes; jtype++)
	if (ijskip[itype][jtype] == 1) skip = 1;

    if (skip) {
      neighbor->requests[i]->skip = 1;
      neighbor->requests[i]->iskip = iskip;
      neighbor->requests[i]->ijskip = ijskip;
    } else {
      delete [] iskip;
      memory->destroy_2d_int_array(ijskip);
    }
  }

  // combine sub-style neigh list requests and create new ones if needed

  modify_requests();
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairHybrid::init_one(int i, int j)
{
  // if I,J is not set explicitly:
  // perform mixing only if I,I sub-style = J,J sub-style
  // also require I,I and J,J both are assigned to single sub-style

  if (setflag[i][j] == 0) {
    if (nmap[i][i] != 1 || nmap[j][j] != 1 || map[i][i][0] != map[j][j][0])
      error->one("All pair coeffs are not set");
    nmap[i][j] = 1;
    map[i][j][0] = map[i][i][0];
  }

  // call init/mixing for all sub-styles of I,J
  // set cutsq in sub-style just as pair::init_one() does
  // sum tail corrections for I,J
  // return max cutoff of all sub-styles assigned to I,J
  // if no sub-styles assigned to I,J (pair_coeff none), cutmax = 0.0 returned

  double cutmax = 0.0;
  if (tail_flag) etail_ij = ptail_ij = 0.0;

  for (int k = 0; k < nmap[i][j]; k++) {
    map[j][i][k] = map[i][j][k];
    double cut = styles[map[i][j][k]]->init_one(i,j);
    styles[map[i][j][k]]->cutsq[i][j] = 
      styles[map[i][j][k]]->cutsq[j][i] = cut*cut;
    if (tail_flag) {
      etail_ij += styles[map[i][j][k]]->etail_ij;
      ptail_ij += styles[map[i][j][k]]->ptail_ij;
    }
    cutmax = MAX(cutmax,cut);
  }

  return cutmax;
}

/* ----------------------------------------------------------------------
   combine sub-style neigh list requests and create new ones if needed
------------------------------------------------------------------------- */

void PairHybrid::modify_requests()
{
  int i,j;

  // if list is skip list, look for non-skip list of same kind
  // if one exists, point at that one
  // else make new non-skip request of same kind and point at that one

  for (i = 0; i < neighbor->nrequest; i++) {
    if (neighbor->requests[i]->skip == 0) continue;

    for (j = 0; j < neighbor->nrequest; j++)
      if (neighbor->requests[i]->same_kind(neighbor->requests[j]) &&
	  neighbor->requests[j]->skip == 0) break;

    if (j < neighbor->nrequest) neighbor->requests[i]->otherlist = j;
    else {
      int irequest = neighbor->request(this);
      neighbor->requests[irequest]->copy_kind(neighbor->requests[i]);
      neighbor->requests[i]->otherlist = irequest;
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairHybrid::write_restart(FILE *fp)
{
  fwrite(&nstyles,sizeof(int),1,fp);

  // each sub-style writes its settings, but no coeff info

  int n;
  for (int m = 0; m < nstyles; m++) {
    n = strlen(keywords[m]) + 1;
    fwrite(&n,sizeof(int),1,fp);
    fwrite(keywords[m],sizeof(char),n,fp);
    styles[m]->write_restart_settings(fp);
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairHybrid::read_restart(FILE *fp)
{
  int me = comm->me;
  if (me == 0) fread(&nstyles,sizeof(int),1,fp);
  MPI_Bcast(&nstyles,1,MPI_INT,0,world);

  styles = new Pair*[nstyles];
  keywords = new char*[nstyles];
  
  // each sub-style is created via new_pair()
  // each reads its settings, but no coeff info

  int n;
  for (int m = 0; m < nstyles; m++) {
    if (me == 0) fread(&n,sizeof(int),1,fp);
    MPI_Bcast(&n,1,MPI_INT,0,world);
    keywords[m] = new char[n];
    if (me == 0) fread(keywords[m],sizeof(char),n,fp);
    MPI_Bcast(keywords[m],n,MPI_CHAR,0,world);
    styles[m] = force->new_pair(keywords[m]);
    styles[m]->read_restart_settings(fp);
  }
}

/* ----------------------------------------------------------------------
   call sub-style to compute single interaction
   since overlay could have multiple sub-styles, sum results explicitly
------------------------------------------------------------------------- */

void PairHybrid::single(int i, int j, int itype, int jtype,
			double rsq, double factor_coul, double factor_lj,
			int eflag, One &one)
{
  if (nmap[itype][jtype] == 0)
    error->one("Invoked pair single on pair style none");

  double fforce = 0.0;
  double eng_vdwl = 0.0;
  double eng_coul = 0.0;

  for (int m = 0; m < nmap[itype][jtype]; m++) {
    styles[map[itype][jtype][m]]->
      single(i,j,itype,jtype,rsq,factor_coul,factor_lj,eflag,one);
    fforce += one.fforce;
    eng_vdwl += one.eng_vdwl;
    eng_coul += one.eng_coul;
  }

  one.fforce = fforce;
  one.eng_vdwl = eng_vdwl;
  one.eng_coul = eng_coul;
}

/* ----------------------------------------------------------------------
   call sub-style to compute single embedding interaction
   since overlay could have multiple sub-styles, sum results explicitly
------------------------------------------------------------------------- */

void PairHybrid::single_embed(int i, int itype, double &phi)
{
  if (nmap[itype][itype] == 0)
    error->one("Invoked pair single on pair style none");
  
  phi = 0.0;
  double phi_single;

  for (int m = 0; m < nmap[itype][itype]; m++) {
    styles[map[itype][itype][m]]->single_embed(i,itype,phi_single);
    phi += phi_single;
  }
}

/* ----------------------------------------------------------------------
   modify parameters of the pair style
   call modify_params of PairHybrid
   also pass command args to each sub-style of hybrid
------------------------------------------------------------------------- */

void PairHybrid::modify_params(int narg, char **arg)
{
  Pair::modify_params(narg,arg);
  for (int m = 0; m < nstyles; m++) styles[m]->modify_params(narg,arg);
}

/* ----------------------------------------------------------------------
   memory usage of each sub-style
------------------------------------------------------------------------- */

double PairHybrid::memory_usage()
{
  double bytes = 0.0;
  for (int m = 0; m < nstyles; m++) bytes += styles[m]->memory_usage();
  return bytes;
}

/* ----------------------------------------------------------------------
   extract a ptr to a particular quantity stored by pair
   pass request thru to sub-styles
   return first non-NULL result except for cut_coul request
   for cut_coul, insure all non-NULL results are equal since required by Kspace
------------------------------------------------------------------------- */

void *PairHybrid::extract(char *str)
{
  void *cutptr = NULL;
  void *ptr;
  double cutvalue;

  for (int m = 0; m < nstyles; m++) {
    ptr = styles[m]->extract(str);
    if (ptr && strcmp(str,"cut_coul") == 0) {
      double *p_newvalue = (double *) ptr;
      double newvalue = *p_newvalue;
      if (cutptr && newvalue != cutvalue)
	error->all("Coulomb cutoffs of pair hybrid sub-styles do not match");
      cutptr = ptr;
      cutvalue = newvalue;
    } else if (ptr) return ptr;
  }

  if (strcmp(str,"cut_coul") == 0) return cutptr;
  return NULL;
}
