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

#include "string.h"
#include "compute_stress_atom.h"
#include "atom.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "modify.h"
#include "comm.h"
#include "update.h"
#include "force.h"
#include "bond.h"
#include "pair.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

ComputeStressAtom::ComputeStressAtom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg < 3) error->all("Illegal compute stress/atom command");

  peratom_flag = 1;
  size_peratom = 6;
  comm_reverse = 6;

  // process args

  kerequest = pairrequest = bondrequest = 1;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"ke") == 0) {
      if (iarg+2 > narg) error->all("Illegal compute ebond/atom command");
      if (strcmp(arg[iarg+1],"yes") == 0) kerequest = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) kerequest = 0;
      iarg += 2;
    } else if (strcmp(arg[iarg],"pair") == 0) {
      if (iarg+2 > narg) error->all("Illegal compute ebond/atom command");
      if (strcmp(arg[iarg+1],"yes") == 0) pairrequest = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) pairrequest = 0;
      iarg += 2;
    } else if (strcmp(arg[iarg],"bond") == 0) {
      if (iarg+2 > narg) error->all("Illegal compute ebond/atom command");
      if (strcmp(arg[iarg+1],"yes") == 0) bondrequest = 1;
      else if (strcmp(arg[iarg+1],"no") == 0) bondrequest = 0;
      iarg += 2;
    } else error->all("Illegal compute ebond/atom command");
  }

  nmax = 0;
  stress = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeStressAtom::~ComputeStressAtom()
{
  memory->destroy_2d_double_array(stress);
}

/* ---------------------------------------------------------------------- */

void ComputeStressAtom::init()
{
  if (pairrequest) {
    if (force->pair == NULL || force->pair->single_enable == 0)
      error->all("Pair style does not support computing per-atom stress");
    pairflag = 1;

    // need an occasional half neighbor list

    int irequest = neighbor->request((void *) this);
    neighbor->requests[irequest]->pair = 0;
    neighbor->requests[irequest]->compute = 1;
    neighbor->requests[irequest]->occasional = 1;

  } else pairflag = 0;

  if (bondrequest && force->bond) bondflag = 1;
  else bondflag = 0;

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"stress/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning("More than one compute stress/atom");
}

/* ---------------------------------------------------------------------- */

void ComputeStressAtom::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeStressAtom::compute_peratom()
{
  int i,j,ii,jj,n,i1,i2,inum,jnum,itype,jtype,iflag,jflag;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq,eng;
  double factor_coul,factor_lj,fforce,rmass;
  int *ilist,*jlist,*numneigh,**firstneigh;
  Pair::One one;

  // grow stress array if necessary

  if (atom->nmax > nmax) {
    memory->destroy_2d_double_array(stress);
    nmax = atom->nmax;
    stress =
      memory->create_2d_double_array(nmax,6,"compute/stress/atom:stress");
    vector_atom = stress;
  }

  // clear stress array
  // n includes ghosts only if newton flag is set

  int nlocal = atom->nlocal;

  if (force->newton) n = nlocal + atom->nghost;
  else n = nlocal;

  for (i = 0; i < n; i++) {
    stress[i][0] = 0.0;
    stress[i][1] = 0.0;
    stress[i][2] = 0.0;
    stress[i][3] = 0.0;
    stress[i][4] = 0.0;
    stress[i][5] = 0.0;
  }

  // compute pairwise stress for atoms via pair->single()
  // if neither atom is in compute group, skip that pair
  // only add stress to atoms in group

  if (pairflag) {

    // invoke half neighbor list (will copy or build if necessary)

    neighbor->build_one(list->index);

    inum = list->inum;
    ilist = list->ilist;
    numneigh = list->numneigh;
    firstneigh = list->firstneigh;

    double *special_coul = force->special_coul;
    double *special_lj = force->special_lj;
    double **cutsq = force->pair->cutsq;
    int newton_pair = force->newton_pair;

    double **x = atom->x;
    int *mask = atom->mask;
    int *type = atom->type;
    int nall = nlocal + atom->nghost;

    for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      itype = type[i];
      iflag = mask[i] & groupbit;
      jlist = firstneigh[i];
      jnum = numneigh[i];

      for (jj = 0; jj < jnum; jj++) {
	j = jlist[jj];
	jflag = mask[j] & groupbit;
	if (iflag == 0 && jflag == 0) continue;

	if (j < nall) factor_coul = factor_lj = 1.0;
	else {
	  factor_coul = special_coul[j/nall];
	  factor_lj = special_lj[j/nall];
	  j %= nall;
	}
	
	delx = xtmp - x[j][0];
	dely = ytmp - x[j][1];
	delz = ztmp - x[j][2];
	rsq = delx*delx + dely*dely + delz*delz;
	jtype = type[j];
	
	if (rsq < cutsq[itype][jtype]) {
	  force->pair->single(i,j,itype,jtype,rsq,factor_coul,factor_lj,0,one);
	  fforce = one.fforce;
	  
	  if (iflag) {
	    stress[i][0] -= delx*delx*fforce;
	    stress[i][1] -= dely*dely*fforce;
	    stress[i][2] -= delz*delz*fforce;
	    stress[i][3] -= delx*dely*fforce;
	    stress[i][4] -= delx*delz*fforce;
	    stress[i][5] -= dely*delz*fforce;
	  }
	  if (jflag && (newton_pair || j < nlocal)) {
	    stress[j][0] -= delx*delx*fforce;
	    stress[j][1] -= dely*dely*fforce;
	    stress[j][2] -= delz*delz*fforce;
	    stress[j][3] -= delx*dely*fforce;
	    stress[j][4] -= delx*delz*fforce;
	    stress[j][5] -= dely*delz*fforce;
	  }
	}
      }
    }
    
  }

  // compute bond stress for atoms via bond->single()
  // if neither atom is in compute group, skip that bond
  // only add stress to atoms in group

  if (bondflag) {
    double **x = atom->x;
    int *mask = atom->mask;
    int **bondlist = neighbor->bondlist;
    int nbondlist = neighbor->nbondlist;
    int newton_bond = force->newton_bond;
    int type;

    for (n = 0; n < nbondlist; n++) {
      i1 = bondlist[n][0];
      i2 = bondlist[n][1];
      type = bondlist[n][2];
      iflag = mask[i1] & groupbit;
      jflag = mask[i2] & groupbit;
      if (iflag == 0 && jflag == 0) continue;
      
      delx = x[i1][0] - x[i2][0];
      dely = x[i1][1] - x[i2][1];
      delz = x[i1][2] - x[i2][2];
      domain->minimum_image(delx,dely,delz);
      
      rsq = delx*delx + dely*dely + delz*delz;
      force->bond->single(type,rsq,i1,i2,0,fforce,eng);
      
      if (iflag) {
	stress[i1][0] -= delx*delx*fforce;
	stress[i1][1] -= dely*dely*fforce;
	stress[i1][2] -= delz*delz*fforce;
	stress[i1][3] -= delx*dely*fforce;
	stress[i1][4] -= delx*delz*fforce;
	stress[i1][5] -= dely*delz*fforce;
      }

      if (jflag && (newton_bond || i2 < nlocal)) {
	stress[i2][0] -= delx*delx*fforce;
	stress[i2][1] -= dely*dely*fforce;
	stress[i2][2] -= delz*delz*fforce;
	stress[i2][3] -= delx*dely*fforce;
	stress[i2][4] -= delx*delz*fforce;
	stress[i2][5] -= dely*delz*fforce;
      }
    }
  }

  // communicate stress between neighbor procs

  if (force->newton) comm->reverse_comm_compute(this);

  // remove double counting

  for (i = 0; i < nlocal; i++) {
    stress[i][0] *= 0.5;
    stress[i][1] *= 0.5;
    stress[i][2] *= 0.5;
    stress[i][3] *= 0.5;
    stress[i][4] *= 0.5;
    stress[i][5] *= 0.5;
  }

  // include kinetic energy term for each atom in group
  // mvv2e converts mv^2 to energy

  if (kerequest) {
    double **v = atom->v;
    double *mass = atom->mass;
    int *mask = atom->mask;
    int *type = atom->type;
    double mvv2e = force->mvv2e;
    
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	rmass = mvv2e * mass[type[i]];
	stress[i][0] -= rmass*v[i][0]*v[i][0];
	stress[i][1] -= rmass*v[i][1]*v[i][1];
	stress[i][2] -= rmass*v[i][2]*v[i][2];
	stress[i][3] -= rmass*v[i][0]*v[i][1];
	stress[i][4] -= rmass*v[i][0]*v[i][2];
	stress[i][5] -= rmass*v[i][1]*v[i][2];
      }
    }
  }

  // convert to pressure units (actually stress/volume = pressure)

  double nktv2p = force->nktv2p;
  for (i = 0; i < nlocal; i++) {
    stress[i][0] *= nktv2p;
    stress[i][1] *= nktv2p;
    stress[i][2] *= nktv2p;
    stress[i][3] *= nktv2p;
    stress[i][4] *= nktv2p;
    stress[i][5] *= nktv2p;
  }
}

/* ---------------------------------------------------------------------- */

int ComputeStressAtom::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = stress[i][0];
    buf[m++] = stress[i][1];
    buf[m++] = stress[i][2];
    buf[m++] = stress[i][3];
    buf[m++] = stress[i][4];
    buf[m++] = stress[i][5];
  }
  return 6;
}

/* ---------------------------------------------------------------------- */

void ComputeStressAtom::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    stress[j][0] += buf[m++];
    stress[j][1] += buf[m++];
    stress[j][2] += buf[m++];
    stress[j][3] += buf[m++];
    stress[j][4] += buf[m++];
    stress[j][5] += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeStressAtom::memory_usage()
{
  double bytes = nmax*6 * sizeof(double);
  return bytes;
}
