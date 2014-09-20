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
#include "stdlib.h"
#include "compute_q6_atom.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "math_const.h"
#include <math.h>
#include <complex>

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

ComputeQ6Atom::ComputeQ6Atom(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if ((narg != 4)) error->all(FLERR,"Illegal compute q6/atom command");
  q6cutoff = force->numeric(FLERR,arg[3]);
  if (q6cutoff <= 0.0) error->all(FLERR,"q6 cutoff must be > 0.0");
  peratom_flag = 1;
  size_peratom_cols = 0;
  nmax = 0;
  q6 = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeQ6Atom::~ComputeQ6Atom()
{
  memory->destroy(q6);
}

/* ---------------------------------------------------------------------- */

void ComputeQ6Atom::init()
{
  if (force->pair == NULL)
    error->all(FLERR,"Compute q6/atom requires a pair style be defined");
    
  if ((q6cutoff > force->pair->cutforce) && comm->me==0)
    error->warning(FLERR,"Compute q6/atom cutoff > force cutoff, force cutoff used for q6 calculation");
  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"q6/atom") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning(FLERR,"More than one compute q6/atom");
    
  // need an occasional full neighbor list

  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->occasional = 1;
}

/* ---------------------------------------------------------------------- */

void ComputeQ6Atom::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */
void ComputeQ6Atom::compute_peratom()
{
  std::complex<double> q6m[7], iim(0,1.0);
  int i,j,ii,jj,n,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz,phi,costheta,sintheta,ssintheta,rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;
  
  invoked_peratom = update->ntimestep;

  // grow q6 array if necessary

  if (atom->nlocal > nmax) {
    memory->destroy(q6);
    nmax = atom->nmax;
    memory->create(q6,nmax,"q6/atom:q6");
    vector_atom = q6;
  }

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list->index);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  // compute q6 parameter for each atom in group
  // use full neighbor list

  double **x = atom->x;
  int *mask = atom->mask;
  double cutsq;
  if (force->pair->cutforce > q6cutoff) cutsq = q6cutoff*q6cutoff;
  else cutsq = force->pair->cutforce*force->pair->cutforce;
  
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    q6[i] = 0;
    if (mask[i] & groupbit) {
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      // loop over list of all neighbors within force cutoff
		
      n = 0;
      for (jj = 0; jj < 7; jj++) {
		  q6m[jj] = 0.0;
	  }
	  
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;

        delx = x[j][0] - xtmp;
        dely = x[j][1] - ytmp;
        delz = x[j][2] - ztmp;
        rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < cutsq) {
          n += 1;
          costheta = delz/sqrt(rsq);
          sintheta = sin(acos(costheta));
          ssintheta = 1.-costheta*costheta;
          if (delx > 0) phi = atan(dely/delx);
          else if (delx < 0) 
		    phi = MY_PI + atan(dely/delx);
          else if (dely > 0) 
            phi = MY_PI2;
          else 
            phi = -MY_PI2;

          q6m[0] += 1./32.*sqrt(13./MY_PI)*(231.*pow(costheta,6)-315.*pow(costheta,4)+105.*pow(costheta,2)-5.);
          q6m[1] += 1./16.*sqrt(273./MY_2PI)*sintheta*(33.*pow(costheta,5)-30.*pow(costheta,3)+5.*costheta)*exp(iim*phi);
          q6m[2] += 1./64.*sqrt(1365./MY_PI)*ssintheta*(33.*pow(costheta,4)-18.*pow(costheta,2)+1.)*exp(2.*iim*phi);
          q6m[3] += 1./32.*sqrt(1365./MY_PI)*sintheta*ssintheta*(11.*pow(costheta,3)-3.*costheta)*exp(3.*iim*phi);
          q6m[4] += 3./32.*sqrt(91./MY_2PI)*ssintheta*ssintheta*(11.*pow(costheta,2)-1.)*exp(4.*iim*phi);
          q6m[5] += 3./32.*sqrt(1001./MY_PI)*ssintheta*ssintheta*sintheta*costheta*exp(5.*iim*phi);
          q6m[6] += 1./64.*sqrt(3003./MY_PI)*ssintheta*ssintheta*ssintheta*exp(6.*iim*phi);
        }
      }
      
      q6m[0] /= (double) n;
      q6[i] += norm(q6m[0]);
      
      for (jj = 1; jj < 7; jj++) {
		  q6m[jj] /= (double) n;
		  q6[i] += norm(q6m[jj])*2.;
	  }
	  
	  q6[i] *= MY_4PI/13.;
	  q6[i] = sqrt(q6[i]);
    }
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeQ6Atom::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}
