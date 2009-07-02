/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-93AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
 ------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Reese Jones (Sandia)
                         Philip Howell (Siemens)
                         Vikas Varsney (Air Force Research Laboratory)
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "compute_heat_flux.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "modify.h"
#include "group.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{DUMMY0,INVOKED_SCALAR,INVOKED_VECTOR,DUMMMY3,INVOKED_PERATOM};

/* ---------------------------------------------------------------------- */

ComputeHeatFlux::ComputeHeatFlux(LAMMPS *lmp, int narg, char **arg) : 
  Compute(lmp, narg, arg)
{
  if (narg != 4) error->all("Illegal compute heat/flux command");

  vector_flag = 1;
  size_vector = 6;
  extvector = 1;

  // store pe/atom ID used by heat flux computation
  // insure it is valid for pe/atom computation

  int n = strlen(arg[3]) + 1;
  id_atomPE = new char[n];
  strcpy(id_atomPE,arg[3]);

  int icompute = modify->find_compute(id_atomPE);
  if (icompute < 0) error->all("Could not find compute heat/flux compute ID");
  if (modify->compute[icompute]->peatomflag == 0)
    error->all("Compute heat/flux compute ID does not compute pe/atom");

  vector = new double[6];
}

/* ---------------------------------------------------------------------- */

ComputeHeatFlux::~ComputeHeatFlux()
{ 
  delete [] vector;
}

/* ---------------------------------------------------------------------- */

void ComputeHeatFlux::init()
{
  // error checks

  if (atom->avec->ghost_velocity == 0)
    error->all("Compute heat/flux requires ghost atoms store velocity");

  if (force->pair == NULL || force->pair->single_enable == 0)
    error->all("Pair style does not support compute heat/flux");

  int icompute = modify->find_compute(id_atomPE);
  if (icompute < 0) 
    error->all("Compute ID for compute heat/flux does not exist");
  atomPE = modify->compute[icompute];

  pair = force->pair;
  cutsq = force->pair->cutsq;

  // need an occasional half neighbor list

  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->occasional = 1;
}

/* ---------------------------------------------------------------------- */

void ComputeHeatFlux::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeHeatFlux::compute_vector()
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double rsq,eng,fpair,factor_coul,factor_lj,factor;
  double fdotv,massone,ke,pe;
  int *ilist,*jlist,*numneigh,**firstneigh;

  invoked_vector = update->ntimestep;

  double **x = atom->x; 
  double **v = atom->v; 
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  // invoke half neighbor list (will copy or build if necessary)

  neighbor->build_one(list->index);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // heat flux J = \sum_i e_i v_i + \sum_{i<j} (f_ij . v_j) x_ij

  // virial-like contribution
  // loop over neighbors of my atoms
  // require either i or j be in compute group

  double Jv[3] = {0.0,0.0,0.0};

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];

      if (j < nall) factor_coul = factor_lj = 1.0;
      else {
	factor_coul = special_coul[j/nall];
	factor_lj   = special_lj[j/nall];
	j %= nall;
      }

      if (!(mask[i] & groupbit) && !(mask[j] & groupbit)) continue;
 
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
	eng = pair->single(i,j,itype,jtype,rsq,factor_coul,factor_lj,fpair);
	
        if (newton_pair || j < nlocal) factor = 1.0;
        else factor = 0.5; 

        // symmetrize velocities

        double vx = 0.5*(v[i][0]+v[j][0]);
        double vy = 0.5*(v[i][1]+v[j][1]);
        double vz = 0.5*(v[i][2]+v[j][2]);
        fdotv = factor * fpair * (delx*vx + dely*vy + delz*vz); 

	Jv[0] += fdotv*delx;
	Jv[1] += fdotv*dely;
	Jv[2] += fdotv*delz;
      }                  
    }                  
  }

  // energy convection contribution
  // uses per-atom potential energy

  if (!(atomPE->invoked_flag & INVOKED_PERATOM)) {
    atomPE->compute_peratom();
    atomPE->invoked_flag |= INVOKED_PERATOM;
  }

  double *mass = atom->mass;
  double *rmass = atom->rmass;
  double mvv2e = force->mvv2e;

  double Jc[3] = {0.0,0.0,0.0};

  for (int i = 0; i < nlocal; i++) { 
    if (mask[i] & groupbit) {
      massone =  (rmass) ? rmass[i] : mass[type[i]];
      ke = mvv2e * 0.5 * massone *
	(v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
      pe = atomPE->scalar_atom[i]; 
      eng = pe + ke;
      Jc[0] += v[i][0]*eng; 
      Jc[1] += v[i][1]*eng;
      Jc[2] += v[i][2]*eng;
    }
  } 

  // total flux

  double data[6] = {Jv[0]+Jc[0],Jv[1]+Jc[1],Jv[2]+Jc[2],
		    Jc[0],Jc[1],Jc[2]};
  MPI_Allreduce(data,vector,6,MPI_DOUBLE,MPI_SUM,world);
}
