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

#include "pair_grid.h"
#include "pair_sna_grid.h"
#include "sna.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "pair.h"
#include "domain.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "tokenizer.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

PairSNAGrid::PairSNAGrid(LAMMPS *lmp) :
  PairGrid(lmp),
  radelem(nullptr), wjelem(nullptr)
{
  snaptr = nullptr;
}

/* ---------------------------------------------------------------------- */

PairSNAGrid::~PairSNAGrid()
{
  memory->destroy(radelem);
  memory->destroy(wjelem);
  memory->destroy(cutsq);
  delete snaptr;

  if (chemflag) memory->destroy(map);
}

/* ---------------------------------------------------------------------- */

void PairSNAGrid::init_style()
{
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style sna/grid requires newton pair on");

  // // need a full neighbor list

  // int irequest = neighbor->request(this,instance_me);
  // neighbor->requests[irequest]->half = 0;
  // neighbor->requests[irequest]->full = 1;

  snaptr = new SNA(lmp, rfac0, twojmax,
                   rmin0, switchflag, bzeroflag,
                   chemflag, bnormflag, wselfallflag,
		   nelements, switchinnerflag);  
  ncoeff = snaptr->ncoeff;
  ndesc = ndesc_base + ncoeff;
  snaptr->init();

}

/* ---------------------------------------------------------------------- */

void PairSNAGrid::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void PairSNAGrid::compute(int eflag, int vflag)
{
  double fij[3];
  
  ev_init(eflag,vflag);

  // compute sna for each gridpoint

  double** const x = atom->x;
  double **f = atom->f;
  const int* const mask = atom->mask;
  int * const type = atom->type;
  const int ntotal = atom->nlocal + atom->nghost;

  // insure rij, inside, and typej are of size ntotal
  
  snaptr->grow_rij(ntotal);

  // first generate fingerprint,
  // which allows calculation of beta
  
  for (int iz = nzlo; iz <= nzhi; iz++)
    for (int iy = nylo; iy <= nyhi; iy++)
      for (int ix = nxlo; ix <= nxhi; ix++) {
	double xgrid[3];
	grid2x(ix, iy, iz, xgrid);
	const double xtmp = xgrid[0];
	const double ytmp = xgrid[1];
	const double ztmp = xgrid[2];

	// currently, all grid points are type 1
	
	const int itype = 1;
	int ielem = 0;
	if (chemflag)
	  ielem = map[itype];
	const double radi = radelem[itype];

	// rij[][3] = displacements between atom I and those neighbors
	// inside = indices of neighbors of I within cutoff
	// typej = types of neighbors of I within cutoff

	int ninside = 0;
	for (int j = 0; j < ntotal; j++) {

	  const double delx = xtmp - x[j][0];
	  const double dely = ytmp - x[j][1];
	  const double delz = ztmp - x[j][2];
	  const double rsq = delx*delx + dely*dely + delz*delz;
	  int jtype = type[j];
	  int jelem = 0;
	  if (chemflag)
	    jelem = map[jtype];
	  if (rsq < cutsq[jtype][jtype] && rsq > 1e-20) {
	    snaptr->rij[ninside][0] = delx;
	    snaptr->rij[ninside][1] = dely;
	    snaptr->rij[ninside][2] = delz;
	    snaptr->inside[ninside] = j;
	    snaptr->wj[ninside] = wjelem[jtype];
	    snaptr->rcutij[ninside] = 2.0*radelem[jtype]*rcutfac;
	    snaptr->element[ninside] = jelem; // element index for multi-element snap
	    ninside++;
	  }
	}

	snaptr->compute_ui(ninside, ielem);
	snaptr->compute_zi();
	snaptr->compute_bi(ielem);

	// linear contributions

	for (int icoeff = 0; icoeff < ncoeff; icoeff++)
	  gridlocal[ndesc_base+icoeff][iz][iy][ix] = 
	    snaptr->blist[icoeff];

	// quadratic contributions
	// untested

	if (quadraticflag) {
	  int ncount = ncoeff;
	  for (int icoeff = 0; icoeff < ncoeff; icoeff++) {
	    double bveci = snaptr->blist[icoeff];
	    gridlocal[ndesc_base+ncount++][iz][iy][ix] = 
	      0.5*bveci*bveci;
	    for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++)
	      gridlocal[ndesc_base+ncount++][iz][iy][ix] = 
		bveci*snaptr->blist[jcoeff];
	  }
	}
      }

  // this is a proxy for a call to the energy model
  // beta is dE/dB^i, the derivative of the total
  // energy w.r.t. to descriptors of grid point i
  
  compute_beta();
  
  // second compute forces using beta
  
  int igrid = 0;
  for (int iz = nzlo; iz <= nzhi; iz++)
    for (int iy = nylo; iy <= nyhi; iy++)
      for (int ix = nxlo; ix <= nxhi; ix++) {
	double xgrid[3];
	grid2x(ix, iy, iz, xgrid);
	const double xtmp = xgrid[0];
	const double ytmp = xgrid[1];
	const double ztmp = xgrid[2];

	// currently, all grid points are type 1
	
	const int itype = 1;
	int ielem = 0;
	if (chemflag)
	  ielem = map[itype];
	const double radi = radelem[itype];

	// rij[][3] = displacements between atom I and those neighbors
	// inside = indices of neighbors of I within cutoff
	// typej = types of neighbors of I within cutoff

	int ninside = 0;
	for (int j = 0; j < ntotal; j++) {

	  const double delx = xtmp - x[j][0];
	  const double dely = ytmp - x[j][1];
	  const double delz = ztmp - x[j][2];
	  const double rsq = delx*delx + dely*dely + delz*delz;
	  int jtype = type[j];
	  int jelem = 0;
	  jelem = map[jtype];

	  if (rsq < cutsq[jtype][jtype] && rsq > 1e-20) {
	    snaptr->rij[ninside][0] = delx;
	    snaptr->rij[ninside][1] = dely;
	    snaptr->rij[ninside][2] = delz;
	    snaptr->inside[ninside] = j;
	    snaptr->wj[ninside] = wjelem[jtype];
	    snaptr->rcutij[ninside] = 2.0*radelem[jtype]*rcutfac;
	    if (switchinnerflag) {
	      snaptr->sinnerij[ninside] = 0.5*(sinnerelem[ielem]+sinnerelem[jelem]);
	      snaptr->dinnerij[ninside] = 0.5*(dinnerelem[ielem]+dinnerelem[jelem]);
	    }
	    if (chemflag) snaptr->element[ninside] = jelem;
	    ninside++;
	  }
	}

	// compute Ui, Yi for atom I
	
	if (chemflag)
	  snaptr->compute_ui(ninside, ielem);
	else
	  snaptr->compute_ui(ninside, 0);

	// for neighbors of I within cutoff:
	// compute Fij = dEi/dRj = -dEi/dRi
	// add to Fi, subtract from Fj
	// scaling is that for type I

	snaptr->compute_yi(beta[igrid]);

	for (int jj = 0; jj < ninside; jj++) {
	  int j = snaptr->inside[jj];
	  snaptr->compute_duidrj(jj);

	  snaptr->compute_deidrj(fij);

	  f[j][0] += fij[0];
	  f[j][1] += fij[1];
	  f[j][2] += fij[2];

	  // tally per-atom virial contribution

	  if (vflag)
	    ev_tally_xyz(-1,j,atom->nlocal,force->newton_pair,0.0,0.0,
			 fij[0],fij[1],fij[2],
			 -snaptr->rij[jj][0],-snaptr->rij[jj][1],
			 -snaptr->rij[jj][2]);
	}

	// tally energy contribution

	if (eflag) {

	  // get descriptors again

	  snaptr->compute_zi();
	  snaptr->compute_bi(ielem);
	  
	  // evdwl = energy of atom I, sum over coeffs_k * Bi_k

	  double evdwl = 0.0;

	  // E = beta.B 

	  for (int icoeff = 0; icoeff < ncoeff; icoeff++)
	    evdwl += beta[igrid][icoeff]*snaptr->blist[icoeff];
	  
	  ev_tally_full(-1,2.0*evdwl,0.0,0.0,0.0,0.0,0.0);

	}
	igrid++;
      }

  if (vflag_fdotr) virial_fdotr_compute();
  
}


/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairSNAGrid::settings(int narg, char ** arg)
{

  // call base class first

  PairGrid::settings(narg, arg);
    
  // skip over arguments used by base class
  // so that argument positions are identical to
  // regular per-atom compute
  
  arg += nargbase;
  narg -= nargbase;

  int ntypes = atom->ntypes;
  int nargmin = 3+2*ntypes;

  if (narg < nargmin) error->all(FLERR,"Illegal pair sna/grid command");

  // default values

  rmin0 = 0.0;
  switchflag = 1;
  bzeroflag = 1;
  quadraticflag = 0;
  chemflag = 0;
  bnormflag = 0;
  wselfallflag = 0;
  switchinnerflag = 0;
  nelements = 1;
  
  // process required arguments

  memory->create(radelem,ntypes+1,"pair:sna/grid:radelem"); // offset by 1 to match up with types
  memory->create(wjelem,ntypes+1,"pair:sna/grid:wjelem");

  rcutfac = atof(arg[0]);
  rfac0 = atof(arg[1]);
  twojmax = atoi(arg[2]);

  for(int i = 0; i < ntypes; i++)
    radelem[i+1] = atof(arg[3+i]);
  for(int i = 0; i < ntypes; i++)
    wjelem[i+1] = atof(arg[3+ntypes+i]);

  // construct cutsq

  double cut;
  cutmax = 0.0;
  memory->create(cutsq,ntypes+1,ntypes+1,"pair:sna/grid:cutsq");
  for(int i = 1; i <= ntypes; i++) {
    cut = 2.0*radelem[i]*rcutfac;
    if (cut > cutmax) cutmax = cut;
    cutsq[i][i] = cut*cut;
    for(int j = i+1; j <= ntypes; j++) {
      cut = (radelem[i]+radelem[j])*rcutfac;
      cutsq[i][j] = cutsq[j][i] = cut*cut;
    }
  }

  // set local input checks

  int sinnerflag = 0;
  int dinnerflag = 0;

  // process optional args

  int iarg = nargmin;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"rmin0") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal pair sna/grid command");
      rmin0 = atof(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"switchflag") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal pair sna/grid command");
      switchflag = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"bzeroflag") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal pair sna/grid command");
      bzeroflag = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"quadraticflag") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal pair sna/grid command");
      quadraticflag = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"chem") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal pair sna/grid command");
      chemflag = 1;
      memory->create(map,ntypes+1,"pair:sna/grid:map");
      nelements = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      for (int i = 0; i < ntypes; i++) {
        int jelem = utils::inumeric(FLERR,arg[iarg+2+i],false,lmp);
        if (jelem < 0 || jelem >= nelements)
          error->all(FLERR,"Illegal pair sna/grid command");
        map[i+1] = jelem;
      }
      iarg += 2+ntypes;
    } else if (strcmp(arg[iarg],"bnormflag") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal pair sna/grid command");
      bnormflag = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"wselfallflag") == 0) {
      if (iarg+2 > narg)
        error->all(FLERR,"Illegal pair sna/grid command");
      wselfallflag = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"switchinnerflag") == 0) {
      if (iarg+2 > narg)
	error->all(FLERR,"Illegal pair sna/grid command");
      switchinnerflag = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"sinner") == 0) {
      iarg++;
      if (iarg+ntypes > narg)
	error->all(FLERR,"Illegal pair sna/grid command");
      memory->create(sinnerelem,ntypes+1,"snap:sinnerelem");
      for (int i = 0; i < ntypes; i++)
        sinnerelem[i+1] = utils::numeric(FLERR,arg[iarg+i],false,lmp);
      sinnerflag = 1;
      iarg += ntypes;
    } else if (strcmp(arg[iarg],"dinner") == 0) {
      iarg++;
      if (iarg+ntypes > narg)
	error->all(FLERR,"Illegal pair sna/grid command");
      memory->create(dinnerelem,ntypes+1,"snap:dinnerelem");
      for (int i = 0; i < ntypes; i++)
        dinnerelem[i+1] = utils::numeric(FLERR,arg[iarg+i],false,lmp);
      dinnerflag = 1;
      iarg += ntypes;
    } else error->all(FLERR,"Illegal pair sna/grid command");

  }

}
  
/* ----------------------------------------------------------------------
   memory usage
------------------------------------------------------------------------- */

double PairSNAGrid::memory_usage()
{
  double nbytes = snaptr->memory_usage();    // SNA object
  int n = atom->ntypes+1;
  nbytes += (double)n*sizeof(int);            // map

  return nbytes;
}

